"""SAC preparation — attribution and provenance.

This module was compiled with Claude (Anthropic) assistance. The processing
chain follows the mttime example notebooks by Andrea Chiang (LLNL),
specifically 01_Data_Processing and 02_Prepare_Data_and_Synthetics_For_
Inversion: https://github.com/LLNL/mttime/tree/master/examples/notebooks
(mttime: https://github.com/LLNL/mttime, LLNL-CODE-814839), implemented
with ObsPy (https://github.com/obspy/obspy).

Deviations from the original notebooks: GeoNet NRT/archive FDSN sources
with retry logic; automated station pooling (magnitude-scaled distance
window with one-shot radius extension, HH?>BH? priority); peak-to-noise
quality measured in the distance-adaptive inversion window;
magnitude-dependent filter-band menu applied per event rather than fixed
corners; fail-loud drop accounting (every rejected station recorded with
a reason).

Selection v4 (2026-09-04): this module DELETES ONLY UNUSABLE DATA — no
waveform, no response, gaps, a dead channel, or a broken-response
amplitude outlier. Everything else enters the pool carrying demotion
TAGS (near_field / weak_signal / cluster_surplus) and earns or loses its
seat by FIT, in the funnel in invert.py. The previous scheme's
pre-filters were removing good data: the 3x-depth rule alone took the
closest station out of every event it touched.

Waveform acquisition and pre-processing to mttime-ready SAC files.

Chain (mttime example notebooks 01+02, matching the EPS207 recipe):
  inventory -> broadband selection -> download -> response removal to
  displacement -> rotate ZNE -> NE->RT -> bandpass -> decimate to 1 sps ->
  trim origin-30s..origin+200s -> m to cm -> SAC "NET.STA.LOC.{Z,R,T}.dat"

Every dropped station/trace is recorded with a reason; nothing is padded or
substituted silently.
"""

from __future__ import annotations

from pathlib import Path

from obspy import Stream, UTCDateTime, read
from obspy.core.util.attribdict import AttribDict
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees

import config
from geonet import Event, fdsn_client


def _cache_dir(public_id: str) -> Path:
    """Per-event raw-download cache (see config.WF_CACHE_DIR: disposable)."""
    d = config.WF_CACHE_DIR / public_id
    d.mkdir(parents=True, exist_ok=True)
    readme = config.WF_CACHE_DIR / "README.txt"
    if not readme.exists():
        readme.write_text(
            "auto_tdmt_NZ raw-waveform download cache.\n"
            "Everything here is re-downloadable from GeoNet FDSN and is\n"
            "DISPOSABLE: delete this directory at any time (rm -rf) —\n"
            "nothing else references it; the next run just re-downloads.\n")
    return d


def select_stations(client, event: Event, origin: UTCDateTime,
                    max_dist_km: float | None = None):
    """Broadband NZ stations within the working distance range.

    Returns (inventory, station_rows) where station_rows is a list of dicts
    with one preferred (channel-band, location) per station, nearest first.
    ``max_dist_km`` overrides the magnitude-scaled radius (radius extension).
    """
    max_dist = max_dist_km or config.station_max_dist_km(event.prelim_mag)
    min_dist = config.station_min_dist_km(event.prelim_mag)
    inv_cache = (_cache_dir(event.public_id)
                 / f"inventory_{int(round(max_dist))}km.xml")
    if inv_cache.exists():
        from obspy import read_inventory
        inv = read_inventory(str(inv_cache))
    else:
        inv = client.get_stations(
            network=config.NETWORK,
            channel=",".join(config.CHANNEL_PRIORITY),
            latitude=event.latitude,
            longitude=event.longitude,
            maxradius=kilometers2degrees(max_dist),
            level="response",
            starttime=origin,
            endtime=origin + config.TIME_AFTER_S,
        )
        inv.write(str(inv_cache), format="STATIONXML")
    rows = []
    for net in inv:
        for sta in net:
            # channels present for this station, grouped by (loc, band)
            groups = {}
            for cha in sta:
                band = cha.code[:2]
                groups.setdefault((cha.location_code, band), set()).add(cha.code)
            # prefer HH over BH; need all three components in one group
            chosen = None
            for prio in config.CHANNEL_PRIORITY:
                band = prio[:2]
                for (loc, b), codes in sorted(groups.items()):
                    if b == band and len(codes) >= 3:
                        chosen = (loc, band)
                        break
                if chosen:
                    break
            if chosen is None:
                continue
            dist_m, az, baz = gps2dist_azimuth(
                event.latitude, event.longitude, sta.latitude, sta.longitude
            )
            dist_km = dist_m / 1000.0
            if not (min_dist <= dist_km <= max_dist):
                continue
            rows.append(
                dict(
                    network=net.code,
                    station=sta.code,
                    location=chosen[0],
                    band=chosen[1],
                    latitude=sta.latitude,
                    longitude=sta.longitude,
                    distance_km=dist_km,
                    azimuth=az,
                    back_azimuth=baz,
                )
            )
    rows.sort(key=lambda r: r["distance_km"])
    assert rows, "no broadband stations found in distance range"
    return inv, rows


def demotion_tags(row: dict, event: Event) -> list[str]:
    """Demotion tags for a pool station — reasons to be SUSPICIOUS, never
    reasons to exclude. The funnel tests every tagged station anyway; the
    tags only shift the burden of proof (and appear on the figures).

    - near_field: inside 3x the source depth, where the point-source
      far-field heuristic is uncomfortable. The CPS Green's functions are
      complete-wavefield (near-field terms included), so the fit decides.
    - weak_signal: peak/noise below PEAK_NOISE_STRONG.
    """
    tags = []
    if (event.depth_km not in config.PLACEHOLDER_DEPTHS_KM
            and event.depth_km <= config.DIST_DEPTH_RULE_MAX_DEPTH_KM
            and row["distance_km"] < config.MIN_DIST_DEPTH_RATIO
            * event.depth_km):
        tags.append("near_field")
    if row.get("pk_n", 0.0) < config.PEAK_NOISE_STRONG:
        tags.append("weak_signal")
    return tags


def fetch_and_process(
    event: Event, workdir: Path, band_hz: tuple[float, float],
    stages: dict | None = None,
) -> tuple[list[dict], list[dict]]:
    """Download and pre-process waveforms for one event.

    Writes SAC files ``<workdir>/NET.STA.LOC.{Z,R,T}.dat`` and returns
    (pool_rows, dropped). The pool is every station with usable data, each
    row carrying pk_n, window_end_s, sector, tier ("trusted"/"demoted")
    and tags; ``dropped`` holds only unusable data, each with a reason
    string in the shared vocabulary (see invert.reason_class).

    If ``stages`` (a dict) is passed, copies of every pool station's stream
    are stored under keys "raw", "displacement", "final" for QC plotting.
    """
    origin = UTCDateTime(event.origin_time)
    client = fdsn_client(origin)
    inv, rows = select_stations(client, event, origin)

    workdir.mkdir(parents=True, exist_ok=True)
    fmin, fmax = band_hz
    if stages is not None:
        stages.update({"raw": Stream(), "displacement": Stream(), "final": Stream()})

    # azimuth-sector-interleaved candidate order (8 x 45 deg sectors,
    # nearest-first within each): the first stations tried span the full
    # compass, codifying quadrant-first manual selection. Pure
    # nearest-first ordering produced one-sided geometries when the close
    # stations clustered on one side of the epicentre.
    sectors: dict[int, list] = {}
    for r in rows:
        sectors.setdefault(int(r["azimuth"] // 45) % 8, []).append(r)
    ordered = []
    while any(sectors.values()):
        for k in sorted(sectors):
            if sectors[k]:
                ordered.append(sectors[k].pop(0))

    import numpy as np

    used, dropped = [], []
    queue = list(ordered)
    extended = False
    i = 0
    while True:
        if i >= len(queue):
            # single radius extension when the usable pool is starved
            # (offshore events: 2026p047833 review)
            base = config.station_max_dist_km(event.prelim_mag)
            if extended or len(used) >= config.MIN_USABLE_BEFORE_EXTEND:
                break
            extended = True
            ext = base + config.RADIUS_EXTEND_KM
            try:
                inv_ext, more = select_stations(
                    client, event, origin, max_dist_km=ext)
            except AssertionError:
                break
            # the annulus stations' responses live in the EXTENDED
            # inventory — merge it, or response removal fails for them
            inv += inv_ext
            annulus = [r for r in more
                       if r["distance_km"] > base + 1e-6]
            if not annulus:
                break
            print(f"only {len(used)} usable stations: radius extended to "
                  f"{ext:g} km ({len(annulus)} more candidates)")
            queue.extend(sorted(annulus, key=lambda r: r["distance_km"]))
            continue
        row = queue[i]
        i += 1
        if len(used) >= config.MAX_POOL_STATIONS:
            break
        sid = f"{row['network']}.{row['station']}.{row['location']}"

        def _drop(reason, stage="data"):
            """Record an UNUSABLE station (the only kind removed here)."""
            dropped.append({
                **row, "station": sid, "reason": reason, "stage": stage,
                "status": "rejected",
                "distance_km": round(row["distance_km"], 1),
            })
        # NOTE (v4): the 3x-depth far-field rule used to drop stations here.
        # It now only TAGS them (demotion_tags) — the audit found it
        # removing the closest station in every event it touched, and 94%
        # of those events ended grade C/D.
        wf_cache = (_cache_dir(event.public_id)
                    / f"{sid}.{row['band']}.mseed")
        try:
            if wf_cache.exists():
                st = read(str(wf_cache))
            else:
                st = client.get_waveforms(
                    network=row["network"],
                    station=row["station"],
                    location=row["location"],
                    channel=f"{row['band']}?",
                    starttime=origin - 5 * config.TIME_BEFORE_S,
                    endtime=origin + config.TIME_AFTER_S
                    + config.TIME_BEFORE_S,
                    attach_response=False,
                )
                st.write(str(wf_cache), format="MSEED")
        except Exception as e:  # noqa: BLE001 - record and move on, loudly
            _drop(f"no data: download failed: {e}")
            continue

        st.merge(method=0)
        if any(hasattr(tr.data, "mask") for tr in st) or len(st) < 3:
            _drop("no data: gaps or <3 components")
            continue

        raw_copy = st.copy() if stages is not None else None

        try:
            st.detrend("linear")
            st.remove_response(
                inventory=inv,
                pre_filt=config.RESPONSE_PRE_FILT,
                output="DISP",
                zero_mean=True,
            )
            st.detrend("linear")
            st.detrend("demean")
            st._rotate_to_zne(inv, components=("ZNE", "Z12"))
        except Exception as e:  # noqa: BLE001
            _drop(f"no data: response/rotation failed: {e}")
            continue

        if len(st.select(component="Z")) != 1 or len(st) != 3:
            _drop("no data: not exactly 3 ZNE components")
            continue

        if stages is not None:
            for src, key in ((raw_copy, "raw"), (st.copy(), "displacement")):
                for tr in src:
                    tr.stats.distance = row["distance_km"] * 1000.0
                stages[key] += src

        for tr in st:
            tr.stats.back_azimuth = row["back_azimuth"]
        st.rotate(method="NE->RT")

        st.filter(
            "bandpass", freqmin=fmin, freqmax=fmax,
            corners=config.FILTER_CORNERS, zerophase=True,
        )
        st.taper(max_percentage=0.05)
        for tr in st:
            factor = int(round(tr.stats.sampling_rate * config.DT))
            assert factor >= 1, f"{sid}: sampling rate {tr.stats.sampling_rate}"
            tr.decimate(factor=factor, strict_length=False, no_filter=True)
            tr.resample(1.0 / config.DT, strict_length=False, no_filter=True)
            tr.trim(
                origin - config.TIME_BEFORE_S,
                origin + config.TIME_AFTER_S,
                nearest_sample=True,
            )
            tr.data = 100.0 * tr.data  # m -> cm (TDMT convention)

        npts = {tr.stats.npts for tr in st}
        expected = config.TIME_BEFORE_S + config.TIME_AFTER_S + 1
        if npts != {expected}:
            _drop(f"no data: trim gave npts {npts}, want {expected}")
            continue

        # peak-to-noise, measured ONLY inside the distance-adaptive window
        # actually inverted (see config PEAK_NOISE_* docs). Dead channels
        # are still written under a rejected_ prefix so the per-band
        # all-station waveform figure can show WHY they were excluded;
        # mtinv.in never references prefixed files.
        tail = config.window_tail_s(event.prelim_mag)
        wlen = int(min(config.INV_NPTS, max(
            config.WINDOW_MIN_S,
            config.TIME_BEFORE_S
            + row["distance_km"] / config.WINDOW_GROUP_VEL_KMS + tail)))
        tend = wlen - config.TIME_BEFORE_S  # s after origin (kinematic MIN)

        # signal-aware window end (REVIEW_LEARNINGS item 13): slow paths
        # (Hikurangi accretionary prism, ~1.2-1.7 km/s effective group
        # velocity) deliver trains the kinematic cut bisects. Extend the
        # end to where the smoothed 3-component envelope decays back
        # toward the pre-event level — deterministic, data-derived —
        # never shorter than the kinematic value, capped by INV_NPTS.
        ref = st[0]
        b0 = -1.0 * (origin - ref.stats.starttime)
        tt = b0 + np.arange(ref.stats.npts) * ref.stats.delta
        env = np.max([np.abs(tr.data) for tr in st], axis=0)
        win = int(round(15.0 / config.DT))  # 15 s smoothing
        env = np.convolve(env, np.ones(win) / win, mode="same")
        noise_lvl = float(np.median(env[(tt > b0 + 5) & (tt < -2)]))
        thresh = max(2.0 * noise_lvl, 0.05 * float(env.max()))
        t_cap = config.INV_NPTS - config.TIME_BEFORE_S  # 120 s after origin
        live = (tt >= tend) & (tt <= t_cap) & (env > thresh)
        if live.any():
            tend = int(min(t_cap, tt[live].max() + 10.0))  # +10 s pad
        ratios = []
        for tr in st:
            b = -1.0 * (origin - tr.stats.starttime)
            t = b + np.arange(tr.stats.npts) * tr.stats.delta
            noise = tr.data[(t > b + 5) & (t < -2)]
            sig = tr.data[(t >= 0) & (t <= tend)]
            assert len(noise) > 10 and len(sig) > 10, \
                f"{sid}: peak/noise windows too short"
            ratios.append(
                float(np.abs(sig).max() / np.sqrt(np.mean(noise ** 2))))
        pk_n = float(np.median(ratios))
        row["pk_n"] = round(pk_n, 2)
        row["window_end_s"] = tend
        row["sector"] = config.sector(row["azimuth"])
        row["tags"] = demotion_tags(row, event)
        row["tier"] = "trusted" if not row["tags"] else "demoted"
        row["status"] = "pool"
        prefix = "rejected_" if pk_n < config.PEAK_NOISE_DEAD else ""

        for tr in st:
            sacd = AttribDict()
            sacd.stla, sacd.stlo = row["latitude"], row["longitude"]
            sacd.evla, sacd.evlo = event.latitude, event.longitude
            sacd.evdp = event.depth_km * 1000.0
            sacd.dist = row["distance_km"]
            sacd.az, sacd.baz = row["azimuth"], row["back_azimuth"]
            sacd.o = 0.0
            sacd.b = -1.0 * (origin - tr.stats.starttime)
            tr.stats.sac = sacd
            comp = tr.stats.channel[-1]
            assert comp in "ZRT", f"{sid}: unexpected component {comp}"
            tr.write(str(workdir / f"{prefix}{sid}.{comp}.dat"),
                     format="SAC")

        if prefix:
            _drop(f"dead channel: peak/noise {pk_n:.1f} < "
                  f"{config.PEAK_NOISE_DEAD:g}")
            continue

        if stages is not None:
            final_copy = st.copy()
            for tr in final_copy:
                tr.stats.distance = row["distance_km"] * 1000.0
            stages["final"] += final_copy

        row["filter_hz"] = [fmin, fmax]
        used.append(row)

    # amplitude-consistency screen: peak x distance should be comparable
    # across the network; a station orders of magnitude off has broken
    # response metadata (e.g. NZ.RDHZ 2026-09) and would single-handedly
    # steer the least-squares moment. Screen BEFORE inversion.
    import numpy as np

    if len(used) >= config.AMP_SCREEN_MIN_STATIONS:
        for r in used:
            tr = read(str(workdir / f"{r['network']}.{r['station']}."
                          f"{r['location']}.Z.dat"), format="SAC")[0]
            r["peak_x_dist"] = float(
                abs(tr.data).max() * r["distance_km"])
        med = float(np.median([r["peak_x_dist"] for r in used]))
        for r in used:
            r["amp_ratio"] = round(r["peak_x_dist"] / med, 3)
        # ONE-SIDED (2026-09-04, D. Lindsay): only a station far ABOVE the
        # network median is evidence of a broken response. A station far
        # BELOW it may simply sit near a nodal plane — "just because it
        # doesn't see signal doesn't mean it's a bad fit" — and its small
        # amplitude is real information about the mechanism. Dead channels
        # are caught by the peak/noise floor instead.
        flagged = [r for r in used
                   if r["peak_x_dist"]
                   > med * config.AMPLITUDE_OUTLIER_FACTOR]
        for r in flagged:
            sid = f"{r['network']}.{r['station']}.{r['location']}"
            for comp in "ZRT":
                p = workdir / f"{sid}.{comp}.dat"
                if p.exists():  # keep for the all-station waveform figure
                    p.rename(workdir / f"rejected_{p.name}")
            dropped.append({
                **r, "station": sid,
                "reason": f"amplitude outlier: peak x dist "
                          f"{r['peak_x_dist']:.2e} vs network median "
                          f"{med:.2e}",
                "stage": "data", "status": "rejected",
                "distance_km": round(r["distance_km"], 1),
            })
            print(f"  amplitude outlier screened: {sid} "
                  f"({r['peak_x_dist']/med:.0f}x median)")
        used = [r for r in used if r not in flagged]

    # station-cluster TAGGING (v4: was thinning): dense sub-networks (the
    # Ruapehu ring is 69% of all historic cluster drops) must not stack
    # near-identical records into one azimuth sector — but the decision
    # belongs to the fit, so beyond CLUSTER_MAX_STATIONS within
    # CLUSTER_RADIUS_KM a station is only tagged "cluster_surplus" (and
    # must show a joint-fit gain to earn its seat in the funnel).
    anchors: list[dict] = []
    for r in sorted(used, key=lambda x: -x.get("pk_n", 0.0)):
        near = [k for k in anchors if gps2dist_azimuth(
            r["latitude"], r["longitude"],
            k["latitude"], k["longitude"])[0] / 1000.0
            <= config.CLUSTER_RADIUS_KM]
        if len(near) >= config.CLUSTER_MAX_STATIONS:
            r.setdefault("tags", []).append("cluster_surplus")
            r["cluster_with"] = [k["station"] for k in near]
            r["tier"] = "demoted"
        else:
            anchors.append(r)

    pool = sorted(used, key=lambda r: r["distance_km"])
    # NOTE: no assert — an empty pool is a legitimate outcome (offshore
    # events with no usable data). run02 turns it into a "no coherent
    # solution" record rather than a hard failure.
    return pool, dropped
