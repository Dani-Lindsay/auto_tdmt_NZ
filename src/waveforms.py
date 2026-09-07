"""Task 2 — station selection and waveform processing.

    pixi run python waveforms.py --event 2026p666955      # run this task alone

Chain (mttime example notebooks 01+02 by Andrea Chiang, LLNL,
https://github.com/LLNL/mttime/tree/master/examples/notebooks, with
ObsPy): inventory -> broadband selection -> download -> response removal
to displacement -> rotate to ZRT -> bandpass -> 1 sps -> trim
origin-30 s .. origin+200 s -> cm -> SAC "NET.STA.LOC.{Z,R,T}.dat".

Selection rules, each with its source (details in auto_tdmt.cfg §2):
  1. distance window scaled by magnitude            (Gisola; Triantafyllis et al. 2022)
  2. unusable data rejected: no waveform/response, gaps, wrong sample count
  3. per-component SNR = RMS(signal) / RMS(noise) over the inverted
     window; station usable when its BEST component >= snrMin
                                                     (BMKG; Halauwet et al. 2024)
  4. amplitude outliers above amplitudeMaxRatio x the network median
     rejected (broken response); the LOW side is kept - a nodal station
     genuinely has small amplitude                   (Duputel et al. 2012)
  5. candidates ranked by SNR and taken ROUND-ROBIN across azimuth
     sectors until maxStations, so coverage is balanced when the geometry
     allows it and never starved when it does not    (Gisola; SCARDEC; INGV)
The fit then decides (task 3). Nothing else is filtered here.

Every station that enters and leaves is recorded with a reason string in
the shared vocabulary (invert.reason_class) so the figures and the
station ledger can show why.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from obspy import Stream, UTCDateTime, read
from obspy.core.util.attribdict import AttribDict
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees

import config
from config import P
from geonet import Event, fdsn_client


def _cache_dir(public_id: str) -> Path:
    """Per-event raw-download cache (config.WF_CACHE_DIR: disposable)."""
    d = config.WF_CACHE_DIR / public_id
    d.mkdir(parents=True, exist_ok=True)
    readme = config.WF_CACHE_DIR / "README.txt"
    if not readme.exists():
        readme.write_text(
            "auto_tdmt_NZ raw-waveform download cache.\n"
            "Everything here is re-downloadable from GeoNet FDSN and is\n"
            "DISPOSABLE: delete this directory at any time (rm -rf).\n")
    return d


def station_id(row: dict) -> str:
    return f"{row['network']}.{row['station']}.{row['location']}"


# ---------------------------------------------------------------------------
# 2a. candidate inventory
# ---------------------------------------------------------------------------

def select_stations(client, event: Event, origin: UTCDateTime,
                    max_dist_km: float | None = None):
    """Broadband NZ stations within the working distance range.

    Returns (inventory, rows): one preferred (channel band, location) per
    station, nearest first. ``max_dist_km`` overrides the magnitude-scaled
    radius (used by the one-shot radius extension).
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
            groups: dict = {}
            for cha in sta:
                groups.setdefault((cha.location_code, cha.code[:2]),
                                  set()).add(cha.code)
            chosen = None
            for prio in config.CHANNEL_PRIORITY:   # HH over BH
                for (loc, b), codes in sorted(groups.items()):
                    if b == prio[:2] and len(codes) >= 3:
                        chosen = (loc, b)
                        break
                if chosen:
                    break
            if chosen is None:
                continue
            dist_m, az, baz = gps2dist_azimuth(
                event.latitude, event.longitude, sta.latitude, sta.longitude)
            dist_km = dist_m / 1000.0
            if not (min_dist <= dist_km <= max_dist):
                continue
            rows.append(dict(
                network=net.code, station=sta.code, location=chosen[0],
                band=chosen[1], latitude=sta.latitude, longitude=sta.longitude,
                distance_km=dist_km, azimuth=az, back_azimuth=baz,
                sector=config.sector(az),
            ))
    rows.sort(key=lambda r: r["distance_km"])
    assert rows, "no broadband stations found in distance range"
    return inv, rows


# ---------------------------------------------------------------------------
# 2b. signal quality
# ---------------------------------------------------------------------------

def component_snr(st: Stream, origin: UTCDateTime, distance_km: float,
                  signal_end_s: float, window_s: float | None = None
                  ) -> dict[str, float]:
    """Per-component SNR = RMS(signal window after P) / RMS(same length
    before P), on the filtered traces (Halauwet et al. 2024, GJI 239,
    BMKG). The signal window is the one actually inverted — P to
    ``signal_end_s`` after origin — capped at snrWindowS: BMKG's fixed
    200 s suits their M >= 5 events, but on an M3.7 the surface-wave
    train is over in 60 s and a 200 s window dilutes it to noise
    (2026p669681: every station scored ~1.2 that way). P is estimated
    kinematically at 6 km/s; the exact pick does not matter here."""
    t_p = origin + distance_km / 6.0
    window_s = min(window_s or P.station.snrWindowS,
                   max(20.0, signal_end_s - distance_km / 6.0))
    out = {}
    for tr in st:
        comp = tr.stats.channel[-1]
        t = tr.times(reftime=t_p)
        noise = tr.data[(t >= -window_s) & (t < 0)]
        signal = tr.data[(t >= 0) & (t < window_s)]
        assert len(noise) > 10 and len(signal) > 10, (
            f"{tr.id}: SNR windows too short ({len(noise)}/{len(signal)})")
        rms_n = float(np.sqrt(np.mean(noise ** 2)))
        rms_s = float(np.sqrt(np.mean(signal ** 2)))
        out[comp] = round(rms_s / rms_n, 2) if rms_n > 0 else float("inf")
    return out


def signal_window_end_s(st: Stream, origin: UTCDateTime, distance_km: float,
                        prelim_mag: float) -> int:
    """End of the inversion window, seconds after origin: the kinematic
    minimum (distance / groupVel + tail), extended to where the smoothed
    3-component envelope decays back toward the pre-event level (slow
    Hikurangi paths), never shorter than the kinematic value and capped
    at maxWindowS."""
    kin = max(P.station.windowMinS - config.TIME_BEFORE_S,
              distance_km / P.station.groupVelKms
              + config.window_tail_s(prelim_mag))
    # the extension exists for slow paths that add tens of seconds, not
    # for noise: it is bounded to windowExtendMaxS beyond the kinematic
    # end (2026p669681: a 23 km station was extended to 187 s of noise
    # and "aligned" at -21 s)
    cap = min(P.station.maxWindowS, kin + P.station.windowExtendMaxS)
    tt = st[0].times(reftime=origin)
    env = np.max([np.abs(tr.data) for tr in st], axis=0)
    win = max(1, int(round(15.0 / config.DT)))          # 15 s smoothing
    env = np.convolve(env, np.ones(win) / win, mode="same")
    pre = env[(tt > tt[0] + 5) & (tt < -2)]
    noise_lvl = float(np.median(pre)) if pre.size else 0.0
    thresh = max(3.0 * noise_lvl, 0.1 * float(env.max()))
    live = (tt >= kin) & (tt <= cap) & (env > thresh)
    tend = tt[live].max() + 10.0 if live.any() else kin    # +10 s pad
    return int(min(cap, max(kin, tend)))


# ---------------------------------------------------------------------------
# 2c. azimuth-balanced selection
# ---------------------------------------------------------------------------

def select_by_sector(rows: list[dict], max_stations: int | None = None,
                     n_sectors: int | None = None) -> tuple[list, list]:
    """Rank by SNR, then take stations round-robin across populated azimuth
    sectors until ``max_stations``. Balanced coverage when the geometry
    allows it; a one-sided (offshore) geometry still fills up from its
    best stations rather than being starved by a per-sector cap.
    Returns (selected, left_out)."""
    max_stations = max_stations or P.station.maxStations
    n_sectors = n_sectors or P.station.sectors
    by_sector: dict[int, list] = {}
    for r in sorted(rows, key=lambda r: -r.get("snr_best", r["snr_med"])):
        by_sector.setdefault(int(r["azimuth"] // (360.0 / n_sectors))
                             % n_sectors, []).append(r)
    selected: list[dict] = []
    while len(selected) < max_stations and any(by_sector.values()):
        for k in sorted(by_sector):
            if by_sector[k] and len(selected) < max_stations:
                selected.append(by_sector[k].pop(0))
    left = [r for r in rows if r not in selected]
    return selected, left


# ---------------------------------------------------------------------------
# 2d. the task
# ---------------------------------------------------------------------------

def fetch_and_process(
    event: Event, workdir: Path, band_hz: tuple[float, float],
    stages: dict | None = None,
) -> tuple[list[dict], list[dict]]:
    """Download and pre-process waveforms for one event and one band.

    Writes SAC files ``<workdir>/NET.STA.LOC.{Z,R,T}.dat`` for the selected
    stations (``rejected_`` prefix for the rest, so the all-station figure
    can show them) and returns (pool, dropped): the selected station rows
    (with snr, snr_med, window_end_s, sector, amp_ratio) and every other
    candidate with a reason.

    If ``stages`` (a dict) is passed, copies of every usable station's
    stream are stored under "raw", "displacement", "final" for QC plots.
    """
    origin = UTCDateTime(event.origin_time)
    client = fdsn_client(origin)
    inv, rows = select_stations(client, event, origin)
    workdir.mkdir(parents=True, exist_ok=True)
    fmin, fmax = band_hz
    if stages is not None:
        stages.update({"raw": Stream(), "displacement": Stream(),
                       "final": Stream()})
    pre_s = P.station.snrWindowS + config.TIME_BEFORE_S + 60.0

    usable, dropped = [], []
    queue = list(rows)
    extended = False
    i = 0
    while True:
        if i >= len(queue):
            # one-shot radius extension when the usable pool is thin
            base = config.station_max_dist_km(event.prelim_mag)
            if extended or len(usable) >= config.MIN_USABLE_BEFORE_EXTEND:
                break
            extended = True
            ext = base + config.RADIUS_EXTEND_KM
            try:
                inv_ext, more = select_stations(client, event, origin,
                                                max_dist_km=ext)
            except AssertionError:
                break
            inv += inv_ext   # responses for the annulus live here
            annulus = [r for r in more if r["distance_km"] > base + 1e-6]
            if not annulus:
                break
            print(f"only {len(usable)} usable stations: radius extended to "
                  f"{ext:g} km ({len(annulus)} more candidates)")
            queue.extend(annulus)
            continue
        row = queue[i]
        i += 1
        sid = station_id(row)

        def _drop(reason: str) -> None:
            dropped.append({**row, "station": sid, "reason": reason,
                            "distance_km": round(row["distance_km"], 1)})

        wf_cache = _cache_dir(event.public_id) / f"{sid}.{row['band']}.mseed"
        try:
            st = read(str(wf_cache)) if wf_cache.exists() else None
            if st is None or st[0].stats.starttime > origin - pre_s + 5:
                st = client.get_waveforms(
                    network=row["network"], station=row["station"],
                    location=row["location"], channel=f"{row['band']}?",
                    starttime=origin - pre_s,
                    endtime=origin + config.TIME_AFTER_S + config.TIME_BEFORE_S,
                    attach_response=False)
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
            st.remove_response(inventory=inv, pre_filt=config.RESPONSE_PRE_FILT,
                               output="DISP", zero_mean=True)
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
        st.filter("bandpass", freqmin=fmin, freqmax=fmax,
                  corners=config.FILTER_CORNERS, zerophase=True)
        st.taper(max_percentage=0.05)

        # rule 3: per-component SNR on the filtered, un-decimated traces,
        # over the window that will be inverted
        kin_end = (row["distance_km"] / P.station.groupVelKms
                   + config.window_tail_s(event.prelim_mag))
        try:
            snr = component_snr(st, origin, row["distance_km"], kin_end)
        except AssertionError as e:
            _drop(f"no data: {e}")
            continue
        row["snr"] = snr
        # a station is usable when its BEST component carries signal: a
        # station near a nodal plane of one component still constrains
        # the mechanism through the others (2026p669681: WLRZ, the
        # best-fitting station at own VR 72, scored Z 1.9 R 1.5 T 3.5 -
        # the median rule threw it away). The fit decides the rest.
        row["snr_med"] = round(float(np.median(list(snr.values()))), 2)
        row["snr_best"] = round(float(max(snr.values())), 2)

        for tr in st:
            factor = int(round(tr.stats.sampling_rate * config.DT))
            assert factor >= 1, f"{sid}: sampling rate {tr.stats.sampling_rate}"
            tr.decimate(factor=factor, strict_length=False, no_filter=True)
            tr.resample(1.0 / config.DT, strict_length=False, no_filter=True)
            tr.trim(origin - config.TIME_BEFORE_S, origin + config.TIME_AFTER_S,
                    nearest_sample=True)
            tr.data = 100.0 * tr.data  # m -> cm (TDMT convention)
        expected = config.TIME_BEFORE_S + config.TIME_AFTER_S + 1
        if {tr.stats.npts for tr in st} != {expected}:
            _drop(f"no data: trim gave npts "
                  f"{ {tr.stats.npts for tr in st} }, want {expected}")
            continue

        row["window_end_s"] = signal_window_end_s(
            st, origin, row["distance_km"], event.prelim_mag)
        row["filter_hz"] = [fmin, fmax]
        row["peak_x_dist"] = float(
            abs(st.select(component="Z")[0].data).max() * row["distance_km"])
        ok = row["snr_best"] >= P.station.snrMin
        prefix = "" if ok else "rejected_"
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
            tr.write(str(workdir / f"{prefix}{sid}.{comp}.dat"), format="SAC")
        if not ok:
            _drop(f"low SNR: best component SNR {row['snr_best']:g} < "
                  f"{P.station.snrMin:g} (Z {snr['Z']:g} R {snr['R']:g} "
                  f"T {snr['T']:g})")
            continue
        if stages is not None:
            final_copy = st.copy()
            for tr in final_copy:
                tr.stats.distance = row["distance_km"] * 1000.0
            stages["final"] += final_copy
        usable.append(row)

    def _reject_files(r: dict) -> None:
        for comp in "ZRT":
            p = workdir / f"{station_id(r)}.{comp}.dat"
            if p.exists():
                p.rename(workdir / f"rejected_{p.name}")

    # rule 4: amplitude outliers (high side only)
    if len(usable) >= 3:
        med = float(np.median([r["peak_x_dist"] for r in usable]))
        for r in usable:
            r["amp_ratio"] = round(r["peak_x_dist"] / med, 3)
        flagged = [r for r in usable
                   if r["amp_ratio"] > P.station.amplitudeMaxRatio]
        for r in flagged:
            _reject_files(r)
            dropped.append({**r, "station": station_id(r),
                            "reason": f"amplitude outlier: {r['amp_ratio']:.1f}x "
                                      f"the network median",
                            "distance_km": round(r["distance_km"], 1)})
            print(f"  amplitude outlier: {station_id(r)} {r['amp_ratio']:.0f}x")
        usable = [r for r in usable if r not in flagged]

    # rule 5: azimuth-balanced selection
    pool, left = select_by_sector(usable)
    for r in left:
        _reject_files(r)
        dropped.append({**r, "station": station_id(r),
                        "reason": f"not selected: {P.station.maxStations} "
                                  f"seats filled by sector round-robin "
                                  f"(best SNR {r['snr_best']:g})",
                        "distance_km": round(r["distance_km"], 1)})
    pool.sort(key=lambda r: r["distance_km"])
    return pool, dropped


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--event", required=True, help="GeoNet publicID")
    ap.add_argument("--band", default=None,
                    help="passband in Hz, e.g. 0.02-0.10 (default: first "
                         "band of the magnitude menu)")
    ap.add_argument("--out", default=None,
                    help="working directory (default: <events>/<id>/task2)")
    args = ap.parse_args()
    from geonet import get_event
    ev = get_event(args.event)
    band = (tuple(float(x) for x in args.band.split("-")) if args.band
            else config.band_candidates(ev.prelim_mag)[0])
    out = Path(args.out) if args.out else config.EVENTS_DIR / ev.public_id / "task2"
    print(f"{ev.public_id} M{ev.prelim_mag:.1f} {ev.locality}: band "
          f"{config.band_tag(band)} -> {out}")
    pool, dropped = fetch_and_process(ev, out, band)
    print(f"\n{'station':14s} {'dist':>5s} {'az':>4s} sec {'SNR Z/R/T':>15s} "
          f"{'med':>5s} {'win':>4s} {'amp':>5s}")
    for r in pool:
        print(f"{station_id(r):14s} {r['distance_km']:5.0f} {r['azimuth']:4.0f} "
              f"{r['sector']:3d} {r['snr']['Z']:5.1f}/{r['snr']['R']:5.1f}/"
              f"{r['snr']['T']:5.1f} {r['snr_med']:5.1f} {r['window_end_s']:4d} "
              f"{r.get('amp_ratio', float('nan')):5.2f}")
    print(f"\n{len(pool)} selected, {len(dropped)} not:")
    for d in dropped:
        print(f"  {d['station']:14s} {d['reason']}")
