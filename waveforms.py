"""SAC preparation — attribution and provenance.

This module was compiled with Claude (Anthropic) assistance. The processing
chain follows the mttime example notebooks by Andrea Chiang (LLNL),
specifically 01_Data_Processing and 02_Prepare_Data_and_Synthetics_For_
Inversion: https://github.com/LLNL/mttime/tree/master/examples/notebooks
(mttime: https://github.com/LLNL/mttime, LLNL-CODE-814839), implemented
with ObsPy (https://github.com/obspy/obspy).

Deviations from the original notebooks: GeoNet NRT/archive FDSN sources
with retry logic; automated station selection (distance window, HH?>BH?
priority, 3x-depth rule); per-station SNR gate (signal/pre-event noise in
the inversion passband, threshold in config.py); magnitude-dependent
filter-band menu applied per event rather than fixed corners; fail-loud
drop accounting (every rejected station recorded with a reason).

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

from obspy import Stream, UTCDateTime
from obspy.core.util.attribdict import AttribDict
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees

import config
from geonet import Event, fdsn_client


def select_stations(client, event: Event, origin: UTCDateTime):
    """Broadband NZ stations within the working distance range.

    Returns (inventory, station_rows) where station_rows is a list of dicts
    with one preferred (channel-band, location) per station, nearest first.
    """
    inv = client.get_stations(
        network=config.NETWORK,
        channel=",".join(config.CHANNEL_PRIORITY),
        latitude=event.latitude,
        longitude=event.longitude,
        maxradius=kilometers2degrees(config.MAX_STATION_DIST_KM),
        level="response",
        starttime=origin,
        endtime=origin + config.TIME_AFTER_S,
    )
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
            if not (config.MIN_STATION_DIST_KM <= dist_km <= config.MAX_STATION_DIST_KM):
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


def fetch_and_process(
    event: Event, workdir: Path, band_hz: tuple[float, float],
    stages: dict | None = None,
) -> tuple[list[dict], list[dict]]:
    """Download and pre-process waveforms for one event.

    Writes SAC files ``<workdir>/NET.STA.LOC.{Z,R,T}.dat`` and returns
    (used_rows, dropped) where dropped is a list of {station, reason}.

    If ``stages`` (a dict) is passed, copies of every used station's stream
    are stored under keys "raw", "displacement", "final" for QC plotting.
    """
    origin = UTCDateTime(event.origin_time)
    client = fdsn_client(origin)
    inv, rows = select_stations(client, event, origin)

    workdir.mkdir(parents=True, exist_ok=True)
    fmin, fmax = band_hz
    if stages is not None:
        stages.update({"raw": Stream(), "displacement": Stream(), "final": Stream()})

    used, dropped = [], []
    for row in rows:
        if len(used) >= config.MAX_STATIONS:
            break
        sid = f"{row['network']}.{row['station']}.{row['location']}"

        def _drop(reason):
            dropped.append({
                "station": sid, "reason": reason,
                "latitude": row["latitude"], "longitude": row["longitude"],
                "distance_km": round(row["distance_km"], 1),
            })
        # far-field / point-source guard: distance > 3x depth (skip for
        # placeholder depths, which are meaningless)
        if (
            event.depth_km not in config.PLACEHOLDER_DEPTHS_KM
            and row["distance_km"] < config.MIN_DIST_DEPTH_RATIO * event.depth_km
        ):
            _drop("distance < 3x source depth")
            continue
        try:
            st = client.get_waveforms(
                network=row["network"],
                station=row["station"],
                location=row["location"],
                channel=f"{row['band']}?",
                starttime=origin - 5 * config.TIME_BEFORE_S,
                endtime=origin + config.TIME_AFTER_S + config.TIME_BEFORE_S,
                attach_response=False,
            )
        except Exception as e:  # noqa: BLE001 - record and move on, loudly
            _drop(f"download failed: {e}")
            continue

        st.merge(method=0)
        if any(hasattr(tr.data, "mask") for tr in st) or len(st) < 3:
            _drop("gaps or <3 components")
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
            _drop(f"response/rotation failed: {e}")
            continue

        if len(st.select(component="Z")) != 1 or len(st) != 3:
            _drop("not exactly 3 ZNE components")
            continue

        # SNR gate in the inversion passband: RMS(signal)/RMS(pre-event noise)
        snr_st = st.copy()
        snr_st.filter(
            "bandpass", freqmin=fmin, freqmax=fmax,
            corners=config.FILTER_CORNERS, zerophase=True,
        )
        snr = _snr(snr_st, origin)
        if snr < config.MIN_SNR:
            _drop(f"SNR {snr:.1f} < {config.MIN_SNR}")
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
            _drop(f"trim gave npts {npts}, want {expected}")
            continue

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
            tr.write(str(workdir / f"{sid}.{comp}.dat"), format="SAC")

        if stages is not None:
            final_copy = st.copy()
            for tr in final_copy:
                tr.stats.distance = row["distance_km"] * 1000.0
            stages["final"] += final_copy

        row["snr"] = round(snr, 2)
        row["filter_hz"] = [fmin, fmax]
        used.append(row)

    assert used, f"all {len(rows)} candidate stations dropped: {dropped}"
    return used, dropped


def _snr(st: Stream, origin: UTCDateTime) -> float:
    """Min over components of RMS(signal window)/RMS(noise window)."""
    import numpy as np

    snrs = []
    for tr in st:
        noise = tr.slice(origin - 4 * config.TIME_BEFORE_S, origin - 10).data
        signal = tr.slice(origin, origin + config.TIME_AFTER_S).data
        assert len(noise) > 10 and len(signal) > 10, "SNR windows too short"
        rms_n = float(np.sqrt(np.mean(noise**2)))
        rms_s = float(np.sqrt(np.mean(signal**2)))
        assert rms_n > 0, "zero noise RMS (dead channel?)"
        snrs.append(rms_s / rms_n)
    return min(snrs)
