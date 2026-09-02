"""Inversion driver — attribution and provenance.

This module was compiled with Claude (Anthropic) assistance. The mtinv.in
format and inversion invocation follow the mttime example notebooks by
Andrea Chiang (LLNL): https://github.com/LLNL/mttime (LLNL-CODE-814839).

Deviations: automated depth search over the GF library grid; low-VR
station rejection pass; preferred-solution rule (contiguous VR plateau,
then max %DC) replacing manual inspection; machine-readable solution.json
with full provenance.

Moment tensor inversion via mttime: mtinv.in generation, depth search,
quality gates, and solution serialization.

The mtinv.in file (exactly the format of the mttime example notebooks) is
written into the event directory both to drive the inversion and as
provenance — anyone can rerun the event from its archived directory.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import config
import greens
from geonet import Event


def search_depths(event: Event, model: str) -> list[float]:
    """Depth grid to invert. Full library range when GeoNet's depth is a
    fixed placeholder (their initial depths are unreliable); otherwise the
    library depths within DEPTH_SEARCH_MARGIN_KM of the hypocentre."""
    lib = greens.available_depths(model)
    assert lib, "GF library is empty — build it first (greens.py --build)"
    if event.depth_km in config.PLACEHOLDER_DEPTHS_KM:
        return lib
    lo = event.depth_km - config.DEPTH_SEARCH_MARGIN_KM
    hi = event.depth_km + config.DEPTH_SEARCH_MARGIN_KM
    picked = [d for d in lib if lo <= d <= hi]
    return picked if picked else lib


def write_mtinv(
    event: Event, stations: list[dict], depths: list[float],
    event_dir: Path, green_dir: Path,
) -> Path:
    frame = {
        "station": [
            f"{r['network']}.{r['station']}.{r['location']}" for r in stations
        ],
        "distance": [r["gf_distance_km"] for r in stations],
        "azimuth": [round(r["azimuth"], 2) for r in stations],
        "ts": config.TIME_BEFORE_S,
        "npts": config.INV_NPTS,
        "dt": config.DT,
        "used": 1,
        "longitude": [r["longitude"] for r in stations],
        "latitude": [r["latitude"] for r in stations],
    }
    df = pd.DataFrame(frame)
    headers = dict(
        datetime=event.origin_time,
        longitude=event.longitude,
        latitude=event.latitude,
        depth=",".join(f"{d:.4f}" for d in depths),
        # mttime prepends the mtinv.in directory -> paths must be relative
        path_to_data=".",
        path_to_green=green_dir.name,
        green="herrmann",
        components="ZRT",
        degree=config.INVERSION_DEGREE,
        weight="distance",
        plot=0,
        correlate=1,  # solve for per-station time shift (zcor)
    )
    path = event_dir / "mtinv.in"
    with open(path, "w") as f:
        for key, value in headers.items():
            f.write(f"{key:<15}{value}\n")
        f.write(df.to_string(index=False))
    return path


def pick_preferred(vr_pdc: list[tuple[float, float]],
                   contiguous: bool = True) -> int:
    """Index of the preferred solution.

    Hierarchy: VR first — candidates are the depths within
    PREFER_DC_VR_TOLERANCE of the VR maximum; among those, take the highest
    %DC. With ``contiguous=True`` (the depth search, where entries are
    ordered by depth) candidates are restricted to the CONTIGUOUS plateau
    containing the VR maximum: a bimodal VR curve must not let a
    disconnected deep lobe that grazes the tolerance steal the pick on DC.
    Band selection passes ``contiguous=False`` (few, unordered candidates).
    """
    vrs = [vr for vr, _ in vr_pdc]
    i_max = max(range(len(vrs)), key=lambda i: vrs[i])
    floor = vrs[i_max] - config.PREFER_DC_VR_TOLERANCE
    if contiguous:
        lo = i_max
        while lo > 0 and vrs[lo - 1] >= floor:
            lo -= 1
        hi = i_max
        while hi < len(vrs) - 1 and vrs[hi + 1] >= floor:
            hi += 1
        cands = list(range(lo, hi + 1))
    else:
        cands = [i for i, vr in enumerate(vrs) if vr >= floor]
    return max(cands, key=lambda i: vr_pdc[i][1])


def run_inversion(mtinv_path: Path):
    """Run mttime over the configured depths; returns the Inversion object
    with preferred_tensor_id set by the VR+%DC rule (not mttime's VR-max)."""
    import mttime

    cfg = mttime.Configure(path_to_file=str(mtinv_path))
    inv = mttime.Inversion(config=cfg)
    inv.invert()
    assert inv.moment_tensors, "mttime produced no solutions"
    inv.preferred_tensor_id = pick_preferred(
        [(float(mt.total_VR), float(mt.pdc)) for mt in inv.moment_tensors]
    )
    return inv


def invert_with_rejection(
    event: Event, stations: list[dict], depths: list[float],
    event_dir: Path, green_dir: Path,
):
    """Invert, drop stations whose individual VR at the preferred depth is
    below STATION_VR_FLOOR, and invert once more without them (EPS207 §3.1).

    Returns (inversion, kept_stations, rejected) where rejected is a list of
    {station, reason} entries.
    """
    mtinv_path = write_mtinv(event, stations, depths, event_dir, green_dir)
    inv = run_inversion(mtinv_path)

    pref = inv.moment_tensors[inv.preferred_tensor_id]
    table = pref.station_table
    bad = table[table.VR < config.STATION_VR_FLOOR]
    if len(bad) == 0:
        return inv, stations, []

    bad_ids = list(bad.station)
    keep = [
        r for r in stations
        if f"{r['network']}.{r['station']}.{r['location']}" not in bad_ids
    ]
    assert len(keep) >= config.MIN_STATIONS_USED, (
        f"station rejection would leave {len(keep)} stations "
        f"(< {config.MIN_STATIONS_USED}); rejected: {bad_ids}"
    )
    rejected = [
        {
            "station": row.station,
            "reason": f"station VR {row.VR:.0f} < {config.STATION_VR_FLOOR:g} "
                      f"at first-pass depth {pref.depth:g} km",
        }
        for row in bad.itertuples()
    ]
    print(f"rejecting {len(rejected)} low-VR stations: {bad_ids}; re-inverting")
    mtinv_path = write_mtinv(event, keep, depths, event_dir, green_dir)
    return run_inversion(mtinv_path), keep, rejected


def summarize(inv, event: Event, stations: list[dict], dropped: list[dict],
              model: str) -> dict:
    """Serialize the depth search + preferred solution with provenance."""
    import mttime
    import obspy

    def _mt_row(mt) -> dict:
        # fps: two (strike, dip, rake) nodal planes from the DC part
        fps = np.asarray(mt.fps, dtype=float)
        assert fps.shape == (2, 3), f"unexpected fps shape {fps.shape}"
        return dict(
            depth_km=float(mt.depth),
            mw=float(mt.mw),
            m0_dyne_cm=float(mt.mo),
            vr=float(mt.total_VR),
            pdc=float(mt.pdc),  # already percent
            pclvd=float(mt.pclvd),
            piso=float(mt.piso),
            plane1=dict(zip(("strike", "dip", "rake"), fps[0])),
            plane2=dict(zip(("strike", "dip", "rake"), fps[1])),
            # full tensor per depth so the depth-sensitivity figure can
            # draw the true deviatoric beachball at each trial depth
            tensor_rtp_dyne_cm={
                k: float(v)
                for k, v in mt.get_tensor_elements(basis="RTP").items()
            },
        )

    # per-station solved time shifts (zcor) and individual VR at the
    # preferred depth — the audit trail for velocity-model error absorption
    pref_mt = inv.moment_tensors[inv.preferred_tensor_id]
    table = {r.station: r for r in pref_mt.station_table.itertuples()}
    for st in stations:
        sid = f"{st['network']}.{st['station']}.{st['location']}"
        row = table.get(sid)
        if row is not None:
            st["zcor_s"] = float((row.ts - config.TIME_BEFORE_S) * config.DT)
            st["station_vr"] = float(row.VR)

    rows = [_mt_row(mt) for mt in inv.moment_tensors]
    best_vr = max(rows, key=lambda r: r["vr"])
    best_dc = max(rows, key=lambda r: r["pdc"])
    pref = inv.moment_tensors[inv.preferred_tensor_id]

    solution = {
        "event": event.to_dict(),
        "preferred": {
            **_mt_row(pref),
            "tensor_dyne_cm": {
                k: float(v) for k, v in pref.get_tensor_elements().items()
            },
            # Harvard/GCMT r,theta,phi convention — feeds psmeca -Sz so the
            # plotted beachball is the TRUE deviatoric mechanism (CLVD
            # included), not just the closest double couple
            "tensor_rtp_dyne_cm": {
                k: float(v)
                for k, v in pref.get_tensor_elements(basis="RTP").items()
            },
        },
        "depth_search": rows,
        # EPS207 §3.3: %DC max can be more diagnostic than VR max; flag when
        # the two picks disagree so a human looks before trusting the depth.
        "depth_pick_flags": {
            "vr_max_depth_km": best_vr["depth_km"],
            "dc_max_depth_km": best_dc["depth_km"],
            "vr_dc_agree": best_vr["depth_km"] == best_dc["depth_km"],
        },
        "stations_used": stations,
        "stations_dropped": dropped,
        "provenance": {
            "velocity_model": model,
            "gf_version": config.GF_VERSION,
            "mttime_version": mttime.__version__,
            "obspy_version": obspy.__version__,
            "degree": config.INVERSION_DEGREE,
            "weight": "distance",
            "correlate": True,
            "npts": config.INV_NPTS,
            "dt_s": config.DT,
            "preferred_rule": (
                f"max %DC on the contiguous depth plateau within "
                f"{config.PREFER_DC_VR_TOLERANCE:g} VR points of the maximum"
            ),
        },
    }
    solution["quality"] = quality_gates(solution)
    return solution


def quality_gates(solution: dict) -> dict:
    """Publication gates — fail loud, publish nothing below the bar."""
    pref = solution["preferred"]
    n_used = len(solution["stations_used"])
    azimuths = sorted(s["azimuth"] for s in solution["stations_used"])
    gaps = [
        (azimuths[(i + 1) % len(azimuths)] - a) % 360.0
        for i, a in enumerate(azimuths)
    ]
    az_gap = max(gaps) if gaps else 360.0

    checks = {
        "min_stations": n_used >= config.MIN_STATIONS_USED,
        "vr_floor": pref["vr"] >= config.MIN_VR_PUBLISH,
        "az_gap_ok": az_gap <= config.MAX_AZ_GAP_DEG,
    }
    # warnings are recorded but do not block: a shallow crustal event
    # legitimately prefers the shallowest grid depth
    warnings = {
        "depth_at_grid_edge": pref["depth_km"] in (
            min(config.GF_DEPTHS_KM), max(config.GF_DEPTHS_KM),
        ),
    }
    return {
        "n_stations_used": n_used,
        "azimuthal_gap_deg": round(az_gap, 1),
        "checks": checks,
        "warnings": warnings,
        "passed": all(checks.values()),
    }


def save_solution(solution: dict, event_dir: Path) -> Path:
    out = event_dir / "solution.json"
    with open(out, "w") as f:
        json.dump(solution, f, indent=2)
    return out
