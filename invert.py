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
    """Depth-profile station selection (clean hierarchy):

    1. FULL pool, FULL depth search — one inversion pass yields every
       station's VR at every trial depth (its depth profile);
    2. a station whose BEST VR across all depths never clears the
       abundance-conditional floor (30% while >8 stations, 20% >5, 10%
       below) is consistently bad and dropped — one misaligned depth
       cannot condemn a station, and no provisional depth is ever assumed;
       a sector's sole representative survives unless even its best VR is
       negative;
    3. full depth search with the survivors -> preferred solution;
    4. greedy earn-your-seat pass at the preferred depth: the worst
       station is test-dropped while the joint VR improves by >= 2;
    5. if the greedy pass removed anyone, one final full depth search.

    Returns (inversion, kept_stations, rejected).
    """
    rejected: list[dict] = []

    def _sid(r):
        return f"{r['network']}.{r['station']}.{r['location']}"

    def _sector(r):
        return int(r["azimuth"] // 45) % 8

    def _solve(rows, dd):
        return run_inversion(
            write_mtinv(event, rows, dd, event_dir, green_dir))

    def _floor(n):
        if n > config.RICH_STATION_COUNT:
            return config.STATION_VR_FLOOR_RICH
        if n > config.MID_STATION_COUNT:
            return config.STATION_VR_FLOOR_MID
        return config.STATION_VR_FLOOR

    # 1: full pool, full depth search -> per-station VR-vs-depth profiles
    inv1 = _solve(stations, depths)
    best_vr: dict[str, float] = {}
    for mt in inv1.moment_tensors:
        for r in mt.station_table.itertuples():
            v = float(r.VR)
            if v > best_vr.get(r.station, -1e9):
                best_vr[r.station] = v

    # 2: drop the consistently bad (best-over-depth VR below the floor)
    floor0 = _floor(len(stations))
    sectors = {}
    for r in stations:
        sectors.setdefault(_sector(r), []).append(r)
    candidates = sorted(
        (r for r in stations if best_vr.get(_sid(r), 0.0) < floor0),
        key=lambda r: best_vr.get(_sid(r), 0.0))
    current = list(stations)
    for r in candidates:
        if len(current) <= config.MIN_STATIONS_USED:
            break
        sole = len([x for x in sectors[_sector(r)] if x in current]) == 1
        if sole and best_vr.get(_sid(r), 0.0) >= 0.0:
            continue
        current.remove(r)
        rejected.append({
            "station": _sid(r),
            "reason": f"consistently bad: best VR "
                      f"{best_vr.get(_sid(r), 0.0):.0f} over "
                      f"{len(depths)} depths < floor {floor0:g}",
        })
    if rejected:
        print(f"dropped {len(rejected)} consistently-bad stations: "
              f"{[d['station'] for d in rejected]}")

    # 3: survivors, full depth search
    inv2 = _solve(current, depths) if rejected else inv1
    depth_pref = float(
        inv2.moment_tensors[inv2.preferred_tensor_id].depth)

    # 4: greedy earn-your-seat at the preferred depth
    greedy = 0
    while len(current) > config.MIN_STATIONS_USED:
        inv_d = _solve(current, [depth_pref])
        base = float(
            inv_d.moment_tensors[inv_d.preferred_tensor_id].total_VR)
        vrs = {r.station: float(r.VR)
               for r in inv_d.moment_tensors[0].station_table.itertuples()}
        sectors_left = {}
        for r in current:
            sectors_left.setdefault(_sector(r), []).append(r)
        droppable = [r for r in current
                     if not (len(sectors_left[_sector(r)]) == 1
                             and vrs.get(_sid(r), 0) >= 0.0)]
        if not droppable:
            break
        worst = min(droppable, key=lambda r: vrs.get(_sid(r), 0))
        trial = [r for r in current if r is not worst]
        inv_t = _solve(trial, [depth_pref])
        tot_t = float(
            inv_t.moment_tensors[inv_t.preferred_tensor_id].total_VR)
        if tot_t >= base + config.ELIMINATION_VR_GAIN:
            rejected.append({
                "station": _sid(worst),
                "reason": f"test-drop improved joint fit: VR {base:.0f} "
                          f"-> {tot_t:.0f} at {depth_pref:g} km",
            })
            current = trial
            greedy += 1
        else:
            break
    if greedy:
        print(f"greedy pass removed {greedy} more stations")
        current.sort(key=lambda r: r["distance_km"])
        inv2 = _solve(current, depths)

    return inv2, current, rejected


def dc_tensor(strike: float, dip: float, rake: float) -> np.ndarray:
    """Unit-moment DC tensor in NED (Aki & Richards 1980 eqs 4.84-4.89)."""
    p, d, r = np.radians([strike, dip, rake])
    sp, cp = np.sin(p), np.cos(p)
    sd, cd = np.sin(d), np.cos(d)
    sr, cr = np.sin(r), np.cos(r)
    s2p, c2p = np.sin(2 * p), np.cos(2 * p)
    s2d, c2d = np.sin(2 * d), np.cos(2 * d)
    mnn = -(sd * cr * s2p + s2d * sr * sp**2)
    mee = sd * cr * s2p - s2d * sr * cp**2
    mdd = -(mnn + mee)
    mne = sd * cr * c2p + 0.5 * s2d * sr * s2p
    mnd = -(cd * cr * cp + c2d * sr * sp)
    med = -(cd * cr * sp - c2d * sr * cp)
    return np.array([[mnn, mne, mnd], [mne, mee, med], [mnd, med, mdd]])


def tensor_angle_deg(a, b) -> float:
    """Angle between two DC tensors (nodal-plane-choice independent)."""
    m1, m2 = dc_tensor(*a), dc_tensor(*b)
    cos = np.sum(m1 * m2) / (np.linalg.norm(m1) * np.linalg.norm(m2))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def jackknife(event: Event, stations: list[dict], depth: float,
              band_dir: Path, green_dir: Path, full_plane1: dict) -> dict:
    """Leave-one-station-out stability test at the preferred depth (EPS207
    §3.2 in spirit; single-depth for speed). Mechanism stability is
    measured as the DC-tensor rotation of each subset solution relative to
    the full solution — immune to nodal-plane ordering flips between
    subsets. Rewrites mtinv.in per subset and restores it afterwards."""
    if len(stations) < 4:
        return {"n_subsets": 0,
                "note": "fewer than 4 stations; jackknife skipped"}
    ref = (full_plane1["strike"], full_plane1["dip"], full_plane1["rake"])
    subsets = []
    try:
        for i in range(len(stations)):
            subset = stations[:i] + stations[i + 1:]
            mtinv = write_mtinv(event, subset, [depth], band_dir, green_dir)
            mt = run_inversion(mtinv).moment_tensors[0]
            fps = np.asarray(mt.fps, dtype=float)
            subsets.append({
                "left_out": stations[i]["station"],
                "mw": round(float(mt.mw), 3),
                "pdc": round(float(mt.pdc), 1),
                "vr": round(float(mt.total_VR), 1),
                "tensor_rotation_deg": round(
                    tensor_angle_deg(ref, tuple(fps[0])), 1),
                "plane1": dict(zip(("strike", "dip", "rake"),
                                   [round(float(v), 1) for v in fps[0]])),
                "plane2": dict(zip(("strike", "dip", "rake"),
                                   [round(float(v), 1) for v in fps[1]])),
            })
    finally:
        # restore the full-station mtinv.in for reproducibility
        write_mtinv(event, stations, [depth], band_dir, green_dir)
    mws = [s["mw"] for s in subsets]
    dcs = [s["pdc"] for s in subsets]
    rots = [s["tensor_rotation_deg"] for s in subsets]
    return {
        "n_subsets": len(subsets),
        "mw_std": round(float(np.std(mws)), 3),
        "dc_std": round(float(np.std(dcs)), 1),
        "max_tensor_rotation_deg": round(max(rots), 1),
        "mean_tensor_rotation_deg": round(float(np.mean(rots)), 1),
        "subsets": subsets,
    }


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

    # BSL-style letter grade: one glance tells a scientist what the
    # solution is worth. A/B are emailed; C/D are archived for specialists.
    vr = pref["vr"]
    if vr >= 70 and n_used >= 5 and az_gap <= 180:
        grade = "A"
    elif vr >= 60 and n_used >= 3 and az_gap <= 270:
        grade = "B"
    elif vr >= 50 and n_used >= 2:
        grade = "C"
    else:
        grade = "D"

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
        "grade": grade,
        "checks": checks,
        "warnings": warnings,
        "passed": grade in ("A", "B"),
    }


def save_solution(solution: dict, event_dir: Path) -> Path:
    out = event_dir / "solution.json"
    with open(out, "w") as f:
        json.dump(solution, f, indent=2)
    return out
