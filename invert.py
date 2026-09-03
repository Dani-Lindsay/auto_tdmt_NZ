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
    # INDEPENDENT BY CHOICE (2026-09-03): the full grid is always searched.
    # Bounding the search around GeoNet's depth forces agreement and
    # destroys the catalogue's value as an independent check; GeoNet's
    # depth and uncertainty are recorded for comparison only.
    return lib


def write_mtinv(
    event: Event, stations: list[dict], depths: list[float],
    event_dir: Path, green_dir: Path,
) -> Path:
    # adaptive record length: each station's window ends shortly after
    # its surface-wave train. waveforms.py computes the signal-aware end
    # (kinematic minimum, envelope-extended for slow paths) and stores it
    # as window_end_s; the kinematic formula is the fallback for rows
    # that lack it.
    tail = config.window_tail_s(event.prelim_mag)
    npts_col = [
        int(min(config.INV_NPTS, max(
            config.WINDOW_MIN_S,
            config.TIME_BEFORE_S + r["window_end_s"]
            if "window_end_s" in r else
            config.TIME_BEFORE_S
            + r["distance_km"] / config.WINDOW_GROUP_VEL_KMS + tail)))
        for r in stations
    ]
    frame = {
        "station": [
            f"{r['network']}.{r['station']}.{r['location']}" for r in stations
        ],
        "distance": [r["gf_distance_km"] for r in stations],
        "azimuth": [round(r["azimuth"], 2) for r in stations],
        "ts": config.TIME_BEFORE_S,
        "npts": npts_col,
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
    if vrs[i_max] < config.DC_TIEBREAK_MIN_VR:
        return i_max  # junk-grade fit: DC differences are noise
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
    """Core-plus-admission station selection (calibrated against manual
    keep/toss labelling, 2026-09-03):

    1. the CORE — stations whose peak/noise clears PEAK_NOISE_CORE — gets
       a full depth search on its own, so the reference solution is never
       polluted by marginal data (the VR-only experiment showed noise
       traces earn chance VR against a solution their own noise corrupted);
    2. each remaining candidate ("yellow") is added ALONE at the core's
       preferred depth and admitted only if the core-dominated solution
       predicts its waveform: own station VR >= CANDIDATE_STATION_VR_MIN.
       %DC never enters the decision (a noise station can inflate DC);
    3. full depth search with core + admitted;
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

    # 1: core-only full depth search (top-up by peak/noise if the core is
    # thinner than the minimum station count)
    ranked = sorted(stations, key=lambda r: -r.get("pk_n", 0.0))
    core = [r for r in ranked
            if r.get("tier", "core") == "core"]
    if len(core) < config.MIN_STATIONS_USED:
        for r in ranked:
            if r not in core:
                core.append(r)
            if len(core) >= config.MIN_STATIONS_USED:
                break
    cands = [r for r in ranked if r not in core]
    core.sort(key=lambda r: r["distance_km"])
    inv1 = _solve(core, depths)
    mt0 = inv1.moment_tensors[inv1.preferred_tensor_id]
    depth0 = float(mt0.depth)
    print(f"core: {len(core)} stations, depth {depth0:g} km, "
          f"VR {float(mt0.total_VR):.1f}")

    # 2: one-at-a-time candidate admission at the core's preferred depth.
    # Sparse cores relax the floor: 3-4 stations constrain a mechanism
    # poorly, so extra azimuth coverage is worth a weaker individual fit.
    vr_floor = (config.CANDIDATE_STATION_VR_MIN
                if len(core) >= config.SPARSE_CORE_COUNT
                else config.CANDIDATE_STATION_VR_MIN_SPARSE)
    admitted = []
    for c in cands:
        inv_c = _solve(core + [c], [depth0])
        mt_c = inv_c.moment_tensors[inv_c.preferred_tensor_id]
        own_vr = float(mt_c.station_table.iloc[-1].VR)
        if own_vr >= vr_floor:
            c["admission_vr"] = round(own_vr, 1)
            admitted.append(c)
        else:
            rejected.append({
                "station": _sid(c),
                "reason": f"not predicted by core solution: station VR "
                          f"{own_vr:.0f} at {depth0:g} km < {vr_floor:g} "
                          f"(peak/noise {c.get('pk_n', 0.0):.1f})",
            })
    if cands:
        print(f"candidates: {len(admitted)}/{len(cands)} admitted "
              f"({[r['station'] for r in admitted]})")

    # 3: core + admitted, full depth search
    current = sorted(core + admitted, key=lambda r: r["distance_km"])
    inv2 = _solve(current, depths) if admitted else inv1

    # 3b: azimuth-coverage cap (2026p033598 review: an over-full network
    # dilutes %DC). Keep the best-VR station in each 45 deg sector plus the
    # top-4 VR overall — coverage first, then fit — and re-search.
    if len(current) > config.MAX_USED_STATIONS:
        mt2 = inv2.moment_tensors[inv2.preferred_tensor_id]
        vrs = {r.station: float(r.VR)
               for r in mt2.station_table.itertuples()}
        by_sector: dict[int, list] = {}
        for r in current:
            by_sector.setdefault(_sector(r), []).append(r)
        keep_ids = {_sid(max(rows_, key=lambda r: vrs.get(_sid(r), -1e9)))
                    for rows_ in by_sector.values()}
        for r in sorted(current,
                        key=lambda r: -vrs.get(_sid(r), -1e9))[:4]:
            keep_ids.add(_sid(r))
        for r in current:
            if _sid(r) not in keep_ids:
                rejected.append({
                    "station": _sid(r),
                    "reason": "azimuth-coverage cap: sector-best + top-4 "
                              f"VR kept {len(keep_ids)} of {len(current)}; "
                              f"station VR {vrs.get(_sid(r), 0.0):.0f}",
                })
        if len(keep_ids) < len(current):
            print(f"coverage cap: {len(current)} -> {len(keep_ids)} "
                  "stations")
            current = sorted(
                (r for r in current if _sid(r) in keep_ids),
                key=lambda r: r["distance_km"])
            inv2 = _solve(current, depths)
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
        # anti-fitting cull: NEGATIVE own VR means the station fits worse
        # than silence — it only steers the tensor. Dropped without the
        # joint-gain test (distance weighting can hide the joint cost:
        # 2026p238013 NNZ sat at station VR -41 behind a passing joint VR).
        if vrs.get(_sid(worst), 0.0) < 0.0:
            rejected.append({
                "station": _sid(worst),
                "reason": f"anti-fitting: station VR "
                          f"{vrs.get(_sid(worst), 0.0):.0f} at "
                          f"{depth_pref:g} km",
            })
            current = [r for r in current if r is not worst]
            greedy += 1
            continue
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
    """Angle between two DC tensors via the normalised tensor inner
    product (nodal-plane-choice independent). Retained for continuity;
    the preferred mechanism-comparison metric is
    ``min_rotation_angle_deg`` below (J. Townend's recommendation)."""
    m1, m2 = dc_tensor(*a), dc_tensor(*b)
    cos = np.sum(m1 * m2) / (np.linalg.norm(m1) * np.linalg.norm(m2))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def _mech_rotation_matrix(strike: float, dip: float,
                          rake: float) -> np.ndarray:
    """Focal mechanism as a rotation matrix with respect to geographic
    (NED) coordinates (Walsh et al. 2009, GJI 176, eqs 1-3 sense).
    Columns are the principal axes T, P and null B (right-handed,
    T x P = B), built from the Aki & Richards (1980) slip vector u and
    fault normal n. In THIS frame the double-couple symmetry group is
    exactly the diagonal sign flips, and the conjugate-plane
    parameterisation (u and n swapped) maps to one of them — which is
    why the minimum below is independent of the nodal-plane choice."""
    phi, delta, lam = np.radians([strike, dip, rake])
    n = np.array([-np.sin(delta) * np.sin(phi),
                  np.sin(delta) * np.cos(phi),
                  -np.cos(delta)])
    u = np.array([
        np.cos(lam) * np.cos(phi) + np.cos(delta) * np.sin(lam) * np.sin(phi),
        np.cos(lam) * np.sin(phi) - np.cos(delta) * np.sin(lam) * np.cos(phi),
        -np.sin(lam) * np.sin(delta)])
    t = (u + n) / np.sqrt(2.0)          # tension axis
    p = (u - n) / np.sqrt(2.0)          # pressure axis
    b = np.cross(t, p)                  # null axis; det(+1) by order
    return np.column_stack([t, p, b])


def min_rotation_angle_deg(a, b) -> float:
    """Minimum rotation aligning two double couples (J. Townend's
    recipe, 2026-08-20; Townend et al. 2012 supplement eq. 1; Walsh,
    Arnold & Townend 2009; cf. Kagan 1991):
    angle = arccos((tr(R1^T R2) - 1) / 2), minimised over the DC
    symmetry group (180-degree rotations about slip, null and normal),
    so the result is independent of which nodal plane parameterises
    either mechanism. ``a``/``b`` are (strike, dip, rake) tuples."""
    r1, r2 = _mech_rotation_matrix(*a), _mech_rotation_matrix(*b)
    best = 180.0
    for sym in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)):
        cos = (np.trace(r1.T @ (r2 * np.array(sym))) - 1.0) / 2.0
        best = min(best, float(np.degrees(
            np.arccos(np.clip(cos, -1.0, 1.0)))))
    return best


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
                    min_rotation_angle_deg(ref, tuple(fps[0])), 1),
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

    def _plateau_span(rr):
        srt = sorted(rr, key=lambda r: r["depth_km"])
        vmax = max(r["vr"] for r in srt)
        i_max = max(range(len(srt)), key=lambda i: srt[i]["vr"])
        lo = i_max
        while lo > 0 and srt[lo - 1]["vr"] >= vmax - config.PREFER_DC_VR_TOLERANCE:
            lo -= 1
        hi = i_max
        while hi < len(srt) - 1 and srt[hi + 1]["vr"] >= vmax - config.PREFER_DC_VR_TOLERANCE:
            hi += 1
        return srt[hi]["depth_km"] - srt[lo]["depth_km"]

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
            "plateau_km": round(_plateau_span(rows), 1),
            "depth_unconstrained": _plateau_span(rows) > 10.0,
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
