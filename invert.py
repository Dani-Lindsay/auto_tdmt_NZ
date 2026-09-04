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
    event_dir: Path, green_dir: Path, correlate: bool = True,
) -> Path:
    """Write mttime's control file. ``correlate=False`` forbids per-station
    time shifts — used for the pass-1 survey, where an unbounded shift
    search would let noise traces slide into chance alignments and earn
    undeserved VR (mttime searches 0..(231-npts) samples, so a 60 s window
    can slide >100 s)."""
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
        correlate=1 if correlate else 0,  # per-station time shift (zcor)
    )
    path = event_dir / "mtinv.in"
    with open(path, "w") as f:
        for key, value in headers.items():
            f.write(f"{key:<15}{value}\n")
        f.write(df.to_string(index=False))
    return path


def plateau_indices(vrs: list[float], i_ref: int,
                    tol: float | None = None) -> list[int]:
    """Contiguous run of indices around ``i_ref`` whose VR stays within
    ``tol`` of vrs[i_ref] (the depth plateau)."""
    tol = config.PREFER_DC_VR_TOLERANCE if tol is None else tol
    floor = vrs[i_ref] - tol
    lo = i_ref
    while lo > 0 and vrs[lo - 1] >= floor:
        lo -= 1
    hi = i_ref
    while hi < len(vrs) - 1 and vrs[hi + 1] >= floor:
        hi += 1
    return list(range(lo, hi + 1))


def grid_edge_reference(vrs: list[float],
                        tol: float | None = None) -> tuple[int, bool]:
    """Reference index for the depth pick, with the grid-edge guard.

    A VR maximum on the first or last grid depth whose plateau is a single
    point is an artifact, not a depth: the smoothest Green's functions at
    the ends of the library absorb noise (2026p508890 rode a 58 km edge at
    VR 22.7 over the physical 8 km peak at VR 18.3 with DC 88-96; 520779,
    348732 and 300334 the same). When that happens, prefer the best
    INTERIOR local maximum within ``tol`` VR points; keep the edge only if
    nothing interior comes close. Returns (index, guard_applied)."""
    tol = config.GRID_EDGE_VR_TOLERANCE if tol is None else tol
    i_max = max(range(len(vrs)), key=lambda i: vrs[i])
    if len(vrs) < 3 or i_max not in (0, len(vrs) - 1):
        return i_max, False
    if len(plateau_indices(vrs, i_max)) > 1:
        return i_max, False  # a real plateau that happens to touch the end
    interior = [
        j for j in range(1, len(vrs) - 1)
        if vrs[j] >= vrs[j - 1] and vrs[j] >= vrs[j + 1]
        and vrs[j] >= vrs[i_max] - tol
    ]
    if not interior:
        return i_max, False
    return max(interior, key=lambda j: vrs[j]), True


def pick_preferred_detail(vr_pdc: list[tuple[float, float]],
                          contiguous: bool = True) -> tuple[int, bool]:
    """(index of the preferred solution, grid-edge guard applied?).

    Hierarchy: VR first — candidates are the depths within
    PREFER_DC_VR_TOLERANCE of the VR maximum; among those, take the highest
    %DC. With ``contiguous=True`` (the depth search, where entries are
    ordered by depth) candidates are restricted to the CONTIGUOUS plateau
    containing the reference: a bimodal VR curve must not let a
    disconnected deep lobe that grazes the tolerance steal the pick on DC,
    and the grid-edge guard rejects single-point maxima at the grid ends.
    Band selection passes ``contiguous=False`` (few, unordered candidates).
    """
    vrs = [vr for vr, _ in vr_pdc]
    if contiguous:
        i_ref, guarded = grid_edge_reference(vrs)
    else:
        i_ref, guarded = max(range(len(vrs)), key=lambda i: vrs[i]), False
    if vrs[i_ref] < config.DC_TIEBREAK_MIN_VR:
        return i_ref, guarded  # junk-grade fit: DC differences are noise
    if contiguous:
        # DC breaks a near-tie in VR only (DEPTH_DC_TOLERANCE), inside the
        # wider plateau that defines "the depth is not resolved"
        cands = plateau_indices(vrs, i_ref, config.DEPTH_DC_TOLERANCE)
    else:
        floor = vrs[i_ref] - config.PREFER_DC_VR_TOLERANCE
        cands = [i for i, vr in enumerate(vrs) if vr >= floor]
    return max(cands, key=lambda i: vr_pdc[i][1]), guarded


def pick_preferred(vr_pdc: list[tuple[float, float]],
                   contiguous: bool = True) -> int:
    """Index of the preferred solution (see pick_preferred_detail)."""
    return pick_preferred_detail(vr_pdc, contiguous)[0]


def azimuth_pair_ok(azimuths, min_sep: float | None = None) -> bool:
    """True when two stations are at least ``min_sep`` degrees apart in
    azimuth — the minimum geometry for resolving a mechanism at all."""
    min_sep = config.AZ_PAIR_MIN_DEG if min_sep is None else min_sep
    az = list(azimuths)
    for i, a in enumerate(az):
        for b in az[i + 1:]:
            sep = abs(a - b) % 360.0
            if min(sep, 360.0 - sep) >= min_sep:
                return True
    return False


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


class NoCoherentSolution(Exception):
    """Raised when no station set produces a coherent solution. Carries the
    evidence so the event can be archived honestly instead of publishing a
    junk mechanism."""

    def __init__(self, stage: str, reason: str, selection: dict,
                 pool: list[dict], rejected: list[dict]):
        super().__init__(f"{stage}: {reason}")
        self.stage = stage
        self.reason = reason
        self.selection = selection
        self.pool = pool
        self.rejected = rejected


def station_id(row: dict) -> str:
    return f"{row['network']}.{row['station']}.{row['location']}"


def reason_class(reason: str) -> str:
    """Map a rejection reason to a class for figures and the station
    ledger. ONE definition, imported by figure.py and
    station_performance.py — the reason strings are a contract."""
    r = reason or ""
    if r.startswith("no data:"):
        return "nodata"
    if r.startswith("dead channel:"):
        return "dead"
    if r.startswith("amplitude outlier:"):
        return "amp"
    if r.startswith("not admitted:"):
        return "not_admitted"
    if r.startswith("anti-fitting:"):
        return "antifit"
    if r.startswith("aborted:"):
        return "abort"
    return "other"


def station_reads(mt, rows: list[dict]) -> dict:
    """{station_id: (own_vr, zcor_s, zcor_ok)} from a solved tensor."""
    table = {r.station: r for r in mt.station_table.itertuples()}
    out = {}
    for row in rows:
        sid = station_id(row)
        rec = table.get(sid)
        if rec is None:
            continue
        zcor = float((rec.ts - config.TIME_BEFORE_S) * config.DT)
        out[sid] = (float(rec.VR), zcor, abs(zcor) <= config.ZCOR_MAX_S)
    return out


def rank_pass1(rows: list[dict]) -> list[dict]:
    """Order the pool by pass-1 evidence: compliant time shifts first, then
    own VR. Records ``pass1.rank`` on each row."""
    ranked = sorted(
        rows,
        key=lambda r: (not r.get("pass1", {}).get("zcor_ok", True),
                       -r.get("pass1", {}).get("own_vr", -1e9)))
    for i, r in enumerate(ranked):
        r.setdefault("pass1", {})["rank"] = i
    return ranked


def keep_after_pass1(ranked: list[dict]) -> list[dict]:
    """Keep the majority: the best PASS1_KEEP_N stations that fit at all,
    then fill empty azimuth sectors up to PASS1_KEEP_MAX. Never returns
    fewer than MIN_STATIONS_USED while stations exist."""
    def _ok(r):
        p = r.get("pass1", {})
        return p.get("zcor_ok", True) and \
            p.get("own_vr", -1e9) >= config.PASS1_KEEP_VR_MIN

    keep = [r for r in ranked if _ok(r)][:config.PASS1_KEEP_N]
    if len(keep) < config.MIN_STATIONS_USED:
        for r in ranked:
            if r not in keep:
                keep.append(r)
            if len(keep) >= config.MIN_STATIONS_USED:
                break
    have = {config.sector(r["azimuth"]) for r in keep}
    for r in ranked:
        if len(keep) >= config.PASS1_KEEP_MAX:
            break
        if r in keep or not _ok(r):
            continue
        s = config.sector(r["azimuth"])
        if s not in have:
            keep.append(r)
            have.add(s)
    for r in ranked:
        r.setdefault("pass1", {})["kept"] = r in keep
    return keep


def choose_core(rows: list[dict]) -> tuple[list[dict], str]:
    """The core: the best few stations by pass-2 own VR, with the geometry
    sanity that two of them are >= AZ_PAIR_MIN_DEG apart. Returns
    (core, note); an empty core means nothing coherent to build on."""
    def _own(r):
        return r.get("pass2", {}).get("own_vr", -1e9)

    def _ok(r, floor):
        return r.get("pass2", {}).get("zcor_ok", True) and _own(r) >= floor

    note = ""
    eligible = sorted((r for r in rows if _ok(r, config.CORE_VR_MIN)),
                      key=lambda r: -_own(r))
    if len(eligible) < config.CORE_SIZE_MIN:
        eligible = sorted((r for r in rows
                           if _ok(r, config.ADMIT_VR_MIN_SPARSE)),
                          key=lambda r: -_own(r))
        note = (f"core floor relaxed to {config.ADMIT_VR_MIN_SPARSE:g} "
                "(too few stations above the normal floor)")
    if len(eligible) < config.CORE_SIZE_MIN:
        return [], "no coherent core"
    core = eligible[:config.CORE_SIZE_MAX]
    spare = [r for r in eligible if r not in core]
    swaps = 0
    while not azimuth_pair_ok([r["azimuth"] for r in core]) and swaps < 2:
        fix = next((r for r in spare
                    if azimuth_pair_ok(
                        [x["azimuth"] for x in core[:-1]] + [r["azimuth"]])),
                   None)
        if fix is None:
            break
        dropped = core[-1]
        core = core[:-1] + [fix]
        spare = [r for r in spare if r is not fix] + [dropped]
        swaps += 1
        note = (note + "; " if note else "") + \
            f"swapped {dropped['station']} for {fix['station']} to span " \
            f"{config.AZ_PAIR_MIN_DEG:g} deg of azimuth"
    return core, note


def admission_verdict(own_vr: float, zcor_s: float, zcor_ok: bool,
                      joint_before: float, joint_after: float,
                      new_sector: bool, n_current: int, tags: list[str],
                      sparse: bool) -> tuple[bool, str]:
    """Does this station earn a seat? Fit and geometry only — %DC never
    enters (a noise station can inflate it)."""
    floor = (config.ADMIT_VR_MIN_SPARSE if sparse else config.ADMIT_VR_MIN)
    if not zcor_ok:
        return False, (f"not admitted: chance alignment (zcor {zcor_s:.0f} s "
                       f"> {config.ZCOR_MAX_S:g} s)")
    if own_vr < floor and not (new_sector
                               and own_vr >= config.ADMIT_SECTOR_VR_MIN):
        return False, f"not admitted: own VR {own_vr:.0f} < {floor:g}"
    if n_current >= config.MAX_USED_STATIONS:
        return False, (f"not admitted: station cap "
                       f"{config.MAX_USED_STATIONS} reached "
                       f"(own VR {own_vr:.0f})")
    if "cluster_surplus" in (tags or []) and joint_after < joint_before:
        return False, (f"not admitted: cluster surplus with no joint gain "
                       f"({joint_before:.0f} -> {joint_after:.0f})")
    if not new_sector and \
            joint_before - joint_after > config.ADMIT_MAX_JOINT_VR_DROP:
        return False, (f"not admitted: joint VR {joint_before:.0f} -> "
                       f"{joint_after:.0f} (drop "
                       f"{joint_before - joint_after:.0f} > "
                       f"{config.ADMIT_MAX_JOINT_VR_DROP:g}), no new sector")
    bits = [f"own VR {own_vr:.0f}",
            f"joint {joint_before:.0f} -> {joint_after:.0f}"]
    if new_sector:
        bits.append("new sector")
    return True, "admitted: " + ", ".join(bits)


def invert_with_rejection(
    event: Event, stations: list[dict], depths: list[float],
    event_dir: Path, green_dir: Path,
):
    """Station selection v4 — the funnel (D. Lindsay's design, 2026-09-04).

    Nothing usable was deleted upstream: the pool arrives whole, with
    demotion TAGS (near_field / weak_signal / cluster_surplus). Selection
    is then made by the data:

    1. **Pass 1, everyone in, no time shifts.** A survey inversion over the
       whole pool on a coarse depth grid, run with correlate=0 so no trace
       can slide into a chance alignment. Each station's evidence is the
       MEDIAN of its own VR across the contiguous VR plateau (chance
       alignment is depth-specific; real coherence is not).
    2. **Keep the majority** — the best ~10 that fit at all, plus empty
       azimuth sectors.
    3. **Pass 2** re-searches with those (shifts on now, bounded by
       ZCOR_MAX_S), and the **core** is the best 3-6, forced to span
       >= 90 deg of azimuth. The core alone gets a full depth search: it is
       the reference every other station is judged against, so it must be
       clean (the audit found 69% of the old scheme's vetoes were issued by
       cores with VR < 20).
    4. **Pass 3, earn your seat back.** Every non-core station — tagged and
       pruned ones included — is added at the core depth and kept if the
       core solution predicts its waveform, or if it fills an empty azimuth
       sector and does not spoil the joint fit.
    5. **Pass 4** final full depth search, then the anti-fitting cull
       (own VR < 0 fits worse than silence) with a re-search.

    Raises NoCoherentSolution when nothing coheres. Returns
    (inversion, used_stations, rejected, selection).
    """
    rejected: list[dict] = []
    selection: dict = {"pool_n": len(stations)}

    def _solve(rows, dd, correlate=True):
        return run_inversion(write_mtinv(
            event, rows, dd, event_dir, green_dir, correlate=correlate))

    def _reject(row, reason, stage):
        row["status"] = "rejected"
        row["reason"] = reason
        row["stage"] = stage
        rejected.append({**row, "reason": reason, "stage": stage})

    def _abort(stage, reason):
        raise NoCoherentSolution(stage, reason, selection, stations, rejected)

    pool = sorted(stations, key=lambda r: r["distance_km"])
    if len(pool) < config.MIN_STATIONS_USED:
        _abort("pool", f"{len(pool)} usable stations < "
                       f"{config.MIN_STATIONS_USED}")

    # ---- pass 1: survey, everyone in, ZERO time shifts -------------------
    d1 = depths[::config.PASS1_DEPTH_STRIDE]
    if depths[-1] not in d1:
        d1 = d1 + [depths[-1]]
    inv1 = _solve(pool, d1, correlate=False)
    i1 = inv1.preferred_tensor_id
    vrs1 = [float(mt.total_VR) for mt in inv1.moment_tensors]
    plateau = plateau_indices(vrs1, i1)
    per_station: dict[str, list[float]] = {}
    for j in plateau:
        for sid, (own, _z, _ok) in station_reads(
                inv1.moment_tensors[j], pool).items():
            per_station.setdefault(sid, []).append(own)
    reads1 = station_reads(inv1.moment_tensors[i1], pool)
    for r in pool:
        sid = station_id(r)
        own_at_pref, zcor, _ = reads1.get(sid, (-1e9, 0.0, True))
        vals = per_station.get(sid, [own_at_pref])
        r["pass1"] = {
            "own_vr": round(float(np.median(vals)), 1),
            "own_vr_at_pref": round(own_at_pref, 1),
            "zcor_s": 0.0, "zcor_ok": True,  # correlate was off
        }
    vr1 = vrs1[i1]
    depth1 = float(inv1.moment_tensors[i1].depth)
    selection["pass1"] = {
        "n": len(pool), "depth_km": depth1, "vr": round(vr1, 1),
        "depths_km": [float(d) for d in d1],
    }
    print(f"pass 1: {len(pool)} stations (no shifts), depth {depth1:g} km, "
          f"VR {vr1:.1f}")

    fitting = [r for r in pool
               if r["pass1"]["own_vr"] >= config.PASS1_KEEP_VR_MIN]
    if vr1 < config.NO_SOLUTION_VR and len(fitting) < config.MIN_STATIONS_USED:
        _abort("pass1", f"nothing coheres: survey VR {vr1:.0f}, only "
                        f"{len(fitting)} stations above "
                        f"{config.PASS1_KEEP_VR_MIN:g}")

    kept = keep_after_pass1(rank_pass1(pool))
    selection["pass1"]["kept_ids"] = [station_id(r) for r in kept]
    print(f"pass 1 keeps {len(kept)}/{len(pool)}: "
          f"{[r['station'] for r in kept]}")

    # ---- pass 2: re-search the majority, shifts on -----------------------
    kept.sort(key=lambda r: r["distance_km"])
    inv2 = _solve(kept, depths)
    mt2 = inv2.moment_tensors[inv2.preferred_tensor_id]
    reads2 = station_reads(mt2, kept)
    for r in kept:
        own, zcor, ok = reads2.get(station_id(r), (-1e9, 0.0, True))
        r["pass2"] = {"own_vr": round(own, 1), "zcor_s": round(zcor, 1),
                      "zcor_ok": ok}
    vr2 = float(mt2.total_VR)
    selection["pass2"] = {"n": len(kept), "depth_km": float(mt2.depth),
                          "vr": round(vr2, 1)}
    print(f"pass 2: {len(kept)} stations, depth {float(mt2.depth):g} km, "
          f"VR {vr2:.1f}")
    if vr1 < config.NO_SOLUTION_VR and vr2 < config.NO_SOLUTION_VR:
        _abort("pass2", f"no coherent solution: survey VR {vr1:.0f}, "
                        f"majority VR {vr2:.0f}")

    core, core_note = choose_core(kept)
    if not core:
        _abort("core", f"no station fits well enough to anchor a solution "
                       f"(best own VR "
                       f"{max((r['pass2']['own_vr'] for r in kept), default=0):.0f})")
    for r in kept:
        r["pass2"]["core"] = r in core
    core.sort(key=lambda r: r["distance_km"])
    inv_c = _solve(core, depths)
    mt_c = inv_c.moment_tensors[inv_c.preferred_tensor_id]
    d_core, vr_core = float(mt_c.depth), float(mt_c.total_VR)
    selection["pass2"]["core_ids"] = [station_id(r) for r in core]
    selection["pass2"]["core_note"] = core_note
    selection["core_search"] = {"depth_km": d_core, "vr": round(vr_core, 1),
                                "n": len(core)}
    print(f"core: {len(core)} stations {[r['station'] for r in core]}, "
          f"depth {d_core:g} km, VR {vr_core:.1f}"
          + (f" ({core_note})" if core_note else ""))

    # ---- pass 3: every other station earns its seat back -----------------
    current = list(core)
    joint = vr_core
    sparse = len(core) < config.SPARSE_CORE_COUNT
    order = [r for r in rank_pass1(pool) if r not in core]
    admitted = []
    for k, cand in enumerate(order):
        trial = sorted(current + [cand], key=lambda r: r["distance_km"])
        inv_t = _solve(trial, [d_core])
        mt_t = inv_t.moment_tensors[0]
        own, zcor, zok = station_reads(mt_t, [cand]).get(
            station_id(cand), (-1e9, 0.0, False))
        joint_after = float(mt_t.total_VR)
        new_sector = config.sector(cand["azimuth"]) not in {
            config.sector(r["azimuth"]) for r in current}
        ok, reason = admission_verdict(
            own, zcor, zok, joint, joint_after, new_sector, len(current),
            cand.get("tags", []), sparse)
        cand["admission"] = {
            "order": k, "own_vr": round(own, 1), "zcor_s": round(zcor, 1),
            "zcor_ok": zok, "joint_vr_before": round(joint, 1),
            "joint_vr_after": round(joint_after, 1),
            "new_sector": new_sector, "depth_km": d_core,
            "verdict": "admitted" if ok else "rejected", "reason": reason,
        }
        if ok:
            current = trial
            joint = joint_after
            admitted.append(cand)
        else:
            _reject(cand, reason, "admission")
    selection["admission_order"] = [station_id(r) for r in order]
    print(f"pass 3: {len(admitted)}/{len(order)} earned a seat "
          f"({[r['station'] for r in admitted]})")

    # ---- pass 4: final search + anti-fitting cull ------------------------
    inv4 = None
    for _ in range(3):
        current.sort(key=lambda r: r["distance_km"])
        inv4 = _solve(current, depths)
        mt4 = inv4.moment_tensors[inv4.preferred_tensor_id]
        reads4 = station_reads(mt4, current)
        culled = [r for r in current
                  if reads4.get(station_id(r), (0.0, 0.0, True))[0]
                  < config.ANTIFIT_VR]
        if not culled:
            break
        for r in culled:
            own = reads4[station_id(r)][0]
            _reject(r, f"anti-fitting: own VR {own:.0f} at "
                       f"{float(mt4.depth):g} km", "cull")
        current = [r for r in current if r not in culled]
        print(f"anti-fitting cull removed {len(culled)}: "
              f"{[r['station'] for r in culled]}")
        if len(current) < config.MIN_STATIONS_USED:
            _abort("cull", f"fewer than {config.MIN_STATIONS_USED} coherent "
                           "stations after the anti-fitting cull")

    mt4 = inv4.moment_tensors[inv4.preferred_tensor_id]
    vr4 = float(mt4.total_VR)
    if vr4 < config.NO_SOLUTION_VR:
        _abort("final", f"final VR {vr4:.0f} < {config.NO_SOLUTION_VR:g}")
    reads4 = station_reads(mt4, current)
    for r in current:
        own, zcor, _ = reads4.get(station_id(r), (0.0, 0.0, True))
        r["status"] = "used"
        r["final"] = {"own_vr": round(own, 1), "zcor_s": round(zcor, 1)}
    selection["final"] = {"depth_km": float(mt4.depth), "vr": round(vr4, 1),
                          "n_used": len(current)}
    selection["aborted"] = False
    return inv4, current, rejected, selection


def no_solution_record(event: Event, pool: list[dict], dropped: list[dict],
                       model: str, band: tuple[float, float], stage: str,
                       reason: str, selection: dict) -> dict:
    """Archive record for an event with no coherent solution: the honest
    alternative to publishing a junk mechanism. Carries the full station
    ledger and every pass's evidence so the call can be audited."""
    import mttime
    import obspy

    best_vr = max(
        [selection.get(k, {}).get("vr", 0.0)
         for k in ("pass1", "pass2", "core_search", "final")] or [0.0])
    return {
        "status": config.STATUS_NO_SOLUTION,
        "event": event.to_dict(),
        "abort": {"stage": stage, "reason": reason,
                  "best_vr": round(float(best_vr), 1),
                  "band": config.band_tag(band)},
        "selection": selection,
        "stations_used": [],
        "stations_dropped": dropped + [
            {**r, "reason": f"aborted: {reason}", "stage": "abort"}
            for r in pool if r.get("status") != "rejected"
        ] + [r for r in pool if r.get("status") == "rejected"],
        "filter_band_hz": list(band),
        "chosen_band": config.band_tag(band),
        "quality": {"grade": "X", "passed": False, "checks": {},
                    "warnings": {"no_coherent_solution": True},
                    "n_stations_used": 0, "azimuthal_gap_deg": 360.0},
        "publish_decision": {"publish": False,
                             "reason_tags": ["no_coherent_solution"],
                             "reasons": [reason]},
        "provenance": {
            "velocity_model": model,
            "gf_version": config.GF_VERSION,
            "selection_version": config.SELECTION_VERSION,
            "code_commit": config.code_version(),
            "mttime_version": mttime.__version__,
            "obspy_version": obspy.__version__,
        },
    }


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
              model: str, selection: dict | None = None) -> dict:
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
        vrs = [r["vr"] for r in srt]
        i_max = max(range(len(srt)), key=lambda i: vrs[i])
        span = plateau_indices(vrs, i_max)
        return srt[span[-1]]["depth_km"] - srt[span[0]]["depth_km"]

    best_vr = max(rows, key=lambda r: r["vr"])
    best_dc = max(rows, key=lambda r: r["pdc"])
    pref = inv.moment_tensors[inv.preferred_tensor_id]

    edge_depths = (min(config.GF_DEPTHS_KM), max(config.GF_DEPTHS_KM))
    at_edge = float(pref.depth) in edge_depths
    plateau_km = round(_plateau_span(rows), 1)
    _, guard_applied = pick_preferred_detail(
        [(r["vr"], r["pdc"]) for r in rows], contiguous=len(rows) > 1)

    solution = {
        "status": config.STATUS_SOLVED,
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
            "plateau_km": plateau_km,
            "depth_unconstrained": plateau_km > 10.0,
            # a single-point VR maximum at the end of the grid is an
            # artifact, not a depth (see grid_edge_reference)
            "edge_artifact": bool(at_edge and plateau_km == 0.0),
            "grid_edge_guard_applied": bool(guard_applied),
        },
        "selection": selection or {},
        "stations_used": stations,
        "stations_dropped": dropped,
        "provenance": {
            "velocity_model": model,
            "gf_version": config.GF_VERSION,
            "selection_version": config.SELECTION_VERSION,
            "code_commit": config.code_version(),
            "mttime_version": mttime.__version__,
            "obspy_version": obspy.__version__,
            "degree": config.INVERSION_DEGREE,
            "weight": "distance",
            "correlate": True,
            "npts": config.INV_NPTS,
            "dt_s": config.DT,
            "preferred_rule": (
                f"max %DC on the contiguous depth plateau within "
                f"{config.PREFER_DC_VR_TOLERANCE:g} VR points of the "
                f"maximum, with the grid-edge guard"
            ),
        },
    }
    solution["quality"] = quality_gates(solution)
    return solution


def depth_plausible(our_depth_km: float, geonet_depth_km: float) -> bool:
    """Does our centroid depth agree with the GeoNet hypocentre?

    The search itself stays fully independent — this only judges the
    RESULT. A GeoNet placeholder depth (5/12/33 km) is not a measurement,
    so nothing is judged against it. Threshold: DEPTH_PLAUSIBLE_MAX_KM
    (8 km), close to the typical agreement in the Ristau catalogue, whose
    centroid depths sit a median 4 km from the GeoNet hypocentre."""
    if geonet_depth_km in config.PLACEHOLDER_DEPTHS_KM:
        return True
    return abs(our_depth_km - geonet_depth_km) <= config.DEPTH_PLAUSIBLE_MAX_KM


def grade_v2(vr: float, dc: float, min_own_vr: float,
             jk_rot: float | None, edge_artifact: bool,
             pair_ok: bool, depth_ok: bool = True) -> str:
    """Letter grade from evidence only (v2, 2026-09-04).

    Station COUNT and azimuthal GAP are deliberately absent: three
    well-fitting stations spanning 90 degrees make a good solution (BSL
    practice), while ten stations carrying a passenger do not. What counts:
    how much of the data the solution explains (vr), whether EVERY used
    station fits (min_own_vr — no passengers), whether the mechanism
    survives leaving one out (jk_rot), whether the depth is real
    (edge_artifact) and whether the geometry can resolve a mechanism at all
    (pair_ok). B additionally requires DC >= 60, making it exactly BSL's
    publishability rule.
    """
    ra, rb, rc = (config.GRADE_RUBRIC[k] for k in ("A", "B", "C"))
    if not pair_ok or vr < rc["vr"] or min_own_vr < rc["min_own_vr"]:
        return "D"
    if (edge_artifact or not depth_ok or vr < rb["vr"] or dc < rb["dc"]
            or min_own_vr < rb["min_own_vr"]
            or (jk_rot is not None and jk_rot > rb["jk_rot_max"])):
        return "C"
    if (vr < ra["vr"] or dc < ra["dc"] or min_own_vr < ra["min_own_vr"]
            or jk_rot is None or jk_rot > ra["jk_rot_max"]):
        return "B"
    return "A"


def quality_gates(solution: dict) -> dict:
    """Publication gates — fail loud, publish nothing below the bar.

    Idempotent: called once inside summarize (before the jackknife exists)
    and again after the jackknife is attached, so the stability evidence
    reaches the grade."""
    pref = solution["preferred"]
    used = solution["stations_used"]
    n_used = len(used)
    azimuths = sorted(s["azimuth"] for s in used)
    gaps = [
        (azimuths[(i + 1) % len(azimuths)] - a) % 360.0
        for i, a in enumerate(azimuths)
    ]
    az_gap = max(gaps) if gaps else 360.0
    pair_ok = azimuth_pair_ok(azimuths)

    own_vrs = [s.get("final", {}).get("own_vr", s.get("station_vr"))
               for s in used]
    own_vrs = [v for v in own_vrs if v is not None]
    min_own_vr = min(own_vrs) if own_vrs else -999.0

    jk = solution.get("jackknife", {}) or {}
    jk_rot = (jk.get("max_tensor_rotation_deg")
              if jk.get("n_subsets") else None)
    flags = solution.get("depth_pick_flags", {})
    edge_artifact = bool(flags.get("edge_artifact", False))

    depth_ok = depth_plausible(pref["depth_km"],
                               solution["event"].get("depth_km", 0.0))
    grade = grade_v2(pref["vr"], pref["pdc"], min_own_vr, jk_rot,
                     edge_artifact, pair_ok, depth_ok)
    rb = config.GRADE_RUBRIC["B"]
    checks = {
        "vr_floor": pref["vr"] >= config.GRADE_RUBRIC["C"]["vr"],
        "no_passengers": min_own_vr >= rb["min_own_vr"],
        "az_pair_90": pair_ok,
        "jackknife_stable": jk_rot is None or jk_rot <= rb["jk_rot_max"],
        "depth_interior": not edge_artifact,
        "depth_agrees_with_geonet": depth_ok,
        "dc_floor": pref["pdc"] >= rb["dc"],
    }
    warnings = {
        "depth_at_grid_edge": pref["depth_km"] in (
            min(config.GF_DEPTHS_KM), max(config.GF_DEPTHS_KM),
        ),
        "grid_edge_guard_applied": bool(
            flags.get("grid_edge_guard_applied", False)),
        "jackknife_skipped": jk_rot is None,
        "depth_far_from_geonet": not depth_ok,
    }
    return {
        # informational, no longer thresholds
        "n_stations_used": n_used,
        "azimuthal_gap_deg": round(az_gap, 1),
        "min_own_vr": round(min_own_vr, 1),
        "jackknife_rotation_deg": jk_rot,
        "depth_minus_geonet_km": round(
            pref["depth_km"] - solution["event"].get("depth_km", 0.0), 1),
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
