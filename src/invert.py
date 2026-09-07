"""Task 3 — moment tensor inversion.

    pixi run python invert.py --event 2026p666955      # run this task alone
    (needs task 2's SAC files in <events>/<id>/task2; writes there too)

The inversion itself is mttime (Chiang, LLNL; github.com/LLNL/mttime),
the Python implementation of the time-domain method of Dreger &
Helmberger (1993). This module only writes mttime's control file, runs
it, and applies three published rules around it:

  1. DEPTH: the search grid is the GeoNet hypocentre +/- depthWindowKm
     (F-net searches JMA +/-30 km, Gisola +/-31 km, W-phase and SCARDEC
     +/-50 km); the preferred depth is mttime's own maximum-VR pick, and
     the depth uncertainty is the range within depthUncPct of that
     maximum (Vallee et al. 2011; Bernardi et al. 2004).
  2. STATIONS: the loop of Clinton, Hauksson & Solanki (2006, BSSA 96) —
     invert everyone; drop any station whose own VR is below
     max(overall VR - stationVRDrop, stationVRFloor) or whose solved
     time shift exceeds maxTimeShiftS (NEIC SynDepth clips at 6-10 s);
     re-invert; stop when nothing drops or minStations remain.
  3. GRADE: INGV's table (terremoti.ingv.it/en/help#TDMT), where the VR
     needed for a grade FALLS as the station count rises, so dropping
     stations cannot buy a grade (Triantafyllis et al. 2016 show 2
     stations scoring VR 0.9 with a condition number > 10). A and B also
     need DC >= dcMinPublish (the BSL rule) and jackknife rotation <=
     jackknifeMaxDeg (Fukuyama et al. 1998: the jackknife is the only
     reliable detector of a station gone bad).

Mechanism comparison uses the minimum rotation angle (Townend et al.
2012; Walsh, Arnold & Townend 2009; Kagan 1991). All numbers are in
auto_tdmt.cfg §3.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

import config
import greens
from config import P
from geonet import Event


def station_id(row: dict) -> str:
    return f"{row['network']}.{row['station']}.{row['location']}"


# ---------------------------------------------------------------------------
# 3a. depth grid
# ---------------------------------------------------------------------------

def search_depths(event: Event, model: str) -> list[float]:
    """Library depths inside GeoNet depth +/- depthWindowKm; the full
    library when GeoNet's depth is a fixed placeholder (not a
    measurement). The search is bounded, the answer is never forced."""
    lib = greens.available_depths(model)
    assert lib, "GF library is empty — build it first (greens.py --build)"
    if event.depth_km in config.PLACEHOLDER_DEPTHS_KM:
        return lib
    w = P.invert.depthWindowKm
    inside = [d for d in lib if abs(d - event.depth_km) <= w]
    return inside if len(inside) >= 3 else lib


def depth_uncertainty(rows: list[dict], pct: float | None = None
                      ) -> tuple[float, float]:
    """Range of depths whose VR is within ``pct`` percent of the maximum
    (SCARDEC: Vallee et al. 2011; Bernardi et al. 2004), as (lo, hi) km."""
    pct = P.invert.depthUncPct if pct is None else pct
    vmax = max(r["vr"] for r in rows)
    ok = [r["depth_km"] for r in rows if r["vr"] >= vmax * (1 - pct / 100.0)]
    return (min(ok), max(ok))


# ---------------------------------------------------------------------------
# 3b. mttime
# ---------------------------------------------------------------------------

def write_mtinv(event: Event, stations: list[dict], depths: list[float],
                event_dir: Path, green_dir: Path,
                correlate: bool | None = None) -> Path:
    """Write mttime's control file (exactly the mttime example format;
    mttime reads it by LINE ORDER, 12 header lines then the table)."""
    correlate = bool(P.invert.correlate) if correlate is None else correlate
    npts_col = []
    for r in stations:
        end = r.get("window_end_s")
        if end is None:  # rows without task-2 metadata: kinematic window
            end = (r["distance_km"] / P.station.groupVelKms
                   + config.window_tail_s(event.prelim_mag))
        npts_col.append(int(min(config.INV_NPTS, max(
            P.station.windowMinS, config.TIME_BEFORE_S + end))))
    df = pd.DataFrame({
        "station": [station_id(r) for r in stations],
        "distance": [r["gf_distance_km"] for r in stations],
        "azimuth": [round(r["azimuth"], 2) for r in stations],
        "ts": config.TIME_BEFORE_S,
        "npts": npts_col,
        "dt": config.DT,
        "used": 1,
        "longitude": [r["longitude"] for r in stations],
        "latitude": [r["latitude"] for r in stations],
    })
    headers = dict(
        datetime=event.origin_time, longitude=event.longitude,
        latitude=event.latitude,
        depth=",".join(f"{d:.4f}" for d in depths),
        path_to_data=".", path_to_green=green_dir.name, green="herrmann",
        components=P.invert.components, degree=P.invert.degree,
        weight=P.invert.weight, plot=0, correlate=1 if correlate else 0,
    )
    path = event_dir / "mtinv.in"
    with open(path, "w") as f:
        for key, value in headers.items():
            f.write(f"{key:<15}{value}\n")
        f.write(df.to_string(index=False))
    return path


def run_inversion(mtinv_path: Path):
    """mttime Configure -> Inversion -> invert(). The preferred tensor is
    mttime's own: the maximum total VR over the depth list."""
    from mttime import Configure, Inversion
    cfg = Configure(path_to_file=str(mtinv_path))
    inv = Inversion(config=cfg)
    inv.invert(show=False)
    assert inv.moment_tensors, "mttime produced no solutions"
    return inv


def station_reads(mt, rows: list[dict]) -> dict[str, tuple[float, float]]:
    """{station id: (own VR, time shift s)} from one mttime Tensor."""
    out = {}
    table = {r.station: r for r in mt.station_table.itertuples()}
    for r in rows:
        t = table.get(station_id(r))
        if t is not None:
            out[station_id(r)] = (float(t.VR),
                                  float((t.ts - config.TIME_BEFORE_S) * config.DT))
    return out


# ---------------------------------------------------------------------------
# 3c. the Clinton loop
# ---------------------------------------------------------------------------

class NoCoherentSolution(Exception):
    def __init__(self, stage: str, reason: str, rejected: list[dict],
                 best_vr: float = 0.0):
        super().__init__(f"{stage}: {reason}")
        self.stage, self.reason = stage, reason
        self.rejected, self.best_vr = rejected, best_vr


def station_verdict(own_vr: float, shift_s: float, overall_vr: float
                    ) -> str | None:
    """None to keep, else the reason to drop.

    The floor is a FIXED own-VR bar (gempa's minimumFinalStationFit) plus
    the time-shift cap. Clinton et al. (2006) also drop stations below
    (overall VR - 10), but they then RESELECT new stations to refill the
    azimuth sectors; without that refill the relative rule ratchets —
    every drop raises the overall VR and so the bar — and the smoke test
    took a 16-station M4.9 to 3 stations in three rounds. The relative
    rule is therefore off unless stationVRDrop > 0 in the cfg."""
    floor = P.invert.stationVRFloor
    if P.invert.stationVRDrop > 0:
        floor = max(overall_vr - P.invert.stationVRDrop, floor)
    if abs(shift_s) > P.invert.maxTimeShiftS:
        return (f"time shift: {shift_s:+g} s exceeds "
                f"{P.invert.maxTimeShiftS:g} s (chance alignment)")
    if own_vr < floor:
        return f"poor fit: own VR {own_vr:.0f} < {floor:.0f}"
    return None


def clinton_loop(event: Event, stations: list[dict], depths: list[float],
                 event_dir: Path, green_dir: Path):
    """Invert, drop the stations the solution cannot explain, repeat.

    Returns (inversion, used, rejected, rounds). Raises NoCoherentSolution
    when fewer than minStations remain.
    """
    current = [dict(r) for r in stations]
    rejected: list[dict] = []
    rounds: list[dict] = []
    inv = None
    at_floor = False
    for n in range(1, P.invert.maxRounds + 2):
        if len(current) < P.invert.minStations:
            raise NoCoherentSolution(
                f"round {n}", f"only {len(current)} stations left "
                f"(minimum {P.invert.minStations})", rejected,
                rounds[-1]["vr"] if rounds else 0.0)
        inv = run_inversion(write_mtinv(event, current, depths, event_dir,
                                        green_dir))
        pref = inv.moment_tensors[inv.preferred_tensor_id]
        reads = station_reads(pref, current)
        overall = float(pref.total_VR)
        drops = []
        for r in current:
            own, shift = reads[station_id(r)]
            r["station_vr"], r["zcor_s"] = round(own, 1), round(shift, 1)
            why = station_verdict(own, shift, overall)
            if why:
                drops.append((r, why))
        rounds.append({"round": n, "n_stations": len(current),
                       "depth_km": float(pref.depth), "vr": round(overall, 1),
                       "dropped": [station_id(r) for r, _ in drops]})
        print(f"round {n}: {len(current)} stations, depth {pref.depth:g} km, "
              f"VR {overall:.1f}" + (f", dropping {len(drops)}" if drops else
                                    " — stable"))
        if not drops or at_floor or n > P.invert.maxRounds:
            if drops and at_floor:
                print(f"  at the {P.invert.minStations}-station floor: "
                      f"keeping the rest")
            elif drops:
                print(f"  stopped after {P.invert.maxRounds} rounds (maxRounds)")
            break
        # Backward elimination: chance alignments go at once (a property
        # of the trace, not of the reference), poor fits go at most
        # maxDropPerRound per round, worst first, so the reference
        # solution improves before the marginal stations are judged
        # ("every judgment is only as good as its reference solution":
        # 2026p669681 lost a station at own VR 24 vs a floor of 25 when
        # 16 stations were judged against a joint fit of VR 1.6).
        # Never below the minimum in one go.
        shift_drops = [d for d in drops if d[1].startswith("time shift")]
        fit_drops = sorted([d for d in drops if d not in shift_drops],
                           key=lambda d: d[0]["station_vr"])
        drops = shift_drops + fit_drops[:P.invert.maxDropPerRound]
        keep_n = max(P.invert.minStations, len(current) - len(drops))
        n_drop = len(current) - keep_n
        for r, why in drops[:n_drop]:
            print(f"  drop {station_id(r)}: {why}")
            rejected.append({**r, "station": station_id(r), "reason": why,
                             "round": n})
        dropped_ids = {station_id(r) for r, _ in drops[:n_drop]}
        current = [r for r in current if station_id(r) not in dropped_ids]
        at_floor = n_drop < len(drops)   # one more inversion, then stop
        if len(fit_drops) > P.invert.maxDropPerRound and not at_floor:
            print(f"  ({len(fit_drops) - P.invert.maxDropPerRound} more "
                  f"poor fits deferred to the next round)")
    if inv is None or len(current) < P.invert.minStations:
        raise NoCoherentSolution("final", f"{len(current)} stations",
                                 rejected, rounds[-1]["vr"] if rounds else 0.0)
    # final numbers for the surviving set (last inversion == current set)
    pref = inv.moment_tensors[inv.preferred_tensor_id]
    for r, (own, shift) in zip(current, (station_reads(pref, current)[
            station_id(r)] for r in current)):
        r["station_vr"], r["zcor_s"] = round(own, 1), round(shift, 1)
    return inv, current, rejected, rounds


def reason_class(reason: str) -> str:
    """Shared vocabulary for drop reasons -> class (figures, ledger)."""
    r = (reason or "").lower()
    for prefix, cls in (("no data", "nodata"), ("low snr", "snr"),
                        ("amplitude outlier", "amp"),
                        ("not selected", "not_selected"),
                        ("poor fit", "fit"), ("time shift", "shift"),
                        ("aborted", "abort")):
        if r.startswith(prefix):
            return cls
    return "other"


# ---------------------------------------------------------------------------
# 3d. mechanism comparison (Townend recipe) and the jackknife
# ---------------------------------------------------------------------------

def _mech_rotation_matrix(strike: float, dip: float, rake: float) -> np.ndarray:
    """Mechanism as a rotation matrix in NED (Walsh et al. 2009 sense):
    columns are the T, P and null axes from the Aki & Richards (1980)
    slip vector u and fault normal n. In this frame the double-couple
    symmetry group is the diagonal sign flips, so the minimum below is
    independent of which nodal plane parameterises either mechanism."""
    phi, delta, lam = np.radians([strike, dip, rake])
    n = np.array([-np.sin(delta) * np.sin(phi), np.sin(delta) * np.cos(phi),
                  -np.cos(delta)])
    u = np.array([
        np.cos(lam) * np.cos(phi) + np.cos(delta) * np.sin(lam) * np.sin(phi),
        np.cos(lam) * np.sin(phi) - np.cos(delta) * np.sin(lam) * np.cos(phi),
        -np.sin(lam) * np.sin(delta)])
    t, p = (u + n) / np.sqrt(2.0), (u - n) / np.sqrt(2.0)
    return np.column_stack([t, p, np.cross(t, p)])


def min_rotation_angle_deg(a, b) -> float:
    """Minimum rotation aligning two double couples (Townend et al. 2012
    supplement eq. 1; Walsh, Arnold & Townend 2009; cf. Kagan 1991):
    arccos((tr(R1^T R2 S) - 1) / 2) minimised over the DC symmetry group.
    ``a``/``b`` are (strike, dip, rake)."""
    r1, r2 = _mech_rotation_matrix(*a), _mech_rotation_matrix(*b)
    best = 180.0
    for sym in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)):
        cos = (np.trace(r1.T @ (r2 * np.array(sym))) - 1.0) / 2.0
        best = min(best, float(np.degrees(np.arccos(np.clip(cos, -1, 1)))))
    return best


def jackknife(event: Event, stations: list[dict], depth: float,
              band_dir: Path, green_dir: Path, full_plane1: dict) -> dict:
    """Leave-one-station-out at the preferred depth. Stability = the
    largest minimum-rotation angle of any subset mechanism from the full
    one. Restores the full mtinv.in afterwards."""
    if len(stations) < 4 or not P.invert.jackknife:
        return {"n_subsets": 0, "note": "fewer than 4 stations or disabled"}
    ref = (full_plane1["strike"], full_plane1["dip"], full_plane1["rake"])
    subsets = []
    try:
        for i in range(len(stations)):
            subset = stations[:i] + stations[i + 1:]
            mt = run_inversion(write_mtinv(event, subset, [depth], band_dir,
                                           green_dir)).moment_tensors[0]
            fps = np.asarray(mt.fps, dtype=float)
            subsets.append({
                "left_out": stations[i]["station"],
                "mw": round(float(mt.mw), 3), "pdc": round(float(mt.pdc), 1),
                "vr": round(float(mt.total_VR), 1),
                "tensor_rotation_deg": round(
                    min_rotation_angle_deg(ref, tuple(fps[0])), 1),
                "plane1": dict(zip(("strike", "dip", "rake"),
                                   [round(float(v), 1) for v in fps[0]])),
            })
    finally:
        write_mtinv(event, stations, [depth], band_dir, green_dir)
    rots = [s["tensor_rotation_deg"] for s in subsets]
    return {
        "n_subsets": len(subsets),
        "mw_std": round(float(np.std([s["mw"] for s in subsets])), 3),
        "dc_std": round(float(np.std([s["pdc"] for s in subsets])), 1),
        "max_tensor_rotation_deg": round(max(rots), 1),
        "mean_tensor_rotation_deg": round(float(np.mean(rots)), 1),
        "subsets": subsets,
    }


# ---------------------------------------------------------------------------
# 3e. grade and record
# ---------------------------------------------------------------------------

def grade_ingv(vr: float, n_stations: int, dc: float,
               jk_rot: float | None) -> str:
    """INGV letter from (VR, N stations), then the BSL DC rule and the
    jackknife for A/B. Below the C bar is D."""
    col = P.grade_vr_column(n_stations)
    a, b, c = (P.invert.gradeA_VR[col], P.invert.gradeB_VR[col],
               P.invert.gradeC_VR[col])
    if vr < c:
        return "D"
    if vr < b:
        return "C"
    publishable = (dc >= P.invert.dcMinPublish
                   and (jk_rot is None or jk_rot <= P.invert.jackknifeMaxDeg))
    if not publishable:
        return "C"
    return "A" if vr >= a else "B"


def quality_gates(solution: dict) -> dict:
    """Grade + the checks behind it (same shape for every consumer)."""
    pref = solution["preferred"]
    used = solution["stations_used"]
    n = len(used)
    az = sorted(s["azimuth"] for s in used)
    gaps = [(az[(i + 1) % len(az)] - a) % 360.0 for i, a in enumerate(az)]
    own = [s.get("station_vr") for s in used if s.get("station_vr") is not None]
    jk = solution.get("jackknife") or {}
    jk_rot = jk.get("max_tensor_rotation_deg") if jk.get("n_subsets") else None
    grade = grade_ingv(pref["vr"], n, pref["pdc"], jk_rot)
    col = P.grade_vr_column(n)
    lo, hi = solution.get("depth_pick_flags", {}).get("depth_range_km",
                                                        (None, None))
    depths = [r["depth_km"] for r in solution.get("depth_search", [])]
    at_edge = bool(depths) and pref["depth_km"] in (min(depths), max(depths))
    checks = {
        "vr_grade": pref["vr"] >= P.invert.gradeB_VR[col],
        "dc_floor": pref["pdc"] >= P.invert.dcMinPublish,
        "jackknife_stable": jk_rot is None or jk_rot <= P.invert.jackknifeMaxDeg,
        "min_stations": n >= P.invert.minStations,
    }
    return {
        "n_stations_used": n,
        "azimuthal_gap_deg": round(max(gaps), 1) if gaps else 360.0,
        "min_own_vr": round(min(own), 1) if own else None,
        "jackknife_rotation_deg": jk_rot,
        "depth_range_km": [lo, hi],
        "depth_minus_geonet_km": round(
            pref["depth_km"] - solution["event"].get("depth_km", 0.0), 1),
        "grade": grade,
        "checks": checks,
        "warnings": {"depth_at_window_edge": at_edge,
                     "jackknife_skipped": jk_rot is None},
        "passed": grade <= P.publish.minGrade,   # "A" <= "B"
    }


def _provenance(model: str) -> dict:
    import mttime
    import obspy
    return {
        "velocity_model": model, "gf_version": config.GF_VERSION,
        "selection_version": config.SELECTION_VERSION,
        "code_commit": config.code_version(),
        "mttime_version": mttime.__version__, "obspy_version": obspy.__version__,
        "params": P.as_dict(), "params_file": P.source,
        "preferred_rule": "maximum VR over the depth grid (mttime argmax)",
    }


def summarize(inv, event: Event, stations: list[dict], dropped: list[dict],
              model: str, rounds: list[dict] | None = None) -> dict:
    """Serialize the depth search + preferred solution with provenance."""
    def _mt_row(mt) -> dict:
        fps = np.asarray(mt.fps, dtype=float)
        assert fps.shape == (2, 3), f"unexpected fps shape {fps.shape}"
        return dict(
            depth_km=float(mt.depth), mw=float(mt.mw), m0_dyne_cm=float(mt.mo),
            vr=float(mt.total_VR), pdc=float(mt.pdc), pclvd=float(mt.pclvd),
            piso=float(mt.piso),
            plane1=dict(zip(("strike", "dip", "rake"), fps[0])),
            plane2=dict(zip(("strike", "dip", "rake"), fps[1])),
            tensor_rtp_dyne_cm={k: float(v) for k, v in
                                mt.get_tensor_elements(basis="RTP").items()},
        )

    rows = [_mt_row(mt) for mt in inv.moment_tensors]
    pref = inv.moment_tensors[inv.preferred_tensor_id]
    lo, hi = depth_uncertainty(rows)
    best_dc = max(rows, key=lambda r: r["pdc"])
    solution = {
        "status": config.STATUS_SOLVED,
        "event": event.to_dict(),
        "preferred": {
            **_mt_row(pref),
            "tensor_dyne_cm": {k: float(v) for k, v in
                               pref.get_tensor_elements().items()},
            "tensor_rtp_dyne_cm": {k: float(v) for k, v in
                                   pref.get_tensor_elements(basis="RTP").items()},
        },
        "depth_search": rows,
        "depth_pick_flags": {
            "vr_max_depth_km": float(pref.depth),
            "dc_max_depth_km": best_dc["depth_km"],
            "depth_range_km": [lo, hi],
            "plateau_km": round(hi - lo, 1),
            "depth_unconstrained": (hi - lo) > 10.0,
        },
        "selection": {"rounds": rounds or []},
        "stations_used": stations,
        "stations_dropped": dropped,
        "provenance": _provenance(model),
    }
    solution["quality"] = quality_gates(solution)
    return solution


def no_solution_record(event: Event, pool: list[dict], dropped: list[dict],
                       model: str, band: tuple[float, float], stage: str,
                       reason: str, best_vr: float = 0.0) -> dict:
    """The record for an event the data could not constrain: the full
    station ledger, no mechanism (F-net likewise publishes nothing below
    its station floor, Fukuyama et al. 1998)."""
    return {
        "status": config.STATUS_NO_SOLUTION,
        "event": event.to_dict(),
        "abort": {"stage": stage, "reason": reason,
                  "best_vr": round(float(best_vr), 1),
                  "band": config.band_tag(band)},
        "stations_used": [],
        "stations_dropped": dropped + [
            {**r, "station": station_id(r), "reason": f"aborted: {reason}"}
            for r in pool],
        "filter_band_hz": list(band), "chosen_band": config.band_tag(band),
        "quality": {"grade": "X", "passed": False, "checks": {},
                    "warnings": {"no_coherent_solution": True},
                    "n_stations_used": 0, "azimuthal_gap_deg": 360.0},
        "publish_decision": {"publish": False,
                             "reason_tags": ["no_coherent_solution"],
                             "reasons": [reason]},
        "provenance": _provenance(model),
    }


def save_solution(solution: dict, event_dir: Path) -> Path:
    out = event_dir / "solution.json"
    with open(out, "w") as f:
        json.dump(solution, f, indent=2)
    return out


# ---------------------------------------------------------------------------
# standalone
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--event", required=True, help="GeoNet publicID")
    ap.add_argument("--workdir", default=None,
                    help="directory with task-2 SAC files (default "
                         "<events>/<id>/task2)")
    args = ap.parse_args()
    import waveforms
    from geonet import get_event
    ev = get_event(args.event)
    wd = Path(args.workdir) if args.workdir else \
        config.EVENTS_DIR / ev.public_id / "task2"
    model = config.model_for_event(ev.latitude, ev.longitude)
    band = config.band_candidates(ev.prelim_mag)[0]
    pool, dropped = waveforms.fetch_and_process(ev, wd, band)
    depths = search_depths(ev, model)
    print(f"depth grid {depths[0]:g}-{depths[-1]:g} km ({len(depths)}), "
          f"model {model}, {len(pool)} stations in")
    greens.stage_event_greens(model, pool, depths, band, wd / "greens")
    cwd = os.getcwd()
    os.chdir(wd)
    try:
        inv, used, rejected, rounds = clinton_loop(ev, pool, depths, wd,
                                                   wd / "greens")
        sol = summarize(inv, ev, used, dropped + rejected, model, rounds)
    finally:
        os.chdir(cwd)
    p, q = sol["preferred"], sol["quality"]
    print(f"\nMw {p['mw']:.2f}  depth {p['depth_km']:g} km "
          f"[{q['depth_range_km'][0]:g}-{q['depth_range_km'][1]:g}]  "
          f"VR {p['vr']:.1f}  DC {p['pdc']:.0f}  grade {q['grade']}  "
          f"{q['n_stations_used']} stations")
    print(f"plane1 {p['plane1']}\nplane2 {p['plane2']}")
    save_solution(sol, wd)
