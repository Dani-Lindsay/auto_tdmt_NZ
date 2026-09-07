"""Process one GeoNet event end to end — the five tasks in sequence.

    pixi run python src/process_event.py --event 2026p660242 --debug
    pixi run python src/process_event.py --event 2026p660242 --band 0.02-0.10

  task 1  event metadata                 geonet.get_event
  task 2  stations + waveforms           waveforms.fetch_and_process
  task 3  inversion (mttime)             invert.clinton_loop -> summarize -> jackknife
  task 4  forward model (Okada)          okada_forward.forward_both_planes
  task 5  archive + publish decision     catalogue / trigger / publish / figure

The filter-band menu (auto_tdmt.cfg §2) is an ORDERED preference: the
first band whose solution passes its gate wins; later bands run only when
an earlier one fails. Each task can also be run on its own (see the
module docstrings) and is walked through in docs/task_*.ipynb.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path

import catalogue
import config
import diagnostics
import figure
import greens
import invert
import nisar_dates
import okada_forward
import publish
import trigger
import waveforms
from geonet import get_event


def process_band(event, band: tuple[float, float], band_dir: Path,
                 debug: bool, model: str) -> dict:
    """Tasks 2 and 3 for one filter band inside band_dir."""
    band_dir.mkdir(parents=True, exist_ok=True)
    green_dir = band_dir / "greens"
    print(f"\n=== band {1/band[1]:.0f}-{1/band[0]:.0f} s "
          f"({band[0]:g}-{band[1]:g} Hz) ===")

    # --- task 2 -------------------------------------------------------------
    stages: dict | None = {} if debug else None
    pool, dropped = waveforms.fetch_and_process(event, band_dir, band,
                                                stages=stages)
    print(f"task 2: {len(pool)} stations selected, {len(dropped)} not")
    for d in dropped:
        print(f"  {d['station']:14s} {d['reason']}")
    if debug and pool:
        diag_dir = band_dir / "diagnostics"
        figs = diagnostics.plot_stages(stages, event.to_dict(), band, diag_dir)
        figs.append(diagnostics.plot_station_map(
            event.to_dict(), pool, [d["station"] for d in dropped], diag_dir))
        for f in figs:
            print(f"  diagnostic: {f}")
    if len(pool) < config.MIN_STATIONS_USED:
        solution = invert.no_solution_record(
            event, pool, dropped, model, band, "task 2",
            f"only {len(pool)} usable stations (minimum "
            f"{config.MIN_STATIONS_USED})")
        invert.save_solution(solution, band_dir)
        print("NO COHERENT SOLUTION: too few usable stations")
        return solution

    # --- task 3 -------------------------------------------------------------
    depths = invert.search_depths(event, model)
    print(f"task 3: depth grid {depths[0]:g}-{depths[-1]:g} km "
          f"({len(depths)} depths)")
    greens.stage_event_greens(model, pool, depths, band, green_dir)
    cwd = os.getcwd()
    os.chdir(band_dir)
    try:
        inv, used, rejected, rounds = invert.clinton_loop(
            event, pool, depths, band_dir, green_dir)
        inv.plot(view="waveform", option="preferred", format="jpg", show=False)
    except invert.NoCoherentSolution as e:
        solution = invert.no_solution_record(
            event, pool, dropped, model, band, e.stage, e.reason, e.best_vr)
        solution["stations_dropped"] += e.rejected
        print(f"NO COHERENT SOLUTION ({e.stage}): {e.reason}")
    else:
        solution = invert.summarize(inv, event, used, dropped + rejected,
                                    model, rounds)
        solution["filter_band_hz"] = list(band)
        jk = invert.jackknife(event, used, solution["preferred"]["depth_km"],
                              band_dir, green_dir,
                              solution["preferred"]["plane1"])
        solution["jackknife"] = jk
        solution["quality"] = invert.quality_gates(solution)
        if jk.get("n_subsets"):
            print(f"jackknife (n={jk['n_subsets']}): Mw +/-{jk['mw_std']}, "
                  f"DC +/-{jk['dc_std']}%, max mechanism rotation "
                  f"{jk['max_tensor_rotation_deg']} deg")
    finally:
        os.chdir(cwd)

    invert.save_solution(solution, band_dir)
    try:
        wf = figure.plot_band_waveforms(
            band_dir, solution,
            band_dir.parent / f"{event.public_id}_station_waveforms_"
                              f"{config.band_tag(band)}.jpg")
        print(f"  all-station waveforms: {wf.name}")
    except Exception as e:  # noqa: BLE001 — diagnostic figure, not data
        print(f"  WARNING: all-station waveform figure failed: {e}")
    if config.is_solved(solution):
        p, q = solution["preferred"], solution["quality"]
        print(f"band result: depth {p['depth_km']:g} km "
              f"[{q['depth_range_km'][0]:g}-{q['depth_range_km'][1]:g}], "
              f"Mw {p['mw']:.2f}, VR {p['vr']:.1f}%, DC {p['pdc']:.0f}%, "
              f"{q['n_stations_used']} stations, grade {q['grade']}")
    return solution


def _purge_previous_run(event_dir: Path, pid: str) -> None:
    """A reprocessed event starts clean: the previous run's figures, draft
    email and per-band working directories are removed first, so a run
    that ends with no solution can never sit beside an older run's
    displacement field and waveform fits. The previous solution.json is
    left in place until the new one overwrites it, so a crash mid-run
    does not drop the event from the catalogue."""
    for f in event_dir.glob(f"{pid}_*.jpg"):
        f.unlink()
    draft = event_dir / "draft_email.txt"
    if draft.exists():
        draft.unlink()
    for band_dir_ in event_dir.glob("band_*"):
        shutil.rmtree(band_dir_, ignore_errors=True)


def _cleanup(event_dir: Path) -> None:
    """Staged Green's functions and SAC data are regenerable; mtinv.in,
    per-band solution.json and figures stay as provenance."""
    for band_dir_ in event_dir.glob("band_*"):
        shutil.rmtree(band_dir_ / "greens", ignore_errors=True)
        for dat in band_dir_.glob("*.dat"):
            dat.unlink()


def _rebuild_tables() -> None:
    catalogue.build_catalogue()
    import station_performance
    station_performance.build_station_performance()


def _archive_no_solution(event, event_dir: Path, solutions: dict,
                         best_tag: str, debug: bool) -> dict:
    rec = solutions[best_tag]
    rec["chosen_band"] = best_tag
    rec["band_search"] = {
        tag: {"status": config.STATUS_NO_SOLUTION,
              "stage": s["abort"]["stage"], "best_vr": s["abort"]["best_vr"]}
        for tag, s in solutions.items()}
    # events with no solution share one directory: NOSOL/<pid>.json plus
    # the pid-prefixed station map and all-station waveform figures
    pid = event.public_id
    out = config.no_solution_path(pid)
    out.parent.mkdir(parents=True, exist_ok=True)
    for old in out.parent.glob(f"{pid}_*.jpg"):
        old.unlink()
    try:
        figure.plot_station_ledger_map(
            rec, out.parent / f"{pid}_station_map.jpg")
    except Exception as e:  # noqa: BLE001
        print(f"  WARNING: station ledger map failed: {e}")
    for wf in event_dir.glob(f"{pid}_station_waveforms_*.jpg"):
        shutil.move(str(wf), out.parent / wf.name)
    out.write_text(json.dumps(rec, indent=2))
    if debug:
        _cleanup(event_dir)      # keep the working directory for inspection
    else:
        shutil.rmtree(event_dir, ignore_errors=True)
    _rebuild_tables()
    print(f"archived as {out.relative_to(config.EVENTS_DIR)}")
    print(f"\nNO COHERENT SOLUTION ({best_tag}): {rec['abort']['stage']} — "
          f"{rec['abort']['reason']}\nsolution: {out}")
    return rec


def process_event(public_id: str, debug: bool = False,
                  band: tuple[float, float] | None = None) -> dict:
    # --- task 1 -------------------------------------------------------------
    event = get_event(public_id)
    print(f"{event.public_id}: M{event.prelim_mag:.1f} {event.locality}, "
          f"depth {event.depth_km:g} km, quality={event.quality}")
    assert event.quality != "deleted", f"{public_id} is marked deleted by GeoNet"
    event_dir = (config.find_event_dir(event.public_id)
                 or config.EVENTS_DIR / event.public_id)
    event_dir.mkdir(parents=True, exist_ok=True)
    _purge_previous_run(event_dir, event.public_id)
    model = config.model_for_event(event.latitude, event.longitude)
    print(f"velocity model: {model}")

    # --- tasks 2-3, per band in order of preference -------------------------
    bands = [band] if band else config.band_candidates(event.prelim_mag)
    solutions = {}
    for b in bands:
        tag = config.band_tag(b)
        try:
            solutions[tag] = process_band(event, b, event_dir / tag, debug,
                                          model)
        except AssertionError as e:
            print(f"band {tag} failed: {e}")
            continue
        s = solutions[tag]
        if config.is_solved(s) and s["quality"]["passed"]:
            break   # first band that passes its gate wins
        print(f"band {tag}: {'grade ' + s['quality']['grade'] if config.is_solved(s) else 'no coherent solution'}"
              + (", trying the next" if b is not bands[-1] else ""))
    assert solutions, "every filter band failed"

    tags = [config.band_tag(b) for b in bands if config.band_tag(b) in solutions]
    solved = [t for t in tags if config.is_solved(solutions[t])]
    if not solved:
        best_tag = max(tags, key=lambda t: solutions[t]["abort"]["best_vr"])
        return _archive_no_solution(event, event_dir, solutions, best_tag, debug)
    passing = [t for t in solved if solutions[t]["quality"]["passed"]]
    # a passing band wins outright; otherwise the best VR among the solved
    best_tag = passing[0] if passing else max(
        solved, key=lambda t: solutions[t]["preferred"]["vr"])
    best = solutions[best_tag]
    best["chosen_band"] = best_tag
    best["band_search"] = {
        tag: ({"vr": s["preferred"]["vr"], "mw": s["preferred"]["mw"],
               "depth_km": s["preferred"]["depth_km"], "pdc": s["preferred"]["pdc"],
               "n_stations": s["quality"]["n_stations_used"],
               "grade": s["quality"]["grade"]}
              if config.is_solved(s) else
              {"status": config.STATUS_NO_SOLUTION,
               "stage": s["abort"]["stage"], "best_vr": s["abort"]["best_vr"]})
        for tag, s in solutions.items()}

    # --- task 4 -------------------------------------------------------------
    forward = okada_forward.forward_both_planes(best)
    best["forward_model"] = {
        "peak_abs_cm": forward["peak_abs_m"] * 100.0,
        "detectable": forward["detectable"],
        "plane1_fault": forward["plane1"]["fault"],
        "plane2_fault": forward["plane2"]["fault"],
    }
    try:
        passes = nisar_dates.nisar_passes(event.longitude, event.latitude)
    except Exception as e:  # noqa: BLE001 - CMR outage must not kill the MT
        print(f"NISAR query failed: {e}")
        passes = []
    best["nisar_passes"] = passes

    # --- task 5 -------------------------------------------------------------
    pid = event.public_id
    fig_path = figure.make_share_figure(
        best, forward, passes, event_dir / f"{pid}_stations_displacement_field.jpg")
    depth_fig = figure.plot_depth_sensitivity(
        best, event_dir / f"{pid}_depth_sensitivity.jpg")
    for i, bb in enumerate(sorted((event_dir / best_tag).glob("bbwaves.*.jpg"))):
        shutil.copy(bb, event_dir / f"{pid}_waveform_fits_{i:02d}.jpg")
    print(f"figures: {fig_path}, {depth_fig}")

    history = json.loads(config.STATE_FILE.read_text())["published"] \
        if config.STATE_FILE.exists() else []
    best["publish_decision"] = trigger.publish_decision(best, forward, history)
    subject, body = publish.draft_text(best, forward, passes)
    (event_dir / "draft_email.txt").write_text(f"{subject}\n\n{body}\n")
    out = invert.save_solution(best, event_dir)
    _rebuild_tables()
    # the overview map (events/solutions_map.jpg) is a build product like
    # the tables, regenerated by `pixi run tables` (sweep checkpoints, the
    # CI commit step, sync_repo) rather than after every event
    if not debug:
        _cleanup(event_dir)
    canonical = config.EVENTS_DIR / config.event_dir_name(
        pid, best["preferred"]["mw"], best["preferred"]["depth_km"],
        event.locality, event.origin_time)
    if event_dir != canonical:
        if canonical.exists():
            shutil.rmtree(canonical)
        event_dir.rename(canonical)
        print(f"archived as {canonical.name}")
    # a solved event supersedes any earlier no-solution record
    nosol = config.no_solution_path(pid)
    for old in [nosol, *nosol.parent.glob(f"{pid}_*.jpg")]:
        if old.exists():
            old.unlink()

    p, q = best["preferred"], best["quality"]
    print(f"predicted peak displacement: {forward['peak_abs_m']*100:.2f} cm "
          f"(detectable: {forward['detectable']})\n"
          f"publish decision: {best['publish_decision']['publish']} — "
          f"{best['publish_decision']['reasons']}")
    print(f"\nCHOSEN ({best_tag}): depth {p['depth_km']:g} km "
          f"[{q['depth_range_km'][0]:g}-{q['depth_range_km'][1]:g}], "
          f"Mw {p['mw']:.2f}, VR {p['vr']:.1f}%, DC {p['pdc']:.0f}%, "
          f"grade {q['grade']} ({q['n_stations_used']} stations)\n"
          f"plane1 {p['plane1']}\nplane2 {p['plane2']}\n"
          f"quality: {q}\nsolution: {out}")
    return best


def _parse_band(text: str) -> tuple[float, float]:
    lo, hi = (float(x) for x in text.split("-"))
    assert 0 < lo < hi, f"bad band {text}"
    return (lo, hi)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--event", required=True, help="GeoNet publicID")
    ap.add_argument("--debug", action="store_true",
                    help="save troubleshooting figures of every stage")
    ap.add_argument("--band", type=_parse_band, default=None,
                    help="force one passband in Hz, e.g. 0.02-0.10")
    args = ap.parse_args()
    process_event(args.event, debug=args.debug, band=args.band)
