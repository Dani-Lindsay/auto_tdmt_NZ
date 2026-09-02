"""Process one GeoNet event end-to-end: waveforms -> GFs -> mttime depth
search -> solution.json + figures. Tries the candidate filter bands for the
event magnitude (BSL-style band menu) and keeps the best-VR solution.

    pixi run python run02_process.py --event 2026p660242 --debug
    pixi run python run02_process.py --event 2026p660242 --band 0.02-0.05
"""

from __future__ import annotations

import argparse
import json
import os
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
    """Run the full chain for one filter band inside band_dir."""
    band_dir.mkdir(parents=True, exist_ok=True)
    green_dir = band_dir / "greens"
    print(f"\n=== band {1/band[1]:.0f}-{1/band[0]:.0f} s "
          f"({band[0]:g}-{band[1]:g} Hz) ===")

    stages: dict | None = {} if debug else None
    used, dropped = waveforms.fetch_and_process(
        event, band_dir, band, stages=stages)
    print(f"stations: {len(used)} used, {len(dropped)} dropped")
    for d in dropped:
        print(f"  dropped {d['station']}: {d['reason']}")

    if debug:
        diag_dir = band_dir / "diagnostics"
        figs = diagnostics.plot_stages(stages, event.to_dict(), band, diag_dir)
        figs.append(diagnostics.plot_station_map(
            event.to_dict(), used, [d["station"] for d in dropped], diag_dir))
        for f in figs:
            print(f"  diagnostic: {f}")

    depths = invert.search_depths(event, model)
    print(f"depth search over {len(depths)} depths: "
          f"{depths[0]:g}-{depths[-1]:g} km")
    greens.stage_event_greens(model, used, depths, band, green_dir)

    cwd = os.getcwd()
    os.chdir(band_dir)
    try:
        inv, kept, rejected = invert.invert_with_rejection(
            event, used, depths, band_dir, green_dir)
        inv.plot(view="waveform", option="preferred", format="jpg", show=False)
    finally:
        os.chdir(cwd)

    solution = invert.summarize(inv, event, kept, dropped + rejected, model)
    solution["filter_band_hz"] = list(band)
    invert.save_solution(solution, band_dir)
    pref = solution["preferred"]
    print(f"band result: depth {pref['depth_km']:g} km, Mw {pref['mw']:.2f}, "
          f"VR {pref['vr']:.1f}%, DC {pref['pdc']:.0f}%")
    return solution


def process_event(public_id: str, debug: bool = False,
                  band: tuple[float, float] | None = None) -> dict:
    event = get_event(public_id)
    print(f"{event.public_id}: M{event.prelim_mag:.1f} {event.locality}, "
          f"depth {event.depth_km:g} km, quality={event.quality}")
    assert event.quality != "deleted", f"{public_id} is marked deleted by GeoNet"

    event_dir = config.EVENTS_DIR / event.public_id
    event_dir.mkdir(parents=True, exist_ok=True)

    model = config.model_for_event(event.latitude, event.longitude)
    print(f"velocity model: {model}")
    bands = [band] if band else config.band_candidates(event.prelim_mag)
    solutions = {}
    for b in bands:
        tag = config.band_tag(b)
        try:
            solutions[tag] = process_band(event, b, event_dir / tag, debug, model)
        except AssertionError as e:
            print(f"band {tag} failed: {e}")
    assert solutions, "every filter band failed"

    # choose among gate-passing bands when any exist; otherwise all bands
    passing = [t for t in solutions if solutions[t]["quality"]["passed"]]
    tags = passing if passing else list(solutions)
    best_tag = tags[invert.pick_preferred(
        [(solutions[t]["preferred"]["vr"], solutions[t]["preferred"]["pdc"])
         for t in tags]
    )]
    best = solutions[best_tag]
    best["chosen_band"] = best_tag
    best["band_search"] = {
        tag: {
            "vr": s["preferred"]["vr"],
            "mw": s["preferred"]["mw"],
            "depth_km": s["preferred"]["depth_km"],
            "pdc": s["preferred"]["pdc"],
            "n_stations": s["quality"]["n_stations_used"],
            "gates_passed": s["quality"]["passed"],
        }
        for tag, s in solutions.items()
    }
    # forward model, NISAR timing, share figure, publish decision
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

    pid = event.public_id
    fig_path = figure.make_share_figure(
        best, forward, passes,
        event_dir / f"{pid}_stations_displacement_field.jpg")
    depth_fig = figure.plot_depth_sensitivity(
        best, event_dir / f"{pid}_depth_sensitivity.jpg")
    # waveform-fit pages (mttime output, untouched) copied up with
    # descriptive event-ID names for the email attachments
    import shutil
    for i, bb in enumerate(
            sorted((event_dir / best_tag).glob("bbwaves.*.jpg"))):
        shutil.copy(bb, event_dir / f"{pid}_waveform_fits_{i:02d}.jpg")
    print(f"figures: {fig_path}, {depth_fig}")

    history = json.loads(config.STATE_FILE.read_text())["published"] \
        if config.STATE_FILE.exists() else []
    best["publish_decision"] = trigger.publish_decision(best, forward, history)

    subject, body = publish.draft_text(best, forward, passes)
    (event_dir / "draft_email.txt").write_text(f"{subject}\n\n{body}\n")

    out = invert.save_solution(best, event_dir)
    catalogue.build_catalogue()

    pref = best["preferred"]
    print(
        f"predicted peak displacement: {forward['peak_abs_m']*100:.2f} cm "
        f"(detectable: {forward['detectable']})\n"
        f"publish decision: {best['publish_decision']['publish']} — "
        f"{best['publish_decision']['reasons']}"
    )
    print(
        f"\nCHOSEN ({best_tag}): depth {pref['depth_km']:g} km, "
        f"Mw {pref['mw']:.2f}, VR {pref['vr']:.1f}%, DC {pref['pdc']:.0f}%\n"
        f"plane1 {pref['plane1']}\nplane2 {pref['plane2']}\n"
        f"band search: {json.dumps(best['band_search'], indent=2)}\n"
        f"quality gates: {best['quality']}\n"
        f"velocity model: {best['provenance']['velocity_model']}\n"
        f"solution: {out}"
    )
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
                    help="force one passband in Hz, e.g. 0.02-0.05")
    args = ap.parse_args()
    process_event(args.event, debug=args.debug, band=args.band)
