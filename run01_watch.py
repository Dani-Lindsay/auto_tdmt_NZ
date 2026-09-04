"""Poll GeoNet for new events above the processing floor and process them.

    pixi run python run01_watch.py            # poll + list what would run
    pixi run python run01_watch.py --process  # poll + process + (maybe) email

State lives in events/index.json:
    {"processed": {publicID: {...summary...}}, "published": [...]}
In CI this file is committed back to the repo after each run, which both
dedupes across cron runs and archives every solution publicly.

One polite API call per invocation.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone

import config
import trigger
from geonet import recent_quakes


def load_state() -> dict:
    if config.STATE_FILE.exists():
        return json.loads(config.STATE_FILE.read_text())
    return {"processed": {}, "published": []}


def save_state(state: dict) -> None:
    config.STATE_FILE.parent.mkdir(parents=True, exist_ok=True)
    config.STATE_FILE.write_text(json.dumps(state, indent=2))


MAX_ATTEMPTS = 3


def new_events(state: dict) -> list:
    events = []
    for ev in recent_quakes(mmi=3):
        prior = state["processed"].get(ev.public_id)
        if prior is not None:
            # retry transient failures (NRT hiccups) a bounded number of times
            if not (prior.get("status") == "failed"
                    and prior.get("attempts", 1) < MAX_ATTEMPTS):
                continue
        ok, reason = trigger.passes_processing_floor(ev)
        if ok:
            events.append(ev)
        # below-floor events are not recorded: a magnitude revision upward
        # on a later poll should still trigger processing
    return events


def main(process: bool) -> None:
    state = load_state()
    events = new_events(state)
    if not events:
        print("no new events above the processing floor")
        return

    for ev in events:
        print(f"NEW: {ev.public_id} M{ev.prelim_mag:.1f} {ev.locality} "
              f"depth {ev.depth_km:g} km")
        if not process:
            continue

        from run02_process import process_event
        try:
            solution = process_event(ev.public_id)
        except Exception as e:  # noqa: BLE001 - one bad event must not stop the rest
            print(f"{ev.public_id} FAILED: {e}")
            prior = state["processed"].get(ev.public_id, {})
            state["processed"][ev.public_id] = {
                "status": "failed", "error": str(e)[:500],
                "attempts": prior.get("attempts", 0) + 1,
                "prelim_mag": ev.prelim_mag,
                "run_utc": datetime.now(timezone.utc).isoformat(),
            }
            save_state(state)
            continue

        if not config.is_solved(solution):
            ab = solution["abort"]
            state["processed"][ev.public_id] = {
                "status": "no_solution",
                "stage": ab["stage"], "reason": ab["reason"],
                "best_vr": ab["best_vr"],
                "run_utc": datetime.now(timezone.utc).isoformat(),
            }
            save_state(state)
            print(f"{ev.public_id}: no coherent solution "
                  f"({ab['stage']}) — archived, not published")
            continue

        decision = solution["publish_decision"]
        summary = {
            "status": "ok",
            "mw": solution["preferred"]["mw"],
            "depth_km": solution["preferred"]["depth_km"],
            "vr": solution["preferred"]["vr"],
            "pdc": solution["preferred"]["pdc"],
            "publish": decision["publish"],
            "run_utc": datetime.now(timezone.utc).isoformat(),
        }

        if decision["publish"]:
            from run03_publish import publish_event
            try:
                publish_event(ev.public_id, state=state)
                summary["published"] = True
            except Exception as e:  # noqa: BLE001
                print(f"{ev.public_id} publish FAILED: {e}")
                summary["published"] = False
                summary["publish_error"] = str(e)[:500]

        state["processed"][ev.public_id] = summary
        save_state(state)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--process", action="store_true",
                    help="process new events (default: just list them)")
    args = ap.parse_args()
    main(process=args.process)
