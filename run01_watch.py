"""Poll GeoNet for new events above the processing floor and process them.

    pixi run python run01_watch.py            # poll + list what would run
    pixi run python run01_watch.py --process  # poll + process + (maybe) email
    pixi run python run01_watch.py --process --budget-minutes 300

Runs DAILY in CI (see .github/workflows/watch.yml). A 10-minute cron was
throttled by GitHub to a median 3.5-hour gap — an unreliable schedule is
worse than an honest one, and NISAR/Sentinel passes are days apart, so a
daily sweep loses nothing that matters. Urgent events can always be run
by hand (workflow_dispatch, or run02_process.py locally).

Because the poll returns everything recent and the state file records
what has been done, a run that cannot finish its list is not a problem:
it processes the LARGEST events first and stops starting new ones when
its time budget runs out, leaving the rest for the next run. That is
what keeps a 25-event aftershock day inside one CI job.

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


def main(process: bool, budget_minutes: float | None = None) -> None:
    state = load_state()
    events = new_events(state)
    if not events:
        print("no new events above the processing floor")
        return

    # biggest first: if the budget runs out, the events that matter most
    # for an InSAR response are the ones already done
    events.sort(key=lambda e: -e.prelim_mag)
    print(f"{len(events)} event(s) to process, largest first"
          + (f"; budget {budget_minutes:g} min" if budget_minutes else ""))
    started = datetime.now(timezone.utc)

    for ev in events:
        if budget_minutes is not None:
            spent = (datetime.now(timezone.utc) - started).total_seconds() / 60
            if spent >= budget_minutes:
                left = [e.public_id for e in events[events.index(ev):]]
                print(f"time budget spent ({spent:.0f} min): leaving "
                      f"{len(left)} event(s) for the next run: "
                      f"{', '.join(left[:8])}"
                      + (" ..." if len(left) > 8 else ""))
                break
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
    ap.add_argument("--budget-minutes", type=float, default=None,
                    help="stop starting new events after this long; the "
                         "rest are picked up by the next run")
    ap.add_argument("--process", action="store_true",
                    help="process new events (default: just list them)")
    args = ap.parse_args()
    main(process=args.process, budget_minutes=args.budget_minutes)
