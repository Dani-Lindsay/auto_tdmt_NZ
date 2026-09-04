"""Replay the reviewer's station labels against the current selection code.

The labels in docs/review_labels/station_labels.json are D. Lindsay's
per-station keep / maybe / toss calls, read off the all-station waveform
figures (verbatim notes alongside them in notes_*.md). They are the
ACCEPTANCE TEST for station selection: any change to the funnel must be
judged by whether it reproduces these calls, not by whether VR went up.

    pixi run python run07_replay_labels.py                 # all events
    pixi run python run07_replay_labels.py --event 2026p348732
    pixi run python run07_replay_labels.py --events-dir /tmp/scratch

Scoring, per station:
  keep  -> must be USED          (a miss is a false rejection)
  toss  -> must NOT be used      (a miss is a false acceptance)
  maybe -> free, reported only
Per event: whether an event the reviewer called unsolvable actually
aborts, and whether the preferred depth lands in her plausible range.

Needs network access (waveform download) unless the event is already in
the waveform cache. Not a pytest — it inverts real data.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

LABELS = (Path(__file__).resolve().parent / "docs" / "review_labels"
          / "station_labels.json")


def _station_name(sid: str) -> str:
    parts = sid.split(".")
    return parts[1] if len(parts) >= 2 else sid


def replay(labels: dict, only: str | None = None) -> dict:
    import config
    from run02_process import process_event

    results = []
    for spec in labels["events"]:
        pid = spec["public_id"]
        if only and pid != only:
            continue
        want_status = spec.get("event", {}).get("expect", "either")
        print(f"\n{'=' * 70}\n{pid}  (reviewer expects: {want_status})")
        t0 = time.time()
        try:
            sol = process_event(pid)
        except Exception as e:  # noqa: BLE001 - report, keep going
            print(f"  FAILED: {str(e)[:200]}")
            results.append({"public_id": pid, "error": str(e)[:200]})
            continue
        runtime = time.time() - t0

        solved = config.is_solved(sol)
        used = {_station_name(f"{r['network']}.{r['station']}")
                for r in sol.get("stations_used", [])}
        fates = {}
        for d in sol.get("stations_dropped", []):
            fates[_station_name(d["station"])] = d.get("reason", "")[:70]

        hits = misses = free = 0
        rows = []
        for sta, label in spec.get("stations", {}).items():
            in_use = sta in used
            if label == "keep":
                ok = in_use
            elif label == "toss":
                ok = not in_use
            else:
                ok = None
            rows.append((sta, label, "USED" if in_use else "not used",
                         "" if in_use else fates.get(sta, "not in pool"), ok))
            if ok is True:
                hits += 1
            elif ok is False:
                misses += 1
            else:
                free += 1

        status_ok = (want_status == "either"
                     or (want_status == "solved") == solved)
        depth_ok = None
        rng = spec.get("event", {}).get("depth_range_km")
        if rng and solved:
            depth_ok = rng[0] <= sol["preferred"]["depth_km"] <= rng[1]

        if solved:
            p, q = sol["preferred"], sol["quality"]
            head = (f"  -> Mw {p['mw']:.2f}  depth {p['depth_km']:g} km  "
                    f"VR {p['vr']:.0f}%  DC {p['pdc']:.0f}%  "
                    f"grade {q['grade']}  ({len(used)} stations)")
        else:
            head = (f"  -> NO COHERENT SOLUTION "
                    f"({sol['abort']['stage']}: {sol['abort']['reason'][:60]})")
        print(head + f"  [{runtime:.0f}s]")
        print(f"     status expectation: "
              f"{'MET' if status_ok else 'NOT MET'}"
              + ("" if depth_ok is None else
                 f" | depth {rng[0]:g}-{rng[1]:g} km: "
                 f"{'MET' if depth_ok else 'NOT MET'}"))
        if rows:
            print(f"     {'station':8s} {'label':6s} {'outcome':9s} why")
            for sta, label, outcome, why, ok in sorted(
                    rows, key=lambda r: (r[1], r[0])):
                mark = "  " if ok is None else (" +" if ok else " X")
                print(f"    {mark} {sta:8s} {label:6s} {outcome:9s} {why}")
            print(f"     agreement: {hits}/{hits + misses} "
                  f"({free} free)" if hits + misses else "")

        results.append({
            "public_id": pid, "solved": solved, "runtime_s": round(runtime),
            "status_ok": status_ok, "depth_ok": depth_ok,
            "hits": hits, "misses": misses, "free": free,
            "stations": [{"station": s, "label": l, "outcome": o,
                          "why": w, "ok": k} for s, l, o, w, k in rows],
            "grade": sol["quality"]["grade"] if solved else "X",
            "mw": sol["preferred"]["mw"] if solved else None,
            "depth_km": sol["preferred"]["depth_km"] if solved else None,
            "vr": sol["preferred"]["vr"] if solved else None,
        })
    return {"results": results}


def report(out: dict) -> None:
    res = [r for r in out["results"] if "error" not in r]
    hits = sum(r["hits"] for r in res)
    misses = sum(r["misses"] for r in res)
    free = sum(r["free"] for r in res)
    status = sum(1 for r in res if r["status_ok"])
    depth = [r for r in res if r["depth_ok"] is not None]
    print(f"\n{'=' * 70}\nSUMMARY over {len(res)} events")
    print(f"  station agreement: {hits}/{hits + misses} "
          f"({100 * hits / (hits + misses):.0f}%) with {free} free 'maybe'"
          if hits + misses else "  no labelled stations")
    print(f"  solved/abort expectation met: {status}/{len(res)}")
    if depth:
        ok = sum(1 for r in depth if r["depth_ok"])
        print(f"  depth in reviewer's range: {ok}/{len(depth)}")
    worst = sorted(res, key=lambda r: -r["misses"])[:5]
    if worst and worst[0]["misses"]:
        print("  events with the most disagreements:")
        for r in worst:
            if r["misses"]:
                print(f"    {r['public_id']}: {r['misses']} misses "
                      f"(grade {r['grade']})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--event", default=None, help="replay one publicID")
    ap.add_argument("--events-dir", default=None,
                    help="scratch AUTO_TDMT_EVENTS (keeps the archive clean)")
    ap.add_argument("--out", default=None, help="write JSON results here")
    args = ap.parse_args()
    if args.events_dir:
        os.environ["AUTO_TDMT_EVENTS"] = args.events_dir
    labels = json.loads(LABELS.read_text())
    out = replay(labels, only=args.event)
    report(out)
    if args.out:
        Path(args.out).write_text(json.dumps(out, indent=2))
        print(f"\nwrote {args.out}")
