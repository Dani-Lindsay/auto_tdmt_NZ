"""Re-grade archived solutions under the CURRENT grade table, without
reprocessing anything.

Grades derive entirely from numbers already stored in each solution.json
(VR, %DC, station count, jackknife rotation), so a change to the table in
auto_tdmt.cfg can be applied to the whole archive in seconds and reviewed
before any resweep.

    pixi run python tools/regrade_archive.py            # dry run, report
    pixi run python tools/regrade_archive.py --apply    # rewrite grades
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import config  # noqa: E402
import invert  # noqa: E402


def main(apply: bool, events_dir: Path) -> None:
    paths = sorted(events_dir.glob("*/solution.json"))
    if not paths:
        print(f"no solutions under {events_dir}")
        return
    moves: Counter = Counter()
    changed = []
    unsolved = 0
    for p in paths:
        sol = json.loads(p.read_text())
        if not config.is_solved(sol):
            unsolved += 1
            continue
        before = sol.get("quality", {}).get("grade", "?")
        quality = invert.quality_gates(sol)
        after = quality["grade"]
        moves[(before, after)] += 1
        if before != after:
            why = [k for k, ok in quality["checks"].items() if not ok]
            changed.append((sol["event"]["public_id"], before, after,
                            ",".join(why) or "-"))
        if apply:
            sol["quality"] = quality
            p.write_text(json.dumps(sol, indent=2))
    total = sum(moves.values())
    print(f"{total} solved solutions re-graded"
          + (f" ({unsolved} with no coherent solution skipped)" if unsolved else ""))
    print("\nbefore -> after")
    for (b, a), n in sorted(moves.items()):
        print(f"  {b} -> {a}: {n}{'' if b == a else '   <-- moved'}")
    old = Counter(); new = Counter()
    for (b, a), n in moves.items():
        old[b] += n; new[a] += n
    print(f"\nold distribution: {dict(sorted(old.items()))}")
    print(f"new distribution: {dict(sorted(new.items()))}")
    if changed:
        print(f"\n{len(changed)} events changed grade (first 40):")
        for pid, b, a, why in changed[:40]:
            print(f"  {pid} {b}->{a}  {why}")
    print("\n(dry run — nothing written; pass --apply to rewrite grades)"
          if not apply else "\ngrades rewritten in place")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--events-dir", default=None)
    args = ap.parse_args()
    main(args.apply, Path(args.events_dir) if args.events_dir else config.EVENTS_DIR)
