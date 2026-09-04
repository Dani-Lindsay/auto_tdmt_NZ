"""Re-grade archived solutions under the CURRENT rubric, without
reprocessing anything.

Grades are derived entirely from numbers already stored in each
solution.json (VR, %DC, per-station VRs, jackknife rotation, depth
flags, station azimuths), so a rubric change can be applied to the whole
archive in seconds and reviewed before any resweep.

    pixi run python tools/regrade_archive.py            # dry run, report
    pixi run python tools/regrade_archive.py --apply    # rewrite grades

Old solutions predate some fields; those are recomputed here where
possible (edge_artifact from the depth search, the azimuth pair from
the used stations) and reported as unknown otherwise.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import config  # noqa: E402
import invert  # noqa: E402


def backfill_flags(sol: dict) -> dict:
    """Recompute the depth flags a v3 solution never stored."""
    rows = sorted(sol.get("depth_search", []), key=lambda r: r["depth_km"])
    flags = dict(sol.get("depth_pick_flags", {}))
    if rows and "edge_artifact" not in flags:
        vrs = [r["vr"] for r in rows]
        i_max = max(range(len(vrs)), key=lambda i: vrs[i])
        span = invert.plateau_indices(vrs, i_max)
        plateau_km = rows[span[-1]]["depth_km"] - rows[span[0]]["depth_km"]
        edges = (min(config.GF_DEPTHS_KM), max(config.GF_DEPTHS_KM))
        flags["edge_artifact"] = bool(
            sol["preferred"]["depth_km"] in edges and plateau_km == 0.0)
    return flags


def main(apply: bool, events_dir: Path) -> None:
    paths = sorted(events_dir.glob("*/solution.json"))
    if not paths:
        print(f"no solutions under {events_dir}")
        return
    moves: Counter = Counter()
    changed: list[tuple[str, str, str, str]] = []
    unsolved = 0
    for p in paths:
        sol = json.loads(p.read_text())
        if not config.is_solved(sol):
            unsolved += 1
            continue
        before = sol.get("quality", {}).get("grade", "?")
        sol["depth_pick_flags"] = backfill_flags(sol)
        quality = invert.quality_gates(sol)
        after = quality["grade"]
        moves[(before, after)] += 1
        if before != after:
            why = []
            for name, ok in quality["checks"].items():
                if not ok:
                    why.append(name)
            changed.append((sol["event"]["public_id"], before, after,
                            ",".join(why) or "-"))
        if apply:
            sol["quality"] = quality
            p.write_text(json.dumps(sol, indent=2))

    total = sum(moves.values())
    print(f"{total} solved solutions re-graded"
          + (f" ({unsolved} with no coherent solution skipped)"
             if unsolved else ""))
    print("\nbefore -> after")
    for (b, a), n in sorted(moves.items()):
        mark = "" if b == a else "   <-- moved"
        print(f"  {b} -> {a}: {n}{mark}")
    old = Counter(b for (b, _a), n in moves.items() for _ in range(n))
    new = Counter(a for (_b, a), n in moves.items() for _ in range(n))
    print(f"\nold distribution: {dict(sorted(old.items()))}")
    print(f"new distribution: {dict(sorted(new.items()))}")
    if changed:
        print(f"\n{len(changed)} events changed grade "
              f"(first 40, with the checks they now fail):")
        for pid, b, a, why in changed[:40]:
            print(f"  {pid} {b}->{a}  {why}")
    print("\n(dry run — nothing written; pass --apply to rewrite grades)"
          if not apply else "\ngrades rewritten in place")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true",
                    help="rewrite quality blocks in solution.json")
    ap.add_argument("--events-dir", default=None)
    args = ap.parse_args()
    main(args.apply,
         Path(args.events_dir) if args.events_dir else config.EVENTS_DIR)
