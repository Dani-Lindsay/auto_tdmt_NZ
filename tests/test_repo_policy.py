"""Repository-size policy: per-event figures are committed only for events
from FIGURES_FROM onward (older ones are in the figures-pre-2026-07
release). The rule lives in .gitignore; these tests hold it to every
event's actual origin time, so a wrong pattern fails here rather than
re-committing a gigabyte of figures."""

import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FIGURES_FROM = datetime(2026, 7, 1, tzinfo=timezone.utc)


def _origin(sol: dict) -> datetime:
    return datetime.fromisoformat(sol["event"]["origin_time"].replace("Z", "+00:00"))


def _ignored(paths: list[str]) -> set[str]:
    """Paths .gitignore would ignore (tracked or not)."""
    r = subprocess.run(["git", "check-ignore", "--no-index", "--stdin"],
                       input="\n".join(paths), capture_output=True,
                       text=True, cwd=REPO)
    assert r.returncode in (0, 1), r.stderr
    return set(r.stdout.split())


SECONDS_PER_ID = 31.968   # fitted over every archived event (2021-2026)
CUTOFF_ID = 489190        # first ID on/after 2026-07-01 00:00 UTC


def test_publicid_encodes_origin_time():
    # the patterns rely on publicID = <year>p<~seconds since 1 Jan UTC / 31.968>;
    # the scatter (a few steps) is GeoNet's detection-to-origin lag
    sol_paths = sorted((REPO / "events").glob("*/solution.json")) + \
        sorted((REPO / "events" / "NOSOL").glob("*.json"))
    assert sol_paths, "no events in the repo archive"
    for p in sol_paths:
        sol = json.loads(p.read_text())
        pid, t = sol["event"]["public_id"], _origin(sol)
        year_start = datetime(int(pid[:4]), 1, 1, tzinfo=timezone.utc)
        predicted = (t - year_start).total_seconds() / SECONDS_PER_ID
        assert abs(int(pid[5:]) - predicted) < 20, pid
    cutoff = (FIGURES_FROM - datetime(2026, 1, 1, tzinfo=timezone.utc)).total_seconds()
    assert math.ceil(cutoff / SECONDS_PER_ID) == CUTOFF_ID   # first ID on/after midnight


def test_gitignore_splits_figures_at_the_cutoff():
    expect = {}
    for p in sorted((REPO / "events").glob("*/solution.json")):
        sol = json.loads(p.read_text())
        fig = f"events/{p.parent.name}/{sol['event']['public_id']}_depth_sensitivity.jpg"
        expect[fig] = _origin(sol) < FIGURES_FROM
    for p in sorted((REPO / "events" / "NOSOL").glob("*.json")):
        sol = json.loads(p.read_text())
        fig = f"events/NOSOL/{sol['event']['public_id']}_station_map.jpg"
        expect[fig] = _origin(sol) < FIGURES_FROM
    # the boundary itself, independent of which events happen to exist
    expect[f"events/X/2026p{CUTOFF_ID - 1}_a.jpg"] = True
    expect[f"events/X/2026p{CUTOFF_ID}_a.jpg"] = False
    expect["events/X/2026p479999_a.jpg"] = True
    expect["events/X/2026p490000_a.jpg"] = False
    expect["events/X/2025p999999_a.jpg"] = True
    expect["events/X/2027p000001_a.jpg"] = False
    expect["events/solutions_map.jpg"] = False
    expect["events/latest_stations_displacement_field.jpg"] = False
    ignored = _ignored(list(expect))
    wrong = [f for f, pre in expect.items() if (f in ignored) != pre]
    assert not wrong, f"gitignore disagrees with the cutoff for: {wrong[:5]}"


def test_no_tracked_figure_before_the_cutoff():
    tracked = subprocess.run(["git", "ls-files", "events"], capture_output=True,
                             text=True, cwd=REPO, check=True).stdout.split()
    old = _ignored([f for f in tracked if f.endswith(".jpg")])
    assert not old, f"{len(old)} pre-cutoff figures are tracked, e.g. {sorted(old)[:3]}"
