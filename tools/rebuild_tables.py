"""Rebuild every build product derived from the archive: the three tables
(catalogue.csv, not_published.csv, station_ledger.jsonl), the per-station
aggregate (station_performance.csv) and the overview map.

    pixi run tables

Used after a `git pull --rebase` by the CI watcher and the local sweep,
so the products are regenerated from the merged archive rather than
merged as files (they are never merged by hand).
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import catalogue  # noqa: E402
import config  # noqa: E402
import figure  # noqa: E402
import station_performance  # noqa: E402

def latest_event_figure() -> None:
    """Copy the stations + displacement figure of the most recent solved
    event (by origin time) to events/latest_stations_displacement_field.jpg,
    the fixed path the README embeds, so the home page always shows the
    latest automated solution."""
    import json
    import shutil

    newest, newest_t = None, ""
    for path in config.solution_paths():
        if path.parent.name == config.NOSOL_DIR_NAME:
            continue
        sol = json.loads(path.read_text())
        if not config.is_solved(sol):
            continue
        t = sol["event"]["origin_time"]
        if t > newest_t:
            newest, newest_t = path.parent, t
    if newest is None:
        return
    pid = newest.name.split("_")[0]
    src = newest / f"{pid}_stations_displacement_field.jpg"
    if src.exists():
        shutil.copy(src, config.EVENTS_DIR / "latest_stations_displacement_field.jpg")
        print(f"latest event figure: {newest.name}")


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    catalogue.build_catalogue(config.EVENTS_DIR)
    station_performance.build_station_performance(config.EVENTS_DIR)
    try:
        figure.make_overview_map(config.EVENTS_DIR,
                                 config.EVENTS_DIR / "solutions_map.jpg")
    except Exception as e:  # noqa: BLE001 - map is cosmetic, tables are not
        print(f"overview map failed: {e}")
    latest_event_figure()
