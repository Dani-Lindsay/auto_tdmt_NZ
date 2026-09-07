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

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    catalogue.build_catalogue(config.EVENTS_DIR)
    station_performance.build_station_performance(config.EVENTS_DIR)
    try:
        figure.make_overview_map(config.EVENTS_DIR,
                                 config.EVENTS_DIR / "solutions_map.jpg")
    except Exception as e:  # noqa: BLE001 - map is cosmetic, tables are not
        print(f"overview map failed: {e}")
