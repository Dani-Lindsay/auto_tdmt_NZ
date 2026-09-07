"""Rebuild events_human/catalogue_human.csv from every reviewed
solution.json — a BUILD PRODUCT, never hand-edited, so student Pull
Requests can never conflict on it.

Same columns as the automated catalogue (catalogue.COLUMNS: origin,
planes, Mw, depth range, quality, the full moment tensor, provenance)
plus the review fields, so the two catalogues line up column for column.

Layouts supported:
  events_human/<event_dir>/<reviewer-slug>/solution.json   (current)
  events_human/<event_dir>/solution.json                   (legacy)

    python3 src/catalogue_human.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HUMAN = ROOT / "events_human"
sys.path.insert(0, str(ROOT / "src"))

import catalogue  # noqa: E402

REVIEW_COLUMNS = ["Reviewer", "ReviewDate", "Decision", "Changes", "Notes",
                  "Auto_Mw", "Auto_Depth", "Auto_Grade"]
COLUMNS = [c for c in catalogue.COLUMNS if c != "published"] + REVIEW_COLUMNS


def row_from(sol: dict) -> dict:
    row = catalogue.row_from(sol)
    hr = sol.get("human_review", {})
    ref = sol.get("automated_reference", {})
    if not row.get("Band"):
        b = sol.get("filter_band_hz")
        if b:
            row["Band"] = f"{round(1 / b[1])}-{round(1 / b[0])}s"
    row.update({
        "Reviewer": hr.get("reviewer", ""), "ReviewDate": hr.get("date", ""),
        "Decision": hr.get("decision", ""), "Changes": hr.get("changes", ""),
        "Notes": hr.get("notes", ""),
        "Auto_Mw": ref.get("mw", ""), "Auto_Depth": ref.get("depth_km", ""),
        "Auto_Grade": ref.get("grade", ""),
    })
    return row


def build() -> Path:
    rows = []
    for sol_path in sorted(HUMAN.glob("*/solution.json")) \
            + sorted(HUMAN.glob("*/*/solution.json")):
        try:
            rows.append(row_from(json.loads(sol_path.read_text())))
        except Exception as e:  # noqa: BLE001 - name the bad file, keep going
            print(f"SKIPPED {sol_path}: {e}")
    rows.sort(key=lambda r: (r["PublicID"], r["Reviewer"]))
    out = HUMAN / "catalogue_human.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in COLUMNS})
    print(f"wrote {out} ({len(rows)} reviews)")
    return out


if __name__ == "__main__":
    build()
