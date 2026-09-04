"""Machine-readable catalogue of all archived solutions.

Regenerates events/catalogue.csv from every events/<publicID>/solution.json
after each processed event — the CSV is always derived from the archived
sidecars (single source of truth), never edited by hand.

Event detection and original hypocentres (time, location, preliminary
magnitude, initial depth) are sourced from GeoNet; the GeoNet_M and
GeoNet_depth columns carry that original solution alongside our inverted
Mw and centroid depth (CD) so the revision is explicit per event.

Column conventions follow the published NZ regional CMT solutions CSV
(GeoNet/data moment-tensor/GeoNet_CMT_solutions.csv) where they overlap so
the two are directly comparable; note their (and our) MT elements are in
units of 1e20 dyne-cm.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import config

_FLAG_NAMES = {
    "vr_floor": "low_VR",
    "no_passengers": "station_not_fitting",
    "az_pair_90": "no_90deg_azimuth_pair",
    "jackknife_stable": "unstable_mechanism",
    "depth_interior": "grid_edge_depth",
    "dc_floor": "low_DC",
    # v3 names, kept so old archived solutions still read sensibly
    "min_stations": "few_stations",
    "az_gap_ok": "wide_az_gap",
}


def _quality_flag(quality: dict) -> str:
    if quality.get("passed"):
        return "True"
    failed = [
        _FLAG_NAMES.get(name, name)
        for name, ok in quality.get("checks", {}).items() if not ok
    ]
    grade = quality.get("grade", "?")
    return ";".join(failed) if failed else f"grade_{grade}"


def _publish_flag(decision: dict) -> str:
    if decision.get("publish"):
        return "True"
    tags = decision.get("reason_tags")
    return ";".join(tags) if tags else "unpublished"


COLUMNS = [
    "PublicID", "Date", "Latitude", "Longitude",
    "strike1", "dip1", "rake1", "strike2", "dip2", "rake2",
    "GeoNet_M", "GeoNet_depth", "GeoNet_depth_unc", "Mw", "Depth",
    "Depth_VRmax", "Depth_DCmax", "Plateau_km", "Mo",
    "NS", "AzGap", "Grade", "DC", "CLVD", "VR",
    "Jk_n", "Jk_Mw_std", "Jk_DC_std", "Jk_rot_deg",
    "PredDisp_cm", "Detectable",
    "Mxx", "Mxy", "Mxz", "Myy", "Myz", "Mzz",
    "Band", "Model", "Status", "Selection", "Code",
    "quality_flag", "publish_flag", "published",
]


def build_catalogue(events_dir: Path | None = None) -> Path | None:
    """Scan <events_dir>/*/solution.json -> <events_dir>/catalogue.csv.

    Returns the CSV path, or None if there are no solutions yet.
    """
    events_dir = events_dir or config.EVENTS_DIR
    solutions = sorted(events_dir.glob("*/solution.json"))
    if not solutions:
        return None

    published_ids: set[str] = set()
    if config.STATE_FILE.exists():
        state = json.loads(config.STATE_FILE.read_text())
        published_ids = {p["public_id"] for p in state.get("published", [])}

    rows = []
    for path in solutions:
        s = json.loads(path.read_text())
        ev = s["event"]
        if not config.is_solved(s):
            # no coherent solution: the row carries the origin and the
            # reason, and NO mechanism numbers at all
            rows.append({
                "PublicID": ev["public_id"],
                "Date": ev["origin_time"],
                "Latitude": round(ev["latitude"], 4),
                "Longitude": round(ev["longitude"], 4),
                "GeoNet_M": round(ev["prelim_mag"], 2),
                "GeoNet_depth": round(ev["depth_km"], 1),
                "NS": 0,
                "Grade": "X",
                "Band": s.get("chosen_band", "").replace("band_", ""),
                "Model": s.get("provenance", {}).get("velocity_model", ""),
                "Status": config.STATUS_NO_SOLUTION,
                "Selection": s.get("provenance", {}).get(
                    "selection_version", ""),
                "Code": s.get("provenance", {}).get("code_commit", ""),
                "quality_flag": (f"no_coherent_solution:"
                                 f"{s['abort']['stage']}"),
                "publish_flag": "no_coherent_solution",
                "published": False,
            })
            continue
        p = s["preferred"]
        mt = s["preferred"]["tensor_dyne_cm"]  # XYZ basis, dyne-cm
        rows.append({
            "PublicID": ev["public_id"],
            "Date": ev["origin_time"],
            "Latitude": round(ev["latitude"], 4),
            "Longitude": round(ev["longitude"], 4),
            "strike1": round(p["plane1"]["strike"], 1),
            "dip1": round(p["plane1"]["dip"], 1),
            "rake1": round(p["plane1"]["rake"], 1),
            "strike2": round(p["plane2"]["strike"], 1),
            "dip2": round(p["plane2"]["dip"], 1),
            "rake2": round(p["plane2"]["rake"], 1),
            "GeoNet_M": round(ev["prelim_mag"], 2),
            "GeoNet_depth": round(ev["depth_km"], 1),
            "GeoNet_depth_unc": ev.get("depth_unc_km", ""),
            "Mw": round(p["mw"], 2),
            # our centroid depth beside our Mw: a large Mw revision is
            # often explained by the depth revision visible in the adjacent
            # GeoNet_depth column
            "Depth": p["depth_km"],
            "Depth_VRmax": s.get("depth_pick_flags", {}).get(
                "vr_max_depth_km", ""),
            "Depth_DCmax": s.get("depth_pick_flags", {}).get(
                "dc_max_depth_km", ""),
            "Plateau_km": s.get("depth_pick_flags", {}).get(
                "plateau_km", ""),
            # Mo in dyne-cm; MT elements in 1e20 dyne-cm (GeoNet convention)
            "Mo": f"{p['m0_dyne_cm']:.3e}",
            "NS": s["quality"]["n_stations_used"],
            "AzGap": s["quality"].get("azimuthal_gap_deg", ""),
            "Grade": s["quality"].get("grade", ""),
            "DC": round(p["pdc"], 1),
            "CLVD": round(p["pclvd"], 1),
            "VR": round(p["vr"], 1),
            "Jk_n": s.get("jackknife", {}).get("n_subsets", ""),
            "Jk_Mw_std": s.get("jackknife", {}).get("mw_std", ""),
            "Jk_DC_std": s.get("jackknife", {}).get("dc_std", ""),
            "Jk_rot_deg": s.get("jackknife", {}).get(
                "max_tensor_rotation_deg", ""),
            "PredDisp_cm": round(
                s.get("forward_model", {}).get("peak_abs_cm", 0), 3)
                if s.get("forward_model") else "",
            "Detectable": s.get("forward_model", {}).get("detectable", ""),
            **{
                k.capitalize().replace("m", "M", 1): round(v / 1e20, 4)
                for k, v in mt.items()
            },
            "Band": s.get("chosen_band", "").replace("band_", ""),
            "Model": s["provenance"]["velocity_model"],
            # which code produced this row (the catalogue is
            # self-describing when it mixes vintages)
            "Status": s.get("status", config.STATUS_SOLVED),
            "Selection": s["provenance"].get("selection_version", "v3"),
            "Code": s["provenance"].get("code_commit", ""),
            # "True" when the solution passed, else the actual exit
            # reason(s) so a reader can see WHY it fell short at a glance
            "quality_flag": _quality_flag(s["quality"]),
            "publish_flag": _publish_flag(s.get("publish_decision", {})),
            "published": ev["public_id"] in published_ids,
        })

    rows.sort(key=lambda r: r["Date"])
    out = events_dir / "catalogue.csv"
    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in COLUMNS})
    print(f"catalogue: {len(rows)} solutions -> {out}")
    return out


if __name__ == "__main__":
    build_catalogue()
