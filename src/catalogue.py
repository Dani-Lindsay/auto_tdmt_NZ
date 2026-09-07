"""Task 5 — the archive tables, regenerated from every solution.json.

    catalogue.csv          one row per event (solved or not)
    not_published.csv      events with NO solution: attempted but no
                           coherent solution (stage, reason, best VR), and
                           events seen but not attempted (below the floor:
                           too deep, too small, outside the box)
    station_ledger.jsonl   one JSON line per event: a per-station dictionary
                           of numeric values (distance, azimuth, SNR per
                           component, amplitude ratio, own VR, time shift,
                           used / drop class) — the raw material for a
                           station quality index

The tables sit at the REPOSITORY ROOT when the archive is the repo's
events/ (so they are the first thing on the GitHub page); for a scratch
archive elsewhere they sit next to it (tables_dir()).

All three are BUILD PRODUCTS derived from the archived sidecars, never
edited by hand. Column conventions follow the published NZ regional CMT
CSV (GeoNet/data moment-tensor) where they overlap; MT elements are in
1e20 dyne-cm like theirs. The human-review catalogue reuses COLUMNS so
the two catalogues line up column for column.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import config

_FLAG_NAMES = {
    "vr_grade": "VR_below_grade_B",
    "dc_floor": "low_DC",
    "jackknife_stable": "unstable_mechanism",
    "min_stations": "few_stations",
    # older vintages, kept so archived rows still read sensibly
    "vr_floor": "low_VR", "no_passengers": "station_not_fitting",
    "az_pair_90": "no_90deg_azimuth_pair", "depth_interior": "grid_edge_depth",
    "depth_agrees_with_geonet": "depth_far_from_geonet",
    "az_gap_ok": "wide_az_gap",
}


def _quality_flag(quality: dict) -> str:
    if quality.get("passed"):
        return "True"
    failed = [_FLAG_NAMES.get(k, k)
              for k, ok in quality.get("checks", {}).items() if not ok]
    return ";".join(failed) if failed else f"grade_{quality.get('grade', '?')}"


def _publish_flag(decision: dict) -> str:
    if decision.get("publish"):
        return "True"
    tags = decision.get("reason_tags")
    return ";".join(tags) if tags else "unpublished"


COLUMNS = [
    "PublicID", "Date", "Latitude", "Longitude",
    "strike1", "dip1", "rake1", "strike2", "dip2", "rake2",
    "GeoNet_M", "GeoNet_depth", "GeoNet_depth_unc", "Mw", "Depth",
    "Depth_lo", "Depth_hi", "dZ_GeoNet", "Depth_DCmax", "Mo",
    "NS", "AzGap", "Grade", "DC", "CLVD", "VR", "MinStaVR",
    "Jk_n", "Jk_Mw_std", "Jk_DC_std", "Jk_rot_deg",
    "PredDisp_cm", "Detectable",
    "Mxx", "Mxy", "Mxz", "Myy", "Myz", "Mzz",
    "Band", "Model", "Status", "Selection", "Code",
    "quality_flag", "publish_flag", "published",
]

NO_SOLUTION_COLUMNS = [
    "PublicID", "Date", "Latitude", "Longitude", "GeoNet_M", "GeoNet_depth",
    "Outcome", "Stage", "Reason", "Best_VR", "Selection",
]


def row_from(s: dict, published_ids: set[str] | None = None) -> dict:
    """Catalogue row for one solution.json (shared with catalogue_human)."""
    ev = s["event"]
    prov = s.get("provenance", {})
    base = {
        "PublicID": ev["public_id"], "Date": ev["origin_time"],
        "Latitude": round(ev["latitude"], 4),
        "Longitude": round(ev["longitude"], 4),
        "GeoNet_M": round(ev["prelim_mag"], 2),
        "GeoNet_depth": round(ev["depth_km"], 1),
        "GeoNet_depth_unc": ev.get("depth_unc_km", ""),
        "Band": s.get("chosen_band", "").replace("band_", ""),
        "Model": prov.get("velocity_model", ""),
        "Status": s.get("status", config.STATUS_SOLVED),
        "Selection": prov.get("selection_version", "v3"),
        "Code": prov.get("code_commit", ""),
        "published": (ev["public_id"] in published_ids
                      if published_ids is not None else ""),
    }
    if not config.is_solved(s):
        return {**base, "NS": 0, "Grade": "X",
                "quality_flag": f"no_coherent_solution:{s['abort']['stage']}",
                "publish_flag": "no_coherent_solution"}
    p, q = s["preferred"], s["quality"]
    mt = p["tensor_dyne_cm"]  # XYZ basis, dyne-cm
    jk = s.get("jackknife", {}) or {}
    lo, hi = (s.get("depth_pick_flags", {}).get("depth_range_km")
              or q.get("depth_range_km") or ("", ""))
    return {
        **base,
        "strike1": round(p["plane1"]["strike"], 1),
        "dip1": round(p["plane1"]["dip"], 1),
        "rake1": round(p["plane1"]["rake"], 1),
        "strike2": round(p["plane2"]["strike"], 1),
        "dip2": round(p["plane2"]["dip"], 1),
        "rake2": round(p["plane2"]["rake"], 1),
        "Mw": round(p["mw"], 2), "Depth": p["depth_km"],
        "Depth_lo": lo, "Depth_hi": hi,
        "dZ_GeoNet": round(p["depth_km"] - ev["depth_km"], 1),
        "Depth_DCmax": s.get("depth_pick_flags", {}).get("dc_max_depth_km", ""),
        "Mo": f"{p['m0_dyne_cm']:.3e}",
        "NS": q["n_stations_used"], "AzGap": q.get("azimuthal_gap_deg", ""),
        "Grade": q.get("grade", ""), "DC": round(p["pdc"], 1),
        "CLVD": round(p["pclvd"], 1), "VR": round(p["vr"], 1),
        "MinStaVR": q.get("min_own_vr", ""),
        "Jk_n": jk.get("n_subsets", ""), "Jk_Mw_std": jk.get("mw_std", ""),
        "Jk_DC_std": jk.get("dc_std", ""),
        "Jk_rot_deg": jk.get("max_tensor_rotation_deg", ""),
        "PredDisp_cm": (round(s["forward_model"]["peak_abs_cm"], 3)
                        if s.get("forward_model") else ""),
        "Detectable": s.get("forward_model", {}).get("detectable", ""),
        **{k.capitalize().replace("m", "M", 1): round(v / 1e20, 4)
           for k, v in mt.items()},
        "quality_flag": _quality_flag(q),
        "publish_flag": _publish_flag(s.get("publish_decision", {})),
    }


def _num(x, nd: int = 1):
    try:
        return round(float(x), nd)
    except (TypeError, ValueError):
        return None


def _station_key(r: dict) -> str:
    st = str(r.get("station", ""))
    return st if "." in st else f"{r.get('network', '')}.{st}"


def ledger_line(s: dict) -> dict:
    """One JSON line per event: the event's headline numbers and a
    per-station dictionary of numeric values, used or not. Keys per station:
    dist_km, az, sector, snr_Z/R/T, snr_med, amp_ratio, window_s, used (1/0),
    drop (reason class when not used), own_vr, shift_s. Absent values are
    omitted rather than written as blanks."""
    import invert
    ev = s["event"]
    solved = config.is_solved(s)
    p = s.get("preferred", {})

    def _entry(r: dict, used: bool) -> dict:
        snr = r.get("snr")
        if not isinstance(snr, dict):   # v3/v4 rows stored a single number
            snr = {}
        vals = {
            "dist_km": _num(r.get("distance_km")),
            "az": _num(r.get("azimuth")),
            "sector": r.get("sector"),
            "snr_Z": _num(snr.get("Z")), "snr_R": _num(snr.get("R")),
            "snr_T": _num(snr.get("T")), "snr_med": _num(r.get("snr_med")),
            "amp_ratio": _num(r.get("amp_ratio"), 2),
            "window_s": _num(r.get("window_end_s"), 0),
            "used": int(used),
            "drop": None if used else invert.reason_class(r.get("reason", "")),
            "own_vr": _num(r.get("station_vr")),
            "shift_s": _num(r.get("zcor_s")),
        }
        return {k: v for k, v in vals.items() if v is not None}

    stations: dict[str, dict] = {}
    for r in s.get("stations_used", []):
        stations[_station_key(r)] = _entry(r, True)
    for r in s.get("stations_dropped", []):
        if "network" in r or "." in str(r.get("station", "")):
            stations.setdefault(_station_key(r), _entry(r, False))
    return {
        "PublicID": ev["public_id"], "Date": ev["origin_time"][:10],
        "Band": s.get("chosen_band", "").replace("band_", ""),
        "Depth": p.get("depth_km"),
        "Mw": round(p["mw"], 2) if solved else None,
        "VR": round(p["vr"], 1) if solved else None,
        "Grade": s.get("quality", {}).get("grade", "X"),
        "Selection": s.get("provenance", {}).get("selection_version", ""),
        "stations": stations,
    }


def no_solution_row(s: dict) -> dict:
    """not_published.csv row for an attempted event with no coherent solution."""
    ev, ab = s["event"], s["abort"]
    return {
        "PublicID": ev["public_id"], "Date": ev["origin_time"],
        "Latitude": round(ev["latitude"], 4),
        "Longitude": round(ev["longitude"], 4),
        "GeoNet_M": round(ev["prelim_mag"], 2),
        "GeoNet_depth": round(ev["depth_km"], 1),
        "Outcome": "no_coherent_solution", "Stage": ab.get("stage", ""),
        "Reason": ab.get("reason", ""), "Best_VR": ab.get("best_vr", ""),
        "Selection": s.get("provenance", {}).get("selection_version", ""),
    }


def skipped_row(pid: str, rec: dict) -> dict:
    """not_published.csv row for an event seen but not attempted (below the
    processing floor), from the state file's "skipped" record."""
    return {
        "PublicID": pid, "Date": rec.get("origin_time", ""),
        "Latitude": rec.get("latitude", ""), "Longitude": rec.get("longitude", ""),
        "GeoNet_M": rec.get("prelim_mag", ""), "GeoNet_depth": rec.get("depth_km", ""),
        "Outcome": "not_attempted", "Stage": "processing floor",
        "Reason": rec.get("reason", ""), "Best_VR": "", "Selection": "",
    }


def tables_dir(events_dir: Path | None = None) -> Path:
    """Where the CSV tables live for this archive."""
    events_dir = Path(events_dir or config.EVENTS_DIR).resolve()
    return (config.REPO_DIR if events_dir == (config.REPO_DIR / "events").resolve()
            else events_dir)


def build_catalogue(events_dir: Path | None = None) -> Path | None:
    """Scan <events_dir>/*/solution.json -> the three tables."""
    events_dir = events_dir or config.EVENTS_DIR
    out_dir = tables_dir(events_dir)
    solutions = sorted(events_dir.glob("*/solution.json"))
    if not solutions:
        return None
    published_ids: set[str] = set()
    state: dict = {}
    if config.STATE_FILE.exists():
        state = json.loads(config.STATE_FILE.read_text())
        published_ids = {p["public_id"] for p in state.get("published", [])}

    rows, no_solution, ledger = [], [], []
    for path in solutions:
        s = json.loads(path.read_text())
        rows.append(row_from(s, published_ids))
        ledger.append(ledger_line(s))
        if not config.is_solved(s):
            no_solution.append(no_solution_row(s))
    # events seen but never attempted (below the floor); an event that was
    # attempted later is listed once, from its solution.json
    attempted = {r["PublicID"] for r in rows}
    no_solution += [skipped_row(pid, rec)
                    for pid, rec in state.get("skipped", {}).items()
                    if pid not in attempted]

    rows.sort(key=lambda r: r["Date"])
    out = out_dir / "catalogue.csv"
    _write(out, COLUMNS, rows)
    no_solution.sort(key=lambda r: str(r["Date"]))
    _write(out_dir / "not_published.csv", NO_SOLUTION_COLUMNS, no_solution)
    ledger.sort(key=lambda r: (r["Date"], r["PublicID"]))
    _write_jsonl(out_dir / "station_ledger.jsonl", ledger)
    n_st = sum(len(r["stations"]) for r in ledger)
    print(f"catalogue: {len(rows)} events ({len(no_solution)} with no "
          f"solution, {n_st} station entries) -> {out_dir}")
    return out


def _write(path: Path, columns: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in columns})


def _write_jsonl(path: Path, rows: list[dict]) -> None:
    with open(path, "w") as f:
        for row in rows:
            f.write(json.dumps(row, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    build_catalogue()
