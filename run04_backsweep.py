"""Backsweep: process every NZ M>=4.0 earthquake in a date range.

    pixi run python run04_backsweep.py --start 2026-01-01

Resumable: events with an existing solution.json are skipped, so the sweep
can be interrupted and relaunched freely. Never emails (processing only).
NZ spans the antimeridian, so the FDSN archive is queried in two longitude
windows and merged.
"""

from __future__ import annotations

import argparse
import time
from datetime import datetime, timezone

import requests

import config
from run02_process import process_event

EVENT_URL = (
    "{base}/fdsnws/event/1/query?starttime={start}&endtime={end}"
    "&minmagnitude={mag}"
    "&minlatitude={lat0}&maxlatitude={lat1}"
    "&minlongitude={lon0}&maxlongitude={lon1}&format=text"
)


def list_events(start: str, end: str = "2100-01-01") -> list[tuple[str, str, float]]:
    """(publicID, origin_time, magnitude) for NZ M>=floor earthquakes
    shallower than the processing depth ceiling."""
    b = config.NZ_BBOX
    windows = [(b["lon_min"], 180.0), (-180.0, b["lon_max"] - 360.0)]
    rows = []
    for lon0, lon1 in windows:
        url = EVENT_URL.format(
            base=config.FDSN_ARCHIVE, start=start, end=end,
            mag=config.PROCESS_MIN_PRELIM_MAG,
            lat0=b["lat_min"], lat1=b["lat_max"], lon0=lon0, lon1=lon1,
        )
        r = requests.get(url, headers={"User-Agent": config.USER_AGENT},
                         timeout=120)
        if r.status_code == 204:  # no content = no events in window
            continue
        r.raise_for_status()
        for line in r.text.splitlines()[1:]:
            f = line.split("|")
            if f[-1].strip() != "earthquake":
                continue
            depth = float(f[4])
            if (depth > config.MAX_PROCESS_DEPTH_KM
                    and depth not in config.PLACEHOLDER_DEPTHS_KM):
                continue
            rows.append((f[0], f[1], float(f[10])))
    rows.sort(key=lambda x: x[1])
    return rows


def main(start: str, end: str = "2100-01-01") -> None:
    events = list_events(start, end)
    print(f"{len(events)} candidate events since {start}", flush=True)
    done = failed = skipped = nosol = 0
    for i, (pid, when, mag) in enumerate(events, 1):
        existing = config.find_event_dir(pid)
        if existing is not None and (existing / "solution.json").exists():
            skipped += 1
            continue
        t0 = time.time()
        try:
            sol = process_event(pid)
            if not config.is_solved(sol):
                print(
                    f"[{i}/{len(events)}] {pid} {when[:10]} M{mag:.1f} -> "
                    f"NO COHERENT SOLUTION ({sol['abort']['stage']}, "
                    f"best VR {sol['abort']['best_vr']:.0f}) "
                    f"({time.time() - t0:.0f}s)", flush=True,
                )
                nosol += 1
            else:
                p = sol["preferred"]
                print(
                    f"[{i}/{len(events)}] {pid} {when[:10]} M{mag:.1f} -> "
                    f"Mw {p['mw']:.2f} depth {p['depth_km']:g} km "
                    f"VR {p['vr']:.0f}% DC {p['pdc']:.0f}% "
                    f"grade {sol['quality']['grade']} "
                    f"({time.time() - t0:.0f}s)", flush=True,
                )
                done += 1
        except Exception as e:  # noqa: BLE001 - sweep must survive bad events
            print(f"[{i}/{len(events)}] {pid} {when[:10]} M{mag:.1f} "
                  f"FAILED: {str(e)[:160]}", flush=True)
            failed += 1
            # a failed event leaves no solution: remove its empty shell
            shell = config.find_event_dir(pid)
            if shell is not None and not (shell / "solution.json").exists():
                import shutil as _sh
                _sh.rmtree(shell, ignore_errors=True)
        time.sleep(2)  # politeness between events
    print(f"sweep complete {datetime.now(timezone.utc):%Y-%m-%dT%H:%MZ}: "
          f"{done} solved, {nosol} no coherent solution, "
          f"{failed} failed, {skipped} already done",
          flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--start", default="2026-01-01")
    ap.add_argument("--end", default="2100-01-01")
    args = ap.parse_args()
    main(args.start, args.end)
