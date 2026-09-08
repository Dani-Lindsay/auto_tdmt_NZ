"""Spider plot of an earthquake sequence: epicentres on a local map, every
solved event's deviatoric beachball pushed out to a ring around the
cluster (time order, clockwise from the top) with a leader line back to
its epicentre, labelled with Mw and date.

    pixi run python tools/plot_sequence.py --start 2026p660160
    pixi run python tools/plot_sequence.py --start 2026p660160 --days 30 \
        --km 40 --out docs/examples/fiordland_2026-09_sequence.jpg

Events are taken from catalogue.csv: everything from the start event's
origin time onward (``--days``) within ``--km`` of it. Grey land, blue
sea and the NZ Active Faults Database, as the README map.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import config  # noqa: E402
import map_style  # noqa: E402


def sequence(start: str, days: float, km: float) -> pd.DataFrame:
    cat = pd.read_csv(ROOT / "catalogue.csv")
    cat["Date"] = pd.to_datetime(cat.Date, utc=True)
    m = cat[cat.PublicID == start]
    assert len(m), f"{start} is not in catalogue.csv"
    m = m.iloc[0]
    dlat = (cat.Latitude - m.Latitude) * 111.2
    dlon = (cat.Longitude - m.Longitude) * 111.2 * math.cos(math.radians(m.Latitude))
    sel = cat[(cat.Date >= m.Date)
              & (cat.Date <= m.Date + pd.Timedelta(days=days))
              & (np.hypot(dlat, dlon) <= km)]
    return sel.sort_values("Date").reset_index(drop=True)


def tensor_row(pid: str) -> dict:
    """mt-convention row for pygmt.meca from the archived solution."""
    d = config.find_event_dir(pid, ROOT / "events")
    s = json.loads((d / "solution.json").read_text())
    rtp = s["preferred"]["tensor_rtp_dyne_cm"]
    m = max(abs(v) for v in rtp.values()) or 1.0
    exp = int(math.floor(math.log10(m)))
    k = 10.0 ** exp
    return {"mrr": rtp["MRR"] / k, "mtt": rtp["MTT"] / k, "mff": rtp["MPP"] / k,
            "mrt": rtp["MRT"] / k, "mrf": rtp["MRP"] / k, "mtf": rtp["MTP"] / k,
            "exponent": exp}


def main(start: str, days: float, km: float, out: Path) -> Path:
    import pygmt

    ev = sequence(start, days, km)
    n = len(ev)
    clat, clon = ev.Latitude.mean(), ev.Longitude.mean()
    cosl = math.cos(math.radians(clat))
    # a square-ish 1-degree frame with the cluster in the middle
    half_lat = 0.36
    half_lon = half_lat / cosl
    region = [clon - half_lon, clon + half_lon, clat - half_lat, clat + half_lat]
    ring_lat = 0.2                      # ring radius, degrees of latitude
    ring_lon = ring_lat / cosl

    fig = pygmt.Figure()
    pygmt.config(FONT="10p", FONT_TITLE="11p,Helvetica", MAP_FRAME_TYPE="plain",
                 FORMAT_GEO_MAP="ddd.xF", MAP_GRID_PEN_PRIMARY="0.25p,gray60,.")
    t0, t1 = ev.Date.iloc[0], ev.Date.iloc[-1]
    # grey land, blue sea and dark-red faults: the house style of the
    # README map and the per-event figures
    sea, land = "#dce9f5", "#e4e4e4"
    fig.basemap(region=region, projection="M14c",
                frame=[f"WSen+t{start} sequence, {t0:%d %b} to {t1:%d %b %Y}",
                       "xa0.25f0.05g0.25", "ya0.25f0.05g0.25"])
    fig.coast(land=land, water=sea, lakes=sea, shorelines="0.4p,gray25",
              resolution="f", map_scale="jBR+w20k+o0.5c/0.5c+f+u")
    xs, ys = [], []
    for lons, lats in map_style.load_faults():
        if (lons.max() < region[0] or lons.min() > region[1]
                or lats.max() < region[2] or lats.min() > region[3]):
            continue
        xs += list(lons) + [np.nan]
        ys += list(lats) + [np.nan]
    if xs:
        fig.plot(x=np.array(xs), y=np.array(ys), pen="0.7p,#8B3A3A")

    # ring positions: clockwise from the top, in time order
    angles = [90.0 - 360.0 * i / n for i in range(n)]
    ball_lon = [clon + ring_lon * math.cos(math.radians(a)) for a in angles]
    ball_lat = [clat + ring_lat * math.sin(math.radians(a)) for a in angles]

    for i, r in ev.iterrows():                 # leader lines, then epicentres
        fig.plot(x=[r.Longitude, ball_lon[i]], y=[r.Latitude, ball_lat[i]],
                 pen="0.9p,black")
    for i, r in ev.iterrows():
        fig.plot(x=[r.Longitude], y=[r.Latitude],
                 style=f"c{0.10 + 0.05 * max(0.0, r.Mw - 3.5):.2f}c",
                 fill="black", pen="0.4p,white")

    # beachballs on the ring: red compression, size by Mw, faded for C/D
    opacity = {"A": 0, "B": 0, "C": 35, "D": 55}
    for i, r in ev.iterrows():
        spec = {**tensor_row(r.PublicID), "longitude": r.Longitude,
                "latitude": r.Latitude, "depth": r.Depth}
        fig.meca(spec=spec, convention="mt", component="deviatoric",
                 scale="1.0c", plot_longitude=ball_lon[i],
                 plot_latitude=ball_lat[i], compression_fill="#d94a2b",
                 extension_fill="white", pen="0.5p,black",
                 transparency=opacity.get(r.Grade, 55))
        # magnitude and date above every ball, offset in plot units so it
        # sits just clear of the ball whatever the map scale
        fig.text(x=ball_lon[i], y=ball_lat[i], offset="0/0.62c",
                 text=f"M@-w@- {r.Mw:.1f}, {r.Date:%d %b}",
                 font="8.5p,Helvetica-Bold", justify="BC", no_clip=True)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out), dpi=200)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--start", required=True, help="publicID of the first event")
    ap.add_argument("--days", type=float, default=30.0)
    ap.add_argument("--km", type=float, default=40.0)
    ap.add_argument("--out", default=None, help="output image path")
    a = ap.parse_args()
    out = Path(a.out) if a.out else ROOT / "docs" / "examples" / f"{a.start}_sequence.jpg"
    print(main(a.start, a.days, a.km, out))
