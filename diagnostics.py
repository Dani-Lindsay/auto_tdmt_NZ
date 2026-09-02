"""QC figures for every intermediate stage of event processing, so a human
can verify the data before trusting any solution.

Written into <event_dir>/diagnostics/:
  01_raw_counts.jpg            raw downloaded traces (counts, as delivered)
  02_displacement.jpg          response-removed displacement, unfiltered
  03_final_zrt.jpg             filtered/decimated ZRT record sections that
                               enter the inversion (with GF band annotated)
  04_station_map.jpg           quick station/epicentre geometry check
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy import Stream, UTCDateTime

import config

# A4-landscape / slide-friendly figures: readable fonts at full-page size
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "figure.dpi": 150})
A4_LANDSCAPE = (11.7, 8.3)


def _record_section(
    st: Stream, origin: UTCDateTime, ax, title: str, xmax: float
) -> None:
    for tr in st:
        t = tr.times() - (origin - tr.stats.starttime)
        dist_km = tr.stats.distance / 1000.0
        data = tr.data.astype(float)
        peak = np.max(np.abs(data))
        if peak > 0:
            data = data / peak * 8.0  # +/- 8 km plot amplitude
        ax.plot(t, data + dist_km, color="black", linewidth=0.5)
        ax.text(
            xmax, dist_km, f"{tr.stats.station}.{tr.stats.channel}",
            fontsize=8, va="bottom", ha="right",
        )
    ax.set_xlim(-config.TIME_BEFORE_S, xmax)
    ax.set_xlabel("Time from origin [s]")
    ax.set_ylabel("Distance [km]")
    ax.set_title(title, fontsize=11)


def plot_stages(
    stages: dict, event: dict, band_hz: tuple[float, float], out_dir: Path
) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    origin = UTCDateTime(event["origin_time"])
    written = []

    per_component_figs = [
        ("raw", "01_raw_counts.jpg", "Raw (counts, as downloaded)", None),
        ("displacement", "02_displacement.jpg",
         "Displacement, response removed, unfiltered [m]", None),
        ("final", "03_final_zrt.jpg",
         f"Inversion input: {1/band_hz[1]:.0f}-{1/band_hz[0]:.0f} s, "
         f"{config.DT:g} sps, cm", "ZRT"),
    ]
    for key, fname, title, comps in per_component_figs:
        st = stages.get(key, Stream())
        if len(st) == 0:
            continue
        components = list(comps) if comps else sorted(
            {tr.stats.channel[-1] for tr in st}
        )
        fig, axes = plt.subplots(
            1, len(components), figsize=A4_LANDSCAPE, squeeze=False
        )
        for comp, ax in zip(components, axes[0]):
            sel = st.select(component=comp)
            _record_section(
                sel, origin, ax, f"{title}\ncomponent {comp}",
                xmax=config.TIME_AFTER_S,
            )
        fig.suptitle(
            f"{event['public_id']} M{event['prelim_mag']:.1f} "
            f"{event['locality']}", fontsize=13,
        )
        fig.tight_layout()
        path = out_dir / fname
        fig.savefig(path, dpi=150)
        plt.close(fig)
        written.append(path)

    return written


def plot_station_map(
    event: dict, used: list[dict], dropped_ids: list[str], out_dir: Path
) -> Path:
    """PyGMT map with coastline: epicentre + stations used."""
    import pygmt

    out_dir.mkdir(parents=True, exist_ok=True)
    lons = [r["longitude"] for r in used] + [event["longitude"]]
    lats = [r["latitude"] for r in used] + [event["latitude"]]
    pad = 0.6
    region = [
        min(lons) - pad, max(lons) + pad, min(lats) - pad, max(lats) + pad,
    ]
    fig = pygmt.Figure()
    fig.basemap(region=region, projection="M16c", frame=["af", "WSen"])
    fig.coast(
        land="grey92", water="lightblue", shorelines="1/0.5p,black",
        resolution="i",
    )
    fig.plot(
        x=[r["longitude"] for r in used], y=[r["latitude"] for r in used],
        style="t0.45c", fill="darkgreen", pen="0.5p,black",
    )
    for r in used:
        fig.text(
            x=r["longitude"], y=r["latitude"], text=r["station"],
            font="9p,Helvetica-Bold,black", justify="TL", offset="0.2c/-0.2c",
        )
    fig.plot(
        x=event["longitude"], y=event["latitude"], style="a0.8c",
        fill="red", pen="1p,black",
    )
    title = (
        f"{event['public_id']} M{event['prelim_mag']:.1f} "
        f"{len(used)} stations used, {len(dropped_ids)} dropped"
    )
    fig.text(
        position="TC", text=title, font="12p,Helvetica-Bold,black",
        offset="0/0.5c", no_clip=True,
    )
    path = out_dir / "04_station_map.jpg"
    fig.savefig(str(path), dpi=150)
    return path
