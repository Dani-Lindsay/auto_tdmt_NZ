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
    """Cartopy station-geometry QC map in the house style."""
    import cartopy.crs as ccrs

    import map_style

    out_dir.mkdir(parents=True, exist_ok=True)
    region = map_style.event_region(event, used, pad_deg=0.5)

    fig = plt.figure(figsize=(8.5, 8.5))
    ax = map_style.geo_axes(fig, [0.07, 0.05, 0.88, 0.88], region)
    ax.plot(
        [r["longitude"] for r in used], [r["latitude"] for r in used],
        "^", color="forestgreen", markeredgecolor="black",
        markeredgewidth=0.4, markersize=10, transform=ccrs.PlateCarree(),
    )
    for r in used:
        ax.annotate(
            f"{r['station']} ({r['distance_km']:.0f} km)",
            (r["longitude"], r["latitude"]),
            xytext=(5, -5), textcoords="offset points", fontsize=7,
        )
    ax.plot(event["longitude"], event["latitude"], "*", color="firebrick",
            markeredgecolor="black", markersize=20,
            transform=ccrs.PlateCarree())
    ax.set_title(
        f"{event['public_id']} M{event['prelim_mag']:.1f}: "
        f"{len(used)} stations used, {len(dropped_ids)} dropped"
    )
    path = out_dir / "04_station_map.jpg"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
