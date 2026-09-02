"""Shared matplotlib/cartopy map styling for all outward-facing figures.

House rules: cartopy geographic axes with OCEAN, LAND, COASTLINE features
and labelled gridlines; cmcrameri 'vik' for diverging fields; fonts sized
for an A4 page or a slide. No PyGMT — matplotlib keeps the toolchain to one
plotting stack and matches the mttime waveform figures.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 11,
    "figure.dpi": 150,
})


def event_region(event: dict, stations: list[dict], pad_deg: float = 0.4) -> list:
    """[lon0, lon1, lat0, lat1] covering epicentre + stations, padded, with
    a minimum 1.5 deg span so sparse geometries stay plottable."""
    lons = [r["longitude"] for r in stations] + [event["longitude"]]
    lats = [r["latitude"] for r in stations] + [event["latitude"]]
    lon0, lon1 = min(lons) - pad_deg, max(lons) + pad_deg
    lat0, lat1 = min(lats) - pad_deg, max(lats) + pad_deg
    if lon1 - lon0 < 1.5:
        mid = 0.5 * (lon0 + lon1)
        lon0, lon1 = mid - 0.75, mid + 0.75
    if lat1 - lat0 < 1.5:
        mid = 0.5 * (lat0 + lat1)
        lat0, lat1 = mid - 0.75, mid + 0.75
    return [lon0, lon1, lat0, lat1]


def square_region(region: list) -> list:
    """Expand the shorter side so the region draws square in Mercator
    (width ~ dlon, height ~ dlat/cos(mid_lat)) — used so panels with
    different content still render at identical drawn size."""
    lon0, lon1, lat0, lat1 = region
    mid = 0.5 * (lat0 + lat1)
    w = lon1 - lon0
    h = (lat1 - lat0) / np.cos(np.radians(mid))
    if w > h:
        extra = (w - h) * np.cos(np.radians(mid)) / 2.0
        lat0, lat1 = lat0 - extra, lat1 + extra
    else:
        extra = (h - w) / 2.0
        lon0, lon1 = lon0 - extra, lon1 + extra
    return [lon0, lon1, lat0, lat1]


def geo_axes(fig, rect, region, labels=True):
    """Cartopy Mercator axes at an explicit figure rect [x, y, w, h]
    (manual placement: cartopy's fixed aspect fights layout engines)."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    proj = ccrs.Mercator(central_longitude=0.5 * (region[0] + region[1]))
    ax = fig.add_axes(rect, projection=proj)
    ax.set_extent(region, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#ffffff",
                   zorder=0)
    ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="0.85",
                   zorder=0)
    ax.add_feature(cfeature.LAKES.with_scale("10m"), facecolor="#ffffff",
                   edgecolor="0.25", linewidth=0.3, zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.6,
                   edgecolor="0.15", zorder=2)
    gl = ax.gridlines(draw_labels=labels, linewidth=0.3, color="0.7",
                      alpha=0.6, linestyle="--")
    if labels:
        gl.top_labels = False
        gl.right_labels = False
    return ax


def vik():
    """cmcrameri vik: scientifically-derived, colourblind-safe diverging."""
    from cmcrameri import cm

    return cm.vik


def panel_label(ax, text: str) -> None:
    """Label inside the panel's top-left corner — cartopy's aspect handling
    makes outside titles unreliable with manual axes placement."""
    ax.text(0.02, 0.985, text, transform=ax.transAxes, va="top", ha="left",
            fontsize=9, zorder=20,
            bbox=dict(facecolor="white", edgecolor="0.4", linewidth=0.5,
                      boxstyle="square,pad=0.25"))


def add_beachball(ax, lon, lat, tensor_rtp: dict, width_frac: float = 0.16,
                  facecolor: str = "firebrick") -> None:
    """Full deviatoric moment tensor beachball (obspy/mopad rendering — the
    same as the mttime waveform figures) at a map position. Uses the actual
    tensor elements so any CLVD is represented, not just the closest DC.
    Position/width are computed in the axes' projected coordinates."""
    import cartopy.crs as ccrs
    from obspy.imaging.beachball import beach

    x, y = ax.projection.transform_point(lon, lat, ccrs.PlateCarree())
    x0, x1, _, _ = ax.get_extent(crs=ax.projection)
    mt = [tensor_rtp["MRR"], tensor_rtp["MTT"], tensor_rtp["MPP"],
          tensor_rtp["MRT"], tensor_rtp["MRP"], tensor_rtp["MTP"]]
    ball = beach(mt, xy=(x, y), width=width_frac * (x1 - x0), linewidth=0.6,
                 facecolor=facecolor, zorder=10)
    ax.add_collection(ball)


# ---------------------------------------------------------------------------
# Bundled context datasets (see data/README.md for provenance)
# ---------------------------------------------------------------------------
_FAULTS_CACHE = None


def load_faults():
    """NZ Active Faults (NZAFD 250K, GNS) as a list of (lons, lats) arrays."""
    global _FAULTS_CACHE
    if _FAULTS_CACHE is None:
        import json
        from pathlib import Path

        path = Path(__file__).parent / "data" / "nz_active_faults.geojson"
        segments = []
        for feat in json.loads(path.read_text())["features"]:
            geom = feat["geometry"]
            lines = (geom["coordinates"] if geom["type"] == "MultiLineString"
                     else [geom["coordinates"]])
            for line in lines:
                arr = np.asarray(line, dtype=float)
                segments.append((arr[:, 0], arr[:, 1]))
        _FAULTS_CACHE = segments
    return _FAULTS_CACHE


def load_gnss_marks():
    """Operating GeoNet GNSS marks (delta snapshot) as a DataFrame."""
    import pandas as pd
    from pathlib import Path

    df = pd.read_csv(Path(__file__).parent / "data" / "gnss_marks.csv")
    return df[df["End Date"].str.startswith("9999")]


def scale_bar(ax, region, ccrs):
    """Simple km scale bar, bottom-left."""
    lon0, lon1, lat0, lat1 = region
    mid_lat = 0.5 * (lat0 + lat1)
    span_km = (lon1 - lon0) * 111.32 * np.cos(np.radians(mid_lat))
    km = min([1, 2, 5, 10, 20, 50, 100, 200],
             key=lambda k: abs(k - 0.25 * span_km))
    dlon = km / (111.32 * np.cos(np.radians(mid_lat)))
    x0 = lon0 + 0.05 * (lon1 - lon0)
    y = lat0 + 0.06 * (lat1 - lat0)
    ax.plot([x0, x0 + dlon], [y, y], "-", color="black", linewidth=3,
            transform=ccrs.PlateCarree(), zorder=15,
            solid_capstyle="butt")
    ax.text(x0 + 0.5 * dlon, y + 0.015 * (lat1 - lat0), f"{km} km",
            ha="center", va="bottom", fontsize=8,
            transform=ccrs.PlateCarree(), zorder=15)


def draw_context(ax, region, ccrs, gnss=True, gnss_labels=False):
    """Active faults (+ optionally operating GNSS marks) in the region."""
    lon0, lon1, lat0, lat1 = region
    for lons, lats in load_faults():
        if (lons.max() < lon0 or lons.min() > lon1
                or lats.max() < lat0 or lats.min() > lat1):
            continue
        ax.plot(lons, lats, "-", color="#8B3A3A", linewidth=0.6, alpha=0.8,
                transform=ccrs.PlateCarree(), zorder=4)
    if not gnss:
        return 0
    marks = load_gnss_marks()
    inside = marks[
        (marks.Longitude >= lon0) & (marks.Longitude <= lon1)
        & (marks.Latitude >= lat0) & (marks.Latitude <= lat1)
    ]
    ax.plot(inside.Longitude, inside.Latitude, "s", color="#0173B2",
            markeredgecolor="black", markeredgewidth=0.3, markersize=5,
            transform=ccrs.PlateCarree(), zorder=6)
    if gnss_labels:
        for _, m in inside.iterrows():
            ax.annotate(m.Mark, (m.Longitude, m.Latitude), xytext=(3, 3),
                        textcoords="offset points", fontsize=6,
                        color="#0173B2")
    return len(inside)


def nodal_plane_arc(strike: float, dip: float, n: int = 91):
    """Lower-hemisphere equal-area trace of one nodal plane, unit radius
    (matches the obspy beachball projection r = sqrt(2) sin(takeoff/2))."""
    phi, delta = np.radians(strike), np.radians(dip)
    t = np.linspace(0.0, np.pi, n)
    sv = np.array([np.sin(phi), np.cos(phi), 0.0])
    dv = np.array([np.cos(phi) * np.cos(delta), -np.sin(phi) * np.cos(delta),
                   -np.sin(delta)])
    v = np.outer(np.cos(t), sv) + np.outer(np.sin(t), dv)
    takeoff = np.arccos(np.clip(-v[:, 2], -1.0, 1.0))
    az = np.arctan2(v[:, 0], v[:, 1])
    r = np.sqrt(2.0) * np.sin(takeoff / 2.0)
    return r * np.sin(az), r * np.cos(az)


def inset_dc_ball(fig, rect, plane1: dict, plane2: dict,
                  highlight_color: str = "#0173B2"):
    """Small DC beachball inset with the modelled plane (plane1) traced in
    a highlight colour over the standard nodal lines."""
    from obspy.imaging.beachball import beach

    ax = fig.add_axes(rect)
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.45, 1.25)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_collection(beach(
        [plane1["strike"], plane1["dip"], plane1["rake"]], xy=(0, 0),
        width=2.0, linewidth=0.8, facecolor="firebrick"))
    x, y = nodal_plane_arc(plane1["strike"], plane1["dip"])
    ax.plot(x, y, "-", color=highlight_color, linewidth=2.2, zorder=20,
            solid_capstyle="round")
    ax.text(0, -1.38, "modelled plane", ha="center", va="center", fontsize=7,
            color=highlight_color)
    return ax
