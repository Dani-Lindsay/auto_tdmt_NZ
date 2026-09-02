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
