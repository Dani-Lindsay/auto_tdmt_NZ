"""Shareable summary figure (email/slide friendly): location map with
focal mechanism + predicted surface displacement + NISAR pass table.

Single JPEG, PyGMT, fonts sized for an A4 page / slide.
"""

from __future__ import annotations

import numpy as np
from pathlib import Path

import config


def _local_km_to_geo(x_km, y_km, lon0: float, lat0: float):
    """Small-region approximation: local east/north km -> lon/lat."""
    lat = lat0 + np.asarray(y_km) / 110.57
    lon = lon0 + np.asarray(x_km) / (111.32 * np.cos(np.radians(lat0)))
    return lon, lat


def make_share_figure(
    solution: dict, forward: dict, passes: list[dict], out_path: Path
) -> Path:
    import pygmt
    import xarray as xr

    ev = solution["event"]
    pref = solution["preferred"]
    p1 = pref["plane1"]
    lon0, lat0 = ev["longitude"], ev["latitude"]

    fig = pygmt.Figure()
    pygmt.config(FONT_ANNOT_PRIMARY="11p", FONT_LABEL="12p", FONT_TITLE="14p")

    # --- left panel: regional map with focal mechanism ---------------------
    pad = 1.6
    region = [lon0 - pad, lon0 + pad, lat0 - pad, lat0 + pad]
    fig.basemap(region=region, projection="M11c", frame=["af", "WSen"])
    fig.coast(land="grey92", water="lightblue", shorelines="1/0.5p,black",
              resolution="i")
    for row in solution["stations_used"]:
        fig.plot(x=row["longitude"], y=row["latitude"], style="t0.35c",
                 fill="darkgreen", pen="0.4p,black")
    fig.meca(
        spec=dict(strike=p1["strike"], dip=p1["dip"], rake=p1["rake"],
                  magnitude=pref["mw"]),
        scale="1.4c", longitude=lon0, latitude=lat0,
        depth=pref["depth_km"], compressionfill="red",
    )
    title = (f"{ev['public_id']}  Mw {pref['mw']:.1f}  "
             f"depth {pref['depth_km']:g} km")
    fig.text(position="TC", text=title, font="13p,Helvetica-Bold,black",
             offset="0/0.6c", no_clip=True)

    # --- right panel: predicted vertical displacement (plane 1) ------------
    fw = forward["plane1"]
    lon_axis, _ = _local_km_to_geo(fw["x_km"], 0.0, lon0, lat0)
    _, lat_axis = _local_km_to_geo(0.0, fw["y_km"], lon0, lat0)
    grid = xr.DataArray(
        fw["uz_m"] * 100.0, dims=("lat", "lon"),
        coords={"lat": lat_axis, "lon": lon_axis},
    )
    vmax = max(0.5, float(np.abs(grid).max()))
    fig.shift_origin(xshift="12.5c")
    region2 = [float(lon_axis.min()), float(lon_axis.max()),
               float(lat_axis.min()), float(lat_axis.max())]
    fig.basemap(region=region2, projection="M11c", frame=["af", "wSEn"])
    pygmt.makecpt(cmap="polar", series=[-vmax, vmax])
    fig.grdimage(grid=grid)
    fig.coast(shorelines="1/0.5p,black", resolution="i")
    fig.plot(x=lon0, y=lat0, style="a0.5c", fill="yellow", pen="0.8p,black")
    fig.colorbar(frame="af+lpredicted vertical displacement [cm]")
    label = ("peak |u| %.1f cm - potentially InSAR detectable"
             % (forward["peak_abs_m"] * 100.0)
             if forward["detectable"]
             else "peak |u| %.2f cm - below InSAR detection"
             % (forward["peak_abs_m"] * 100.0))
    fig.text(position="TC", text=label, font="12p,Helvetica-Bold,black",
             offset="0/0.6c", no_clip=True)

    # --- bottom: NISAR pass table + provenance -----------------------------
    prov = solution["provenance"]
    lines = []
    if passes:
        lines.append("NISAR passes at epicentre (last | predicted next):")
        for p in passes[:4]:
            lines.append(
                f"track {p['track']:03d}{p['direction']}:  "
                f"{p['last_utc']}  |  {p['next_utc']}"
            )
    else:
        lines.append("No NISAR coverage found at the epicentre yet.")
    lines.append(
        f"PRELIMINARY automated solution - model {prov['velocity_model']}, "
        f"mttime {prov['mttime_version']}, VR {pref['vr']:.0f}%, "
        f"DC {pref['pdc']:.0f}%"
    )
    fig.shift_origin(xshift="-12.5c", yshift="-2.6c")
    fig.basemap(region=[0, 24, 0, 3], projection="X24c/2.2c", frame=0)
    for i, line in enumerate(lines):
        fig.text(x=0.2, y=2.6 - 0.55 * i, text=line,
                 font="10p,Courier,black", justify="ML")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path), dpi=200)
    return out_path
