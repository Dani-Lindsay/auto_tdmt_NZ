"""Outward figures (matplotlib/cartopy — one plotting stack, matching the
mttime waveform-fit figures):

make_share_figure: two same-region map panels — (a) stations + the full
deviatoric moment tensor beachball, (b) Okada-predicted vertical surface
displacement (cmcrameri 'vik') — with a NISAR pass table + provenance
footer.

compose_outward: stacks mttime's waveform-fit page(s) (DC+CLVD
decomposition, per-station fits) above the map panels into the single
figure the email carries.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import map_style
import okada_forward


def _local_km_to_geo(x_km, y_km, lon0: float, lat0: float):
    """Small-region approximation: local east/north km -> lon/lat."""
    lat = lat0 + np.asarray(y_km) / 110.57
    lon = lon0 + np.asarray(x_km) / (111.32 * np.cos(np.radians(lat0)))
    return lon, lat


def make_share_figure(
    solution: dict, forward: dict, passes: list[dict], out_path: Path
) -> Path:
    import cartopy.crs as ccrs

    ev = solution["event"]
    pref = solution["preferred"]
    lon0, lat0 = ev["longitude"], ev["latitude"]
    roi = map_style.event_region(ev, solution["stations_used"], pad_deg=0.6)

    fig = plt.figure(figsize=(15.0, 6.4))

    # ---- panel a: stations + full-MT beachball over the station ROI -------
    ax = map_style.geo_axes(fig, [0.035, 0.30, 0.24, 0.60], roi)
    ax.plot(
        [r["longitude"] for r in solution["stations_used"]],
        [r["latitude"] for r in solution["stations_used"]],
        "^", color="forestgreen", markeredgecolor="black",
        markeredgewidth=0.4, markersize=8, transform=ccrs.PlateCarree(),
    )
    for r in solution["stations_used"]:
        ax.annotate(
            r["station"], (r["longitude"], r["latitude"]),
            xytext=(4, -4), textcoords="offset points", fontsize=6.5,
            xycoords=ccrs.PlateCarree()._as_mpl_transform(ax),
        )
    map_style.add_beachball(ax, lon0, lat0, pref["tensor_rtp_dyne_cm"])
    map_style.panel_label(
        ax, f"(a) {ev['public_id']}  Mw {pref['mw']:.1f}  "
            f"depth {pref['depth_km']:g} km")

    # ---- panels b-d: E/N/U predicted displacement over the modelled area --
    fw = forward["plane1"]
    fault = fw["fault"]
    outline = okada_forward.fault_outline(fault)
    olon, olat = _local_km_to_geo(
        np.array(outline["outline_x_km"]), np.array(outline["outline_y_km"]),
        lon0, lat0)
    tlon, tlat = _local_km_to_geo(
        np.array(outline["top_x_km"]), np.array(outline["top_y_km"]),
        lon0, lat0)
    glon, glat = _local_km_to_geo(fw["x_km"], fw["y_km"], lon0, lat0)
    model_region = [float(glon.min()), float(glon.max()),
                    float(glat.min()), float(glat.max())]
    comps = [("(b) east", fw["ue_m"]), ("(c) north", fw["un_m"]),
             ("(d) up", fw["uz_m"])]
    vmax = max(
        0.1, max(float(np.abs(u).max()) for _, u in comps) * 100.0)

    pm = None
    for i, (name, u_m) in enumerate(comps):
        u_cm = u_m * 100.0
        u_plot = np.ma.masked_where(np.abs(u_cm) < 0.02 * vmax, u_cm)
        axi = map_style.geo_axes(
            fig, [0.315 + 0.215 * i, 0.30, 0.19, 0.60], model_region,
            labels=(i == 0),
        )
        pm = axi.pcolormesh(
            glon, glat, u_plot, cmap=map_style.vik(), vmin=-vmax, vmax=vmax,
            transform=ccrs.PlateCarree(), alpha=0.9, shading="auto", zorder=3,
        )
        # surface projection of the modelled plane (bold edge = up-dip)
        axi.plot(olon, olat, "--", color="black", linewidth=0.9,
                 transform=ccrs.PlateCarree(), zorder=10)
        axi.plot(tlon, tlat, "-", color="black", linewidth=2.2,
                 transform=ccrs.PlateCarree(), zorder=10)
        axi.plot(lon0, lat0, "*", color="yellow", markeredgecolor="black",
                 markersize=13, transform=ccrs.PlateCarree(), zorder=11)
        map_style.panel_label(
            axi, f"{name}  (peak {np.abs(u_cm).max():.2f} cm)")

    cax = fig.add_axes([0.955, 0.34, 0.011, 0.52])
    cb = fig.colorbar(pm, cax=cax, orientation="vertical")
    cb.set_label("predicted displacement [cm]", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    peak_cm = forward["peak_abs_m"] * 100.0
    detect = ("potentially InSAR detectable" if forward["detectable"]
              else "below InSAR detection")

    # ---- footer: NISAR passes + provenance --------------------------------
    prov = solution["provenance"]
    lines = [
        f"Okada forward model on NODAL PLANE 1 "
        f"(strike/dip/rake {fault['strike']:.0f}/{fault['dip']:.0f}/"
        f"{fault['rake']:.0f}, dashed outline, bold edge = up-dip): "
        f"peak |u| {peak_cm:.2f} cm - {detect}"]
    if passes:
        lines.append("NISAR passes at epicentre (last | predicted next):  "
                     + "   ".join(
                         f"track {p['track']:03d}{p['direction']}: "
                         f"{p['last_utc']} | {p['next_utc']}"
                         for p in passes[:3]))
    else:
        lines.append("No NISAR coverage found at the epicentre yet.")
    lines.append(
        f"PRELIMINARY automated solution - model {prov['velocity_model']}, "
        f"mttime {prov['mttime_version']}, VR {pref['vr']:.0f}%, "
        f"DC {pref['pdc']:.0f}% / CLVD {pref['pclvd']:.0f}%"
    )
    for i, line in enumerate(lines):
        fig.text(0.035, 0.175 - 0.055 * i, line, fontsize=8.5,
                 family="monospace", va="top")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def compose_outward(
    bbwaves_paths: list[Path], maps_path: Path, out_path: Path
) -> Path:
    """Stack the station-map/forward-model panel above mttime's
    waveform-fit page(s) (full DC+CLVD decomposition and per-station fits)
    into one outward figure for the email."""
    from PIL import Image

    assert bbwaves_paths, "no waveform-fit figures to compose"
    # maps on top, waveform modelling below
    images = [Image.open(p) for p in [maps_path, *bbwaves_paths]]
    width = max(im.width for im in images)
    scaled = [
        im.resize((width, int(im.height * width / im.width)))
        if im.width != width else im
        for im in images
    ]
    total_h = sum(im.height for im in scaled)
    canvas = Image.new("RGB", (width, total_h), "white")
    y = 0
    for im in scaled:
        canvas.paste(im, (0, y))
        y += im.height
    canvas.save(out_path, quality=90)
    return out_path
