"""Outward figures (matplotlib/cartopy — one plotting stack, matching the
mttime waveform-fit figures):

make_share_figure: two same-region map panels — (a) stations + the full
deviatoric moment tensor beachball, (b) Okada-predicted vertical surface
displacement (cmcrameri 'vik') — with a NISAR pass table + provenance
footer.

plot_depth_sensitivity: seaborn-styled VR/%DC/Mw vs depth summary
(replaces mttime's depth.bbmw plot). The email carries three separate
figures: the untouched mttime waveform fits, the modelling maps, and the
depth-sensitivity summary.
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


def plot_depth_sensitivity(solution: dict, out_path: Path) -> Path:
    """Depth-sensitivity summary (replaces mttime's depth.bbmw plot):
    VR, %DC (+%CLVD) and Mw vs source depth, seaborn-style, with the
    preferred depth and the VR-tolerance band of the selection rule shown.
    """
    import config

    rows = sorted(solution["depth_search"], key=lambda r: r["depth_km"])
    depths = [r["depth_km"] for r in rows]
    vr = [r["vr"] for r in rows]
    pdc = [r["pdc"] for r in rows]
    pclvd = [r["pclvd"] for r in rows]
    mw = [r["mw"] for r in rows]
    pref = solution["preferred"]
    ev = solution["event"]

    # seaborn "deep" palette values
    c_vr, c_dc, c_clvd, c_mw = "#4C72B0", "#DD8452", "#C44E52", "#55A868"

    with plt.style.context("seaborn-v0_8-whitegrid"):
        fig, axes = plt.subplots(
            3, 1, figsize=(7.5, 8.5), sharex=True,
            gridspec_kw={"hspace": 0.12},
        )
        ax_vr, ax_dc, ax_mw = axes

        vr_max = max(vr)
        ax_vr.axhspan(vr_max - config.PREFER_DC_VR_TOLERANCE, vr_max,
                      color=c_vr, alpha=0.12,
                      label=f"within {config.PREFER_DC_VR_TOLERANCE:g}% of "
                            f"VR max (selection window)")
        ax_vr.plot(depths, vr, "-o", color=c_vr, ms=5, lw=1.8)
        ax_vr.set_ylabel("variance reduction [%]")
        ax_vr.legend(loc="lower right", frameon=True, fontsize=9)

        ax_dc.plot(depths, pdc, "-o", color=c_dc, ms=5, lw=1.8, label="DC")
        ax_dc.plot(depths, pclvd, "--o", color=c_clvd, ms=4, lw=1.2,
                   alpha=0.7, label="CLVD")
        ax_dc.set_ylabel("source type [%]")
        ax_dc.set_ylim(-3, 103)
        ax_dc.legend(loc="upper right", frameon=True, fontsize=9)

        ax_mw.plot(depths, mw, "-o", color=c_mw, ms=5, lw=1.8)
        ax_mw.set_ylabel("Mw")
        ax_mw.set_xlabel("source depth [km]")

        for ax in axes:
            ax.axvline(pref["depth_km"], color="0.25", lw=1.2, ls=":",
                       zorder=1)
        ax_vr.annotate(
            f"preferred: {pref['depth_km']:g} km\n"
            f"VR {pref['vr']:.1f}%  DC {pref['pdc']:.0f}%  "
            f"Mw {pref['mw']:.2f}",
            xy=(pref["depth_km"], pref["vr"]),
            xytext=(10, -35), textcoords="offset points", fontsize=9,
            bbox=dict(facecolor="white", edgecolor="0.6",
                      boxstyle="round,pad=0.35"),
        )
        band = solution.get("chosen_band", "").replace("band_", "")
        ax_vr.set_title(
            f"{ev['public_id']}  depth sensitivity  ({band}, "
            f"{solution['quality']['n_stations_used']} stations)",
            fontsize=11,
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path
