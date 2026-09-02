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
        axi.plot(olon, olat, "--", color="black", linewidth=1.2,
                 transform=ccrs.PlateCarree(), zorder=10)
        axi.plot(tlon, tlat, "-", color="black", linewidth=2.6,
                 transform=ccrs.PlateCarree(), zorder=10)
        axi.plot(lon0, lat0, "*", color="yellow", markeredgecolor="black",
                 markersize=9, transform=ccrs.PlateCarree(), zorder=11)
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
    """Depth-sensitivity summary in the user's classic 2x3 layout
    (replaces mttime's depth.bbmw plot): VR, %DC, Mw / strike, rake, dip of
    nodal plane 1 vs source depth, black dots, with vertical reference
    lines at the GeoNet depth, the max-VR depth, the max-DC depth and (when
    distinct) the preferred depth from the selection rule.
    """
    rows = sorted(solution["depth_search"], key=lambda r: r["depth_km"])
    depths = [r["depth_km"] for r in rows]
    series = [
        ("Variance Reduction", "VR (%)", [r["vr"] for r in rows]),
        ("Percent DC", "DC (%)", [r["pdc"] for r in rows]),
        ("Moment Magnitude", "Mw", [r["mw"] for r in rows]),
        ("Strike", "Strike (deg)", [r["plane1"]["strike"] for r in rows]),
        ("Rake", "Rake (deg)", [r["plane1"]["rake"] for r in rows]),
        ("Dip", "Dip (deg)", [r["plane1"]["dip"] for r in rows]),
    ]
    ev = solution["event"]
    pref = solution["preferred"]
    vr_max_depth = max(rows, key=lambda r: r["vr"])["depth_km"]
    dc_max_depth = max(rows, key=lambda r: r["pdc"])["depth_km"]

    # seaborn colorblind palette
    c_geonet, c_vr, c_dc = "#0173B2", "#DE8F05", "#029E73"
    ref_lines = [
        (ev["depth_km"], c_geonet,
         f"GeoNet depth {ev['depth_km']:g} km"),
        (vr_max_depth, c_vr, f"Max VR = {vr_max_depth:g} km"),
        (dc_max_depth, c_dc, f"Max DC = {dc_max_depth:g} km"),
    ]
    if pref["depth_km"] not in (vr_max_depth, dc_max_depth):
        ref_lines.append((pref["depth_km"], "black",
                          f"Preferred = {pref['depth_km']:g} km"))

    import config
    from obspy.imaging.beachball import beach

    stations = "_".join(r["station"] for r in solution["stations_used"])
    vr_max = max(r["vr"] for r in rows)
    window = [r for r in rows
              if r["vr"] >= vr_max - config.PREFER_DC_VR_TOLERANCE]
    if len(window) > 8:  # keep the strip readable
        idx = np.linspace(0, len(window) - 1, 8).round().astype(int)
        window = [window[i] for i in sorted(set(idx))]
    edge_for = {vr_max_depth: c_vr, dc_max_depth: c_dc,
                pref["depth_km"]: "black"}

    with plt.style.context("default"):
        fig, axes = plt.subplots(2, 3, figsize=(11, 9))
        for (title, ylabel, values), ax in zip(series, axes.ravel()):
            handles = [
                ax.axvline(d, color=c, lw=1.5, label=lab)
                for d, c, lab in ref_lines
            ]
            ax.plot(depths, values, "o", color="black", ms=4)
            ax.set_title(title, fontsize=11)
            ax.set_ylabel(ylabel, fontsize=9)
            ax.tick_params(labelsize=8)
        for ax in axes[1]:
            ax.set_xlabel("Depth (km)", fontsize=9)
        fig.suptitle(
            f"Depth Sensitivity - {ev['public_id']}  "
            f"{ev['origin_time'][:16]},  Stations: {stations}",
            fontsize=11,
        )
        fig.legend(handles=handles, loc="lower center",
                   ncol=len(ref_lines), fontsize=9, frameon=True,
                   bbox_to_anchor=(0.5, 0.185))
        fig.tight_layout(rect=[0, 0.24, 1, 1])

        # mechanism strip: solutions across the selection window (VR within
        # tolerance of max), where VR is flat but the mechanism can swing
        axb = fig.add_axes([0.06, 0.01, 0.88, 0.17])
        axb.set_xlim(0, max(len(window), 1))
        axb.set_ylim(0, 1.5)
        axb.set_aspect("equal")
        axb.axis("off")
        axb.set_anchor("S")
        for i, r in enumerate(window):
            tensor = r.get("tensor_rtp_dyne_cm")
            if tensor:
                fm = [tensor["MRR"], tensor["MTT"], tensor["MPP"],
                      tensor["MRT"], tensor["MRP"], tensor["MTP"]]
            else:  # older archives: DC mechanism only
                p1 = r["plane1"]
                fm = [p1["strike"], p1["dip"], p1["rake"]]
            edge = edge_for.get(r["depth_km"], "0.4")
            axb.add_collection(beach(
                fm, xy=(i + 0.5, 0.95), width=0.75, linewidth=1.4,
                facecolor="firebrick", edgecolor=edge))
            axb.text(i + 0.5, 0.38, f"{r['depth_km']:g} km", ha="center",
                     va="top", fontsize=8)
            axb.text(i + 0.5, 0.22, f"DC {r['pdc']:.0f}%", ha="center",
                     va="top", fontsize=8, color="0.35")
        axb.text(0.0, 1.45,
                 f"mechanism vs depth across the selection window "
                 f"(VR within {config.PREFER_DC_VR_TOLERANCE:g}% of max)",
                 fontsize=9, va="top")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
