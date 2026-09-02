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

    fig = plt.figure(figsize=(17.5, 6.4))

    # ---- panel a: stations + full-MT beachball over the station ROI -------
    ax = map_style.geo_axes(fig, [0.03, 0.30, 0.19, 0.60], roi)
    map_style.draw_context(ax, roi, ccrs, gnss=False)
    map_style.scale_bar(ax, roi, ccrs)
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
            fig, [0.245 + 0.23 * i, 0.30, 0.22, 0.60], model_region,
            labels=(i == 0),
        )
        pm = axi.pcolormesh(
            glon, glat, u_plot, cmap=map_style.vik(), vmin=-vmax, vmax=vmax,
            transform=ccrs.PlateCarree(), alpha=0.9, shading="auto", zorder=3,
        )
        n_gnss = map_style.draw_context(axi, model_region, ccrs,
                                        gnss_labels=(i == 0))
        if i == 0 and n_gnss == 0:
            axi.text(0.02, 0.03, "no operating GNSS marks in frame",
                     transform=axi.transAxes, fontsize=7.5, color="#0173B2")
        map_style.scale_bar(axi, model_region, ccrs)
        # surface projection of the modelled plane (bold edge = up-dip)
        axi.plot(olon, olat, "--", color="black", linewidth=1.2,
                 transform=ccrs.PlateCarree(), zorder=10)
        axi.plot(tlon, tlat, "-", color="black", linewidth=2.6,
                 transform=ccrs.PlateCarree(), zorder=10)
        axi.plot(lon0, lat0, "*", color="yellow", markeredgecolor="black",
                 markersize=9, transform=ccrs.PlateCarree(), zorder=11)
        map_style.panel_label(
            axi, f"{name}  (peak {np.abs(u_cm).max():.2f} cm)")

    # which nodal plane the forward model used, highlighted on a DC ball
    map_style.inset_dc_ball(
        fig, [0.875, 0.015, 0.075, 0.24], pref["plane1"], pref["plane2"])

    cax = fig.add_axes([0.945, 0.34, 0.009, 0.52])
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
    band_hz = solution.get("filter_band_hz")
    band_txt = (f"{1/band_hz[1]:.0f}-{1/band_hz[0]:.0f} s"
                if band_hz else "n/a")
    lines.append(
        f"PRELIMINARY automated solution - model {prov['velocity_model']}, "
        f"filter {band_txt}, mttime {prov['mttime_version']}, "
        f"VR {pref['vr']:.0f}%, DC {pref['pdc']:.0f}% / "
        f"CLVD {pref['pclvd']:.0f}%"
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
    in_window = [r for r in rows
                 if r["vr"] >= vr_max - config.PREFER_DC_VR_TOLERANCE]
    dc_win_depth = max(in_window, key=lambda r: r["pdc"])["depth_km"]
    geonet_depth = min(depths, key=lambda d: abs(d - ev["depth_km"]))

    # strip: ~7 mechanisms pinned to the key depths (GeoNet, max VR,
    # windowed max DC, preferred) with neighbours of the preferred depth
    # filled in, to make the depth trade-off around the solution visual
    key_depths = {geonet_depth, vr_max_depth, dc_win_depth,
                  pref["depth_km"]}
    i_pref = depths.index(pref["depth_km"])
    step = 1
    while len(key_depths) < min(7, len(depths)):
        for j in (i_pref - step, i_pref + step):
            if 0 <= j < len(depths) and len(key_depths) < 7:
                key_depths.add(depths[j])
        step += 1
    window = [r for r in rows if r["depth_km"] in key_depths]
    window.sort(key=lambda r: r["depth_km"])
    # rim colours match the reference lines; preferred wins a collision
    edge_for = {geonet_depth: c_geonet, vr_max_depth: c_vr,
                dc_win_depth: c_dc, pref["depth_km"]: "black"}

    with plt.style.context("default"):
        fig, axes = plt.subplots(1, 3, figsize=(11, 6))
        for (title, ylabel, values), ax in zip(series, axes.ravel()):
            handles = [
                ax.axvline(d, color=c, lw=1.5, label=lab)
                for d, c, lab in ref_lines
            ]
            ax.plot(depths, values, "o", color="black", ms=4)
            ax.set_title(title, fontsize=11)
            ax.set_ylabel(ylabel, fontsize=9)
            ax.tick_params(labelsize=8)
        for ax in axes:
            ax.set_xlabel("Depth (km)", fontsize=9)
        band_hz = solution.get("filter_band_hz")
        band_txt = (f",  Filter: {1/band_hz[1]:.0f}-{1/band_hz[0]:.0f} s"
                    if band_hz else "")
        fig.suptitle(
            f"Depth Sensitivity - {ev['public_id']}  "
            f"{ev['origin_time'][:16]}{band_txt},  Stations: {stations}",
            fontsize=11,
        )
        fig.legend(handles=handles, loc="lower center",
                   ncol=len(ref_lines), fontsize=9, frameon=True,
                   bbox_to_anchor=(0.5, 0.30))
        fig.tight_layout(rect=[0, 0.38, 1, 1])

        # mechanism strip: solutions across the selection window (VR within
        # tolerance of max), where VR is flat but the mechanism can swing
        axb = fig.add_axes([0.06, 0.02, 0.88, 0.26])
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
                 "deviatoric mechanism vs depth around the solution "
                 "(rims: blue=GeoNet, orange=max VR, teal=max DC, "
                 "black=preferred)",
                 fontsize=9, va="top")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
