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
    """Eight-panel figure: (a) station map with the deviatoric MT ball,
    (b-d) plane-1 E/N/U displacement, (e) solution-stability panel (DC
    ball + jackknife nodal-plane fan + uncertainty text), (f-h) plane-2
    E/N/U. Displacement AOI fits the modelled signal (minimum half-width
    12 km, maximum the modelled grid)."""
    import cartopy.crs as ccrs
    from obspy.imaging.beachball import beach

    ev = solution["event"]
    pref = solution["preferred"]
    lon0, lat0 = ev["longitude"], ev["latitude"]
    jk = solution.get("jackknife", {})
    roi = map_style.square_region(
        map_style.event_region(ev, solution["stations_used"], pad_deg=0.6))

    fig = plt.figure(figsize=(17.5, 11.2))
    row_y = {0: 0.565, 1: 0.11}
    panel_h = 0.365

    # ---- (a) station map: clean deviatoric mechanism ----------------------
    ax = map_style.geo_axes(fig, [0.02, row_y[0], 0.215, panel_h], roi)
    map_style.draw_context(ax, roi, ccrs, gnss=False)
    map_style.scale_bar(ax, roi, ccrs)
    import config as _config
    # dropped stations, marker-coded by cause, named so exclusions can be
    # audited from the figure alone
    grey = [d for d in solution["stations_dropped"] if "latitude" in d]

    def _cause(d):
        r = d.get("reason", "")
        if "amplitude outlier" in r:
            return "amp"
        if ("eliminated" in r or "test-drop" in r
                or "consistently bad" in r):
            return "elim"
        return "snr"

    styles = {
        "snr": dict(marker="^", color="0.65", markeredgecolor="0.4",
                    markeredgewidth=0.3, markersize=6),
        "amp": dict(marker="x", color="0.35", markeredgewidth=1.2,
                    markersize=6),
        "elim": dict(marker="^", color="white", markeredgecolor="0.25",
                     markeredgewidth=0.8, markersize=6),
    }
    for cause, st in styles.items():
        pts = [d for d in grey if _cause(d) == cause]
        if pts:
            ax.plot([d["longitude"] for d in pts],
                    [d["latitude"] for d in pts], linestyle="none",
                    transform=ccrs.PlateCarree(), zorder=5, **st)
    for d in grey:
        ax.annotate(
            d["station"].split(".")[1] if "." in d["station"]
            else d["station"],
            (d["longitude"], d["latitude"]), xytext=(3, -3),
            textcoords="offset points", fontsize=5.5, color="0.45",
            xycoords=ccrs.PlateCarree()._as_mpl_transform(ax))
    ax.text(0.985, 0.015,
            "dropped: grey=low SNR  x=bad amplitude  open=eliminated by fit",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.2,
            color="0.35", zorder=20,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8,
                      pad=1.5))
    used = solution["stations_used"]
    dv = [(r["zcor_s"] / (r["distance_km"] / _config.GROUP_VELOCITY_KMS))
          * 100.0 if "zcor_s" in r else None for r in used]
    if any(v is not None for v in dv):
        vals = [v for v in dv if v is not None]
        dvmax = max(5.0, max(abs(v) for v in vals))
        sc = ax.scatter(
            [r["longitude"] for r in used], [r["latitude"] for r in used],
            c=[v if v is not None else 0.0 for v in dv],
            cmap="coolwarm", vmin=-dvmax, vmax=dvmax, marker="^", s=90,
            edgecolors="black", linewidths=0.5,
            transform=ccrs.PlateCarree(), zorder=7,
        )
        caxa = fig.add_axes([0.045, 0.517, 0.15, 0.010])
        cba = fig.colorbar(sc, cax=caxa, orientation="horizontal")
        cba.set_label("dV% from zcor (red = model fast)", fontsize=7)
        cba.ax.tick_params(labelsize=6)
    else:
        ax.plot([r["longitude"] for r in used], [r["latitude"] for r in used],
                "^", color="forestgreen", markeredgecolor="black",
                markeredgewidth=0.4, markersize=8,
                transform=ccrs.PlateCarree(), zorder=7)
    for r in used:
        ax.annotate(
            r["station"], (r["longitude"], r["latitude"]),
            xytext=(4, -4), textcoords="offset points", fontsize=6.5,
            xycoords=ccrs.PlateCarree()._as_mpl_transform(ax),
        )
    map_style.add_beachball(ax, lon0, lat0, pref["tensor_rtp_dyne_cm"])
    map_style.panel_label(
        ax, f"(a) {ev['public_id']}  Mw {pref['mw']:.1f}  "
            f"depth {pref['depth_km']:g} km")

    # ---- displacement AOI: fit the modelled signal ------------------------
    plane_colors = {"plane1": "#0072B2", "plane2": "#E69F00"}
    vmax = max(
        0.1,
        max(float(np.abs(forward[p][c]).max())
            for p in ("plane1", "plane2")
            for c in ("ue_m", "un_m", "uz_m")) * 100.0,
    )
    thr = 0.02 * vmax
    x_km = forward["plane1"]["x_km"]
    y_km = forward["plane1"]["y_km"]
    sig = np.zeros((len(y_km), len(x_km)), dtype=bool)
    for p in ("plane1", "plane2"):
        for c in ("ue_m", "un_m", "uz_m"):
            sig |= np.abs(forward[p][c]) * 100.0 >= thr
    if sig.any():
        iy, ix = np.where(sig)
        half = max(12.0, 1.15 * max(
            abs(x_km[ix.min()]), abs(x_km[ix.max()]),
            abs(y_km[iy.min()]), abs(y_km[iy.max()])))
    else:
        half = 12.0
    half = min(half, float(x_km.max()))
    slon, slat = _local_km_to_geo(
        np.array([-half, half]), np.array([-half, half]), lon0, lat0)
    model_region = [float(slon[0]), float(slon[1]),
                    float(slat[0]), float(slat[1])]

    # ---- (b-d) plane 1 and (f-h) plane 2 ----------------------------------
    letters = {("plane1", 0): "b", ("plane1", 1): "c", ("plane1", 2): "d",
               ("plane2", 0): "f", ("plane2", 1): "g", ("plane2", 2): "h"}
    pm = None
    for row, plane in enumerate(("plane1", "plane2")):
        fw = forward[plane]
        fault = fw["fault"]
        outline = okada_forward.fault_outline(fault)
        olon, olat = _local_km_to_geo(
            np.array(outline["outline_x_km"]),
            np.array(outline["outline_y_km"]), lon0, lat0)
        tlon, tlat = _local_km_to_geo(
            np.array(outline["top_x_km"]), np.array(outline["top_y_km"]),
            lon0, lat0)
        glon, glat = _local_km_to_geo(fw["x_km"], fw["y_km"], lon0, lat0)
        for i, (comp, u_m) in enumerate(
                (("east", fw["ue_m"]), ("north", fw["un_m"]),
                 ("up", fw["uz_m"]))):
            u_cm = u_m * 100.0
            u_plot = np.ma.masked_where(np.abs(u_cm) < 0.02 * vmax, u_cm)
            axi = map_style.geo_axes(
                fig, [0.275 + 0.23 * i, row_y[row], 0.21, panel_h],
                model_region, labels=(i == 0),
            )
            n_gnss = map_style.draw_context(
                axi, model_region, ccrs, gnss_labels=(i == 0))
            if i == 0 and n_gnss == 0:
                axi.text(0.02, 0.03, "no operating GNSS marks in frame",
                         transform=axi.transAxes, fontsize=7.5,
                         color="#0173B2")
            map_style.scale_bar(axi, model_region, ccrs)
            pm = axi.pcolormesh(
                glon, glat, u_plot, cmap=map_style.vik(),
                vmin=-vmax, vmax=vmax, transform=ccrs.PlateCarree(),
                alpha=0.9, shading="auto", zorder=3,
            )
            # USGS finite-fault convention: plane outlined in black,
            # up-dip (shallowest) edge as a distinct red line
            axi.plot(olon, olat, "-", color="black",
                     linewidth=1.0, transform=ccrs.PlateCarree(), zorder=10)
            axi.plot(tlon, tlat, "-", color=plane_colors[plane],
                     linewidth=2.8, transform=ccrs.PlateCarree(), zorder=11)
            if i == 0:  # dip-direction arrow on the first panel of the row
                dlon, dlat = _local_km_to_geo(
                    np.array(outline["dip_arrow_x_km"]),
                    np.array(outline["dip_arrow_y_km"]), lon0, lat0)
                axi.annotate(
                    "", xy=(dlon[1], dlat[1]), xytext=(dlon[0], dlat[0]),
                    xycoords=ccrs.PlateCarree()._as_mpl_transform(axi),
                    textcoords=ccrs.PlateCarree()._as_mpl_transform(axi),
                    arrowprops=dict(arrowstyle="-|>", linewidth=1.4,
                                    color="black"), zorder=12)
                axi.annotate(
                    "dip", (dlon[1], dlat[1]), xytext=(3, 3),
                    textcoords="offset points", fontsize=7,
                    color="black",
                    xycoords=ccrs.PlateCarree()._as_mpl_transform(axi))
            axi.plot(lon0, lat0, "*", color="yellow",
                     markeredgecolor="black", markersize=9,
                     transform=ccrs.PlateCarree(), zorder=11)
            map_style.panel_label(
                axi, f"({letters[plane, i]}) plane {row + 1} {comp}  "
                     f"(peak {np.abs(u_cm).max():.2f} cm)")

    # ---- (e) solution stability: DC ball + jackknife fan + text -----------
    axe = fig.add_axes([0.02, row_y[1], 0.215, panel_h])
    axe.set_xlim(-1.35, 1.35)
    axe.set_ylim(-2.3, 1.35)
    axe.set_aspect("equal")
    axe.axis("off")
    axe.add_collection(beach(
        [pref["plane1"]["strike"], pref["plane1"]["dip"],
         pref["plane1"]["rake"]], xy=(0, 0.05), width=2.0, linewidth=0.8,
        facecolor="firebrick"))
    for sub in jk.get("subsets", []):
        for key in ("plane1", "plane2"):
            p = sub.get(key)
            if p:
                xa, ya = map_style.nodal_plane_arc(p["strike"], p["dip"])
                axe.plot(xa, ya + 0.05, "-", color="0.1", linewidth=1.0,
                         alpha=0.5, zorder=110)
    for key, colr in (("plane1", plane_colors["plane1"]),
                      ("plane2", plane_colors["plane2"])):
        xa, ya = map_style.nodal_plane_arc(
            pref[key]["strike"], pref[key]["dip"])
        axe.plot(xa, ya + 0.05, "-", color=colr, linewidth=2.2, zorder=120,
                 solid_capstyle="round")
    if jk.get("n_subsets"):
        text = (
            f"leave-one-out jackknife (n={jk['n_subsets']}):\n"
            f"Mw {pref['mw']:.2f} $\\pm$ {jk['mw_std']}   "
            f"DC {pref['pdc']:.0f} $\\pm$ {jk['dc_std']}%\n"
            f"mechanism rotation $\\leq$ "
            f"{jk['max_tensor_rotation_deg']:g}$^\\circ$\n"
            f"grey fan = subset nodal planes"
        )
    else:
        text = "jackknife skipped (fewer than 4 stations)"
    axe.text(0, -1.45, text, ha="center", va="top", fontsize=8.5)
    axe.text(0.02, 0.985, "(e) solution stability",
             transform=axe.transAxes, va="top", ha="left", fontsize=9,
             bbox=dict(facecolor="white", edgecolor="0.4", linewidth=0.5,
                       boxstyle="square,pad=0.25"))

    cax = fig.add_axes([0.955, 0.30, 0.009, 0.45])
    cb = fig.colorbar(pm, cax=cax, orientation="vertical")
    cb.set_label("predicted displacement [cm]", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    peak_cm = forward["peak_abs_m"] * 100.0
    detect = ("potentially InSAR detectable" if forward["detectable"]
              else "below InSAR detection")

    # ---- footer -----------------------------------------------------------
    prov = solution["provenance"]
    p1f = forward["plane1"]["fault"]
    p2f = forward["plane2"]["fault"]
    lines = [
        f"Okada forward models (black outline; COLOURED line = up-dip/shallowest edge, matching that plane's arc on the stability ball; arrow points down-dip): "
        f"plane 1 (blue, top row) {p1f['strike']:.0f}/{p1f['dip']:.0f}/"
        f"{p1f['rake']:.0f}, plane 2 (orange, bottom row) {p2f['strike']:.0f}/"
        f"{p2f['dip']:.0f}/{p2f['rake']:.0f}; "
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
        fig.text(0.02, 0.082 - 0.030 * i, line, fontsize=8.5,
                 family="monospace", va="top")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def make_overview_map(events_dir: Path, out_path: Path) -> Path:
    """Standalone NZ map of every archived solution: full-MT beachballs
    (size = Mw, solid = grade A/B, washed = C/D) over the NZ Active Faults
    Database, with an in-figure legend. No lat/lon grid clutter."""
    import json

    import cartopy.crs as ccrs
    from obspy.imaging.beachball import beach

    sols = []
    for p in sorted(events_dir.glob("*/solution.json")):
        try:
            sols.append(json.loads(p.read_text()))
        except Exception:  # noqa: BLE001
            continue
    region = [163.5, 183.0, -50.7, -33.3]
    fig = plt.figure(figsize=(7.5, 8.7))
    ax = map_style.geo_axes(fig, [0.07, 0.05, 0.9, 0.88], region,
                            grid=False, label_size=11)
    map_style.draw_context(ax, region, ccrs, gnss=False)
    dates = []
    for sol in sols:
        ev = sol["event"]
        pref = sol["preferred"]
        lon = ev["longitude"] % 360.0
        grade = sol["quality"].get("grade", "?")
        width = 0.018 + 0.010 * max(0.0, pref["mw"] - 4.0)
        try:
            rtp = pref["tensor_rtp_dyne_cm"]
        except KeyError:
            continue
        x, y = ax.projection.transform_point(lon, ev["latitude"],
                                             ccrs.PlateCarree())
        xe0, xe1, _, _ = ax.get_extent(crs=ax.projection)
        if grade in ("A", "B"):
            # well-constrained: show the mechanism (filled beachball)
            fm = [rtp["MRR"], rtp["MTT"], rtp["MPP"],
                  rtp["MRT"], rtp["MRP"], rtp["MTP"]]
            ax.add_collection(beach(
                fm, xy=(x, y), width=width * (xe1 - xe0),
                linewidth=0.4, facecolor="black"))
        else:
            # poorly constrained: open circle — location and size only,
            # no mechanism the data cannot support
            ax.add_patch(plt.Circle(
                (x, y), radius=0.5 * width * (xe1 - xe0),
                facecolor="none", edgecolor="black", linewidth=0.9,
                zorder=8))
        dates.append(ev["origin_time"][:10])
    ax.set_title(
        f"Catalogue Solutions\n"
        f"{len(sols)} events, {min(dates)} to {max(dates)}",
        fontsize=13,
    )

    # in-figure legend, bottom-right (empty ocean), strike-slip demo balls
    xe0, xe1, ye0, ye1 = ax.get_extent(crs=ax.projection)
    w_ax, h_ax = xe1 - xe0, ye1 - ye0
    lx, ly = xe0 + 0.595 * w_ax, ye0 + 0.045 * h_ax
    import matplotlib.patches as mpatches
    ax.add_patch(mpatches.Rectangle(
        (lx - 0.015 * w_ax, ly - 0.030 * h_ax), 0.40 * w_ax,
        0.27 * h_ax, facecolor="white", edgecolor="0.4",
        linewidth=0.6, zorder=200))
    ss = [0.0, 90.0, 0.0]  # clean strike-slip DC for the demo balls
    for i, (mw, lab) in enumerate([(4.0, "Mw 4"), (5.0, "Mw 5"),
                                   (6.0, "Mw 6")]):
        bx = lx + (0.055 + 0.125 * i) * w_ax
        by = ly + 0.185 * h_ax
        bw = (0.018 + 0.010 * (mw - 4.0)) * w_ax
        ax.add_collection(beach(ss, xy=(bx, by), width=bw,
                                linewidth=0.4, facecolor="black",
                                zorder=210))
        ax.text(bx, by - 0.055 * h_ax, lab, ha="center",
                fontsize=10, zorder=220)
    b1 = beach(ss, xy=(lx + 0.045 * w_ax, ly + 0.065 * h_ax),
               width=0.024 * w_ax, linewidth=0.4, facecolor="black",
               zorder=210)
    ax.add_collection(b1)
    ax.text(lx + 0.075 * w_ax, ly + 0.065 * h_ax, "grade A/B",
            fontsize=10, va="center", zorder=220)
    ax.add_patch(plt.Circle(
        (lx + 0.22 * w_ax, ly + 0.065 * h_ax), radius=0.012 * w_ax,
        facecolor="none", edgecolor="black", linewidth=0.9, zorder=210))
    ax.text(lx + 0.25 * w_ax, ly + 0.065 * h_ax,
            "grade C/D (poorly constrained)", fontsize=10, va="center",
            zorder=220)
    ax.plot([lx + 0.03 * w_ax, lx + 0.075 * w_ax],
            [ly + 0.005 * h_ax] * 2, "-", color="#8B3A3A", linewidth=1.4,
            zorder=210)
    ax.text(lx + 0.09 * w_ax, ly + 0.005 * h_ax,
            "active faults (NZAFD, GNS)", fontsize=10, va="center",
            zorder=220)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
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
    c_geonet, c_vr, c_dc = "#0072B2", "#E69F00", "#CC79A7"
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
    # always include the GLOBAL max-DC depth (the teal reference line must
    # have its ball in the strip), plus GeoNet, max VR and the preferred
    key_depths = {geonet_depth, vr_max_depth, dc_max_depth, dc_win_depth,
                  pref["depth_km"]}
    i_pref = depths.index(pref["depth_km"])
    step = 1
    limit = min(8, len(depths))
    while len(key_depths) < limit:
        for j in (i_pref - step, i_pref + step):
            if 0 <= j < len(depths) and len(key_depths) < limit:
                key_depths.add(depths[j])
        step += 1
    window = [r for r in rows if r["depth_km"] in key_depths]
    window.sort(key=lambda r: r["depth_km"])
    # rim colours match the reference lines; preferred wins a collision
    edge_for = {geonet_depth: c_geonet, vr_max_depth: c_vr,
                dc_max_depth: c_dc, dc_win_depth: c_dc,
                pref["depth_km"]: "black"}

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
            is_pref = r["depth_km"] == pref["depth_km"]
            ball = beach(
                fm, xy=(i + 0.5, 0.95), width=0.75,
                linewidth=1.8 if is_pref else 1.2,
                facecolor="firebrick", edgecolor=edge)
            if not is_pref:  # wash out everything but the selected solution
                ball.set_alpha(0.45)
            axb.add_collection(ball)
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
