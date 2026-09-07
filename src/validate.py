"""Compare archived automated solutions against (1) the published NZ
regional CMT solutions (Ristau catalogue, GeoNet/data repository) and
(2) the Global CMT catalogue (Ekstrom et al., globalcmt.org; monthly +
quick NDK feeds), with GeoNet as the source of every original hypocentre.

    pixi run python src/validate.py

Catalogue-level, not per event: it reads every archived solution.json.
Run by hand after a sweep, and automatically by sync_repo.py so the
validation report is refreshed every time the archive is published.

Writes validation/comparison.csv and per-metric figures plus a
terminal summary. Metrics per common event: dMw, dDepth, and the angle
between the double-couple tensors built from each catalogue's nodal plane
(basis-independent; 0 = identical mechanism).
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config
from geonet import load_geonet_cmt
from invert import min_rotation_angle_deg

USGS_URL = ("https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson"
            "&starttime=2026-01-01&minmagnitude=4.5&minlatitude=-50.5"
            "&maxlatitude=-33&minlongitude=164&maxlongitude=182.5"
            "&producttype=moment-tensor")


def load_usgs():
    """USGS/NEIC events with moment-tensor products for the NZ box.
    Returns list of (UTCDateTime, lat, lon, detail_url); [] on failure."""
    import requests
    from obspy import UTCDateTime

    try:
        r = requests.get(USGS_URL, timeout=120,
                         headers={"User-Agent": config.USER_AGENT})
        r.raise_for_status()
        out = []
        for f in r.json()["features"]:
            lon, lat = f["geometry"]["coordinates"][:2]
            out.append((UTCDateTime(f["properties"]["time"] / 1000.0),
                        lat, lon, f["properties"]["detail"]))
        return out
    except Exception as e:  # noqa: BLE001
        print(f"USGS query failed: {e}")
        return []


def match_usgs(usgs, origin_time, lat, lon):
    """USGS moment-tensor properties for the matching event, else None."""
    import requests
    from obspy import UTCDateTime

    t0 = UTCDateTime(origin_time)
    for t, la, lo, detail in usgs:
        if abs(t - t0) < 90 and abs(la - lat) < 1.5 \
                and abs((lo - lon + 180) % 360 - 180) < 1.5:
            try:
                d = requests.get(detail, timeout=60).json()
                mt = d["properties"]["products"]["moment-tensor"][0]
                p = mt["properties"]
                return {
                    "Mw": float(p["derived-magnitude"]),
                    "depth": float(p.get("derived-depth",
                                         p.get("depth"))),
                    "strike": float(p["nodal-plane-1-strike"]),
                    "dip": float(p["nodal-plane-1-dip"]),
                    "rake": float(p["nodal-plane-1-rake"]),
                }
            except Exception as e:  # noqa: BLE001
                print(f"USGS detail parse failed: {e}")
                return None
    return None


GCMT_BASE = "https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog"


def load_gcmt_2026():
    """Global CMT solutions for 2026 from monthly + quick NDK feeds.
    Returns an obspy Catalog (possibly empty on network failure)."""
    from obspy import read_events

    months = ["jan", "feb", "mar", "apr", "may", "jun",
              "jul", "aug", "sep", "oct", "nov", "dec"]
    cat = None
    for m in months:
        url = f"{GCMT_BASE}/NEW_MONTHLY/2026/{m}26.ndk"
        try:
            c = read_events(url)
        except Exception:
            continue
        cat = c if cat is None else cat + c
    try:
        q = read_events(f"{GCMT_BASE}/NEW_QUICK/qcmt.ndk")
        cat = q if cat is None else cat + q
    except Exception:
        pass
    return cat


def match_gcmt(cat, origin_time, lat, lon):
    """Nearest GCMT event within 90 s and 1.5 deg, else None."""
    from obspy import UTCDateTime

    if cat is None:
        return None
    t0 = UTCDateTime(origin_time)
    best, best_dt = None, 90.0
    for ev in cat:
        o = ev.preferred_origin() or ev.origins[0]
        dt = abs(o.time - t0)
        if dt < best_dt and abs(o.latitude - lat) < 1.5 \
                and abs((o.longitude - lon + 180) % 360 - 180) < 1.5:
            best, best_dt = ev, dt
    return best


def main() -> None:
    usgs = load_usgs()
    print(f"USGS/NEIC: {len(usgs)} candidate events with MT products")
    gcmt = load_gcmt_2026()
    print(f"Global CMT: {len(gcmt) if gcmt else 0} events loaded")
    cmt = load_geonet_cmt()
    id_col = "PublicID" if "PublicID" in cmt.columns else "EVENT_ID"
    cmt = cmt.set_index(cmt[id_col].astype(str))
    print(f"reference catalogue: {len(cmt)} solutions "
          f"(columns: {list(cmt.columns)[:12]}...)")

    rows = []
    origins = {}
    for path in config.solution_paths():
        sol = json.loads(path.read_text())
        if not config.is_solved(sol):
            continue  # no mechanism to compare
        ev = sol["event"]
        pid = ev["public_id"]
        origins[pid] = (ev["origin_time"], ev["latitude"], ev["longitude"])
        pref = sol["preferred"]
        p1 = pref["plane1"]
        base = {
            "PublicID": pid,
            "GeoNet_M": round(ev["prelim_mag"], 2),
            "GeoNet_depth": ev["depth_km"],
            "our_Mw": round(pref["mw"], 2),
            "our_depth": pref["depth_km"],
            "our_VR": round(pref["vr"], 1),
            "our_DC": round(pref["pdc"], 0),
            "grade": sol["quality"].get("grade", "?"),
            "n_stations": sol["quality"]["n_stations_used"],
        }
        if pid in cmt.index:
            ref = cmt.loc[pid]
            rows.append({
                **base, "reference": "NZ_CMT_Ristau",
                "ref_Mw": float(ref["Mw"]),
                "dMw": round(pref["mw"] - float(ref["Mw"]), 2),
                "ref_depth": float(ref["CD"]),
                "dDepth": round(pref["depth_km"] - float(ref["CD"]), 1),
                "rotation_angle_deg": round(min_rotation_angle_deg(
                    (p1["strike"], p1["dip"], p1["rake"]),
                    (float(ref["strike1"]), float(ref["dip1"]),
                     float(ref["rake1"]))), 1),
            })
        u = match_usgs(usgs, ev["origin_time"], ev["latitude"],
                       ev["longitude"])
        if u is not None:
            rows.append({
                **base, "reference": "USGS_NEIC",
                "ref_Mw": round(u["Mw"], 2),
                "dMw": round(pref["mw"] - u["Mw"], 2),
                "ref_depth": round(u["depth"], 1),
                "dDepth": round(pref["depth_km"] - u["depth"], 1),
                "rotation_angle_deg": round(min_rotation_angle_deg(
                    (p1["strike"], p1["dip"], p1["rake"]),
                    (u["strike"], u["dip"], u["rake"])), 1),
            })
        g = match_gcmt(gcmt, ev["origin_time"], ev["latitude"],
                       ev["longitude"])
        if g is not None:
            fm = g.preferred_focal_mechanism() or g.focal_mechanisms[0]
            np1 = fm.nodal_planes.nodal_plane_1
            mag = (g.preferred_magnitude() or g.magnitudes[0]).mag
            dep = None
            for o in g.origins:
                if o.origin_type and "centroid" in str(o.origin_type):
                    dep = o.depth / 1000.0
            if dep is None:
                dep = (g.preferred_origin() or g.origins[0]).depth / 1000.0
            rows.append({
                **base, "reference": "GlobalCMT",
                "ref_Mw": round(float(mag), 2),
                "dMw": round(pref["mw"] - float(mag), 2),
                "ref_depth": round(dep, 1),
                "dDepth": round(pref["depth_km"] - dep, 1),
                "rotation_angle_deg": round(min_rotation_angle_deg(
                    (p1["strike"], p1["dip"], p1["rake"]),
                    (np1.strike, np1.dip, np1.rake)), 1),
            })

    assert rows, ("no common events with the reference catalogue yet — "
                  "rerun after it next updates")
    df = pd.DataFrame(rows).sort_values("PublicID")
    out_dir = config.REPO_DIR / "validation"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "comparison.csv", index=False)

    print(df.to_string(index=False))
    for ref, sub in df.groupby("reference"):
        print(f"\n[{ref}] n = {len(sub)}")
        print(f"  Mw:    mean dMw {sub.dMw.mean():+.2f}, |dMw| median "
              f"{sub.dMw.abs().median():.2f}")
        print(f"  depth: mean dZ {sub.dDepth.mean():+.1f} km, |dZ| median "
              f"{sub.dDepth.abs().median():.1f} km")
        print(f"  mechanism: median min rotation "
              f"{sub.rotation_angle_deg.median():.0f} deg")

    # inter-reference baseline: how far apart the two reference catalogues
    # are from each other on co-matched events (Ristau vs USGS), for
    # context; drawn on the rotation figure
    from invert import min_rotation_angle_deg as _rot
    gsdr = {str(x["PublicID"]): (float(x["strike1"]), float(x["dip1"]),
                                 float(x["rake1"]))
            for _, x in cmt.iterrows()}
    ref_ref = []
    for pid, sub in df.groupby("PublicID"):
        if "USGS_NEIC" not in set(sub.reference) or pid not in gsdr:
            continue
        row = origins.get(pid)
        if row is None:
            continue
        m = match_usgs(usgs, *row)
        if m and "strike" in m:
            ref_ref.append(_rot(gsdr[pid],
                                (m["strike"], m["dip"], m["rake"])))
    baseline = float(np.median(ref_ref)) if ref_ref else None
    if baseline is not None:
        print(f"\ninter-reference baseline (Ristau vs USGS, n={len(ref_ref)}):"
              f" median min rotation {baseline:.0f} deg between the two"
              " reference catalogues")

    colors = {"NZ_CMT_Ristau": "#0072B2", "GlobalCMT": "#E69F00",
              "USGS_NEIC": "#CC79A7"}

    # Density-shaded scatter, as in the author's InSAR validation figures:
    # the pairs are binned on a fine grid and every occupied cell is a small
    # square shaded by its normalised count on a reversed grey ramp (light
    # grey -> black, saturating at 0.3 so the dense core is black), dense
    # cells drawn last; a 1:1 line; RMSE, R2 and slope in the corner. The
    # two small independent references are overlaid as coloured marks.
    # Every panel is split A/B vs C/D because the grade decides publication.
    from scipy import stats as _stats
    rist = df[df.reference == "NZ_CMT_Ristau"]
    tiers = [("grade A/B (published tier)", rist[rist.grade.isin(["A", "B"])],
              df[df.grade.isin(["A", "B"])]),
             ("grade C/D (archive only)", rist[rist.grade.isin(["C", "D"])],
              df[df.grade.isin(["C", "D"])])]

    def _density(ax, x, y, bin_size, xext, yext=None):
        """Scatter shaded by a Gaussian kernel-density estimate of each
        point's neighbourhood (devon ramp: pale = isolated, dark = the dense
        core), densest points drawn last. Reads at a few hundred events
        where a binned heat map is still sparse; bin_size only sets the
        kernel scale relative to the axis range."""
        x = np.asarray(x, float); y = np.asarray(y, float)
        if len(x) < 4:
            ax.scatter(x, y, s=18, c="#4a3f8a", linewidths=0, zorder=3)
            return
        yext = yext or xext
        xs = (x - xext[0]) / (xext[1] - xext[0])    # equalise the axes
        ys = (y - yext[0]) / (yext[1] - yext[0])
        z = _stats.gaussian_kde(np.vstack([xs, ys]), bw_method=0.12)(
            np.vstack([xs, ys]))
        z = z / z.max()
        order = np.argsort(z)
        from cmcrameri import cm as _cmc
        ax.scatter(x[order], y[order], c=z[order], cmap=_cmc.devon_r,
                   vmin=-0.15, vmax=1.0, s=18, linewidths=0, zorder=3)

    def _fit_text(ax, x, y, unit):
        x = np.asarray(x, float); y = np.asarray(y, float)
        if len(x) < 3:
            return
        rmse = float(np.sqrt(np.sum((x - y) ** 2) / (len(x) - 1)))
        lr = _stats.linregress(x, y)
        ax.text(0.97, 0.10, f"RMSE {rmse:.2f} {unit}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=9)
        ax.text(0.97, 0.03, f"R\u00b2 {lr.rvalue ** 2:.2f}, slope {lr.slope:.2f}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=9)

    def _overlay(ax, sub, xcol, ycol):
        for ref in ("USGS_NEIC", "GlobalCMT"):
            s = sub[sub.reference == ref]
            if len(s):
                ax.scatter(s[xcol], s[ycol], s=30, c=colors[ref],
                           edgecolors="white", linewidths=1.0,
                           label=ref, zorder=5)

    # ---- figure 1: magnitude -------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 8.6))
    ext_mw = (3.4, 6.6)
    for row, (title, r, sub) in enumerate(tiers):
        ax = axes[row, 0]
        _density(ax, r.ref_Mw, r.our_Mw, 0.05, ext_mw)
        ax.plot(ext_mw, ext_mw, "-", color="black", linewidth=1, zorder=2)
        _overlay(ax, sub, "ref_Mw", "our_Mw")
        _fit_text(ax, r.ref_Mw, r.our_Mw, "")
        ax.set_xlim(ext_mw); ax.set_ylim(ext_mw); ax.set_aspect("equal")
        ax.set_xlabel("published Mw"); ax.set_ylabel("auto Mw")
        ax.set_title(f"{title}, n = {len(r)} vs Ristau", fontsize=10)
        ax = axes[row, 1]
        _density(ax, r.ref_Mw, r.dMw, 0.05, ext_mw, (-1.0, 1.0))
        ax.axhline(0, color="black", linewidth=1, zorder=2)
        _overlay(ax, sub, "ref_Mw", "dMw")
        ax.set_xlim(ext_mw); ax.set_ylim(-1.0, 1.0)
        ax.set_xlabel("published Mw"); ax.set_ylabel("auto - published Mw")
        ax.set_title(f"residual: median |dMw| {r.dMw.abs().median():.2f}",
                     fontsize=10)
    axes[0, 0].legend(fontsize=8, loc="upper left")
    fig.suptitle("moment magnitude vs published catalogues "
                 "(shading: Ristau NZ CMT density; marks: USGS, GCMT)")
    fig.tight_layout()
    fig.savefig(out_dir / "comparison_mw.jpg", dpi=150)
    plt.close(fig)

    # ---- figure 2: depth ------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.9))
    lim = float(np.ceil(max(df.ref_depth.max(), df.our_depth.max()) / 10) * 10)
    for ax, (title, r, sub) in zip(axes, tiers):
        _density(ax, r.ref_depth, r.our_depth, 1.0, (0, lim))
        ax.plot([0, lim], [0, lim], "-", color="black", linewidth=1, zorder=2)
        _overlay(ax, sub, "ref_depth", "our_depth")
        _fit_text(ax, r.ref_depth, r.our_depth, "km")
        ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal")
        ax.set_xlabel("published depth (km)"); ax.set_ylabel("auto depth (km)")
        ax.set_title(f"{title}, n = {len(r)} vs Ristau: "
                     f"median |dZ| {r.dDepth.abs().median():.0f} km",
                     fontsize=10)
    axes[0].legend(fontsize=8, loc="upper left")
    fig.suptitle("centroid depth vs published catalogues "
                 "(shading: Ristau NZ CMT density; marks: USGS, GCMT)")
    fig.tight_layout()
    fig.savefig(out_dir / "comparison_depth.jpg", dpi=150)
    plt.close(fig)

    # ---- figure 3: mechanism rotation, split by grade -------------------
    # split by grade because C/D dominates the count; the inter-reference
    # baseline gives the context for the A/B tier
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    rist = df[df.reference == "NZ_CMT_Ristau"]
    bins = np.arange(0, 121, 10)
    ax.hist(rist[rist.grade.isin(["A", "B"])].rotation_angle_deg, bins=bins,
            color="#0072B2", alpha=0.85, label="grade A/B (published tier)")
    ax.hist(rist[rist.grade.isin(["C", "D"])].rotation_angle_deg, bins=bins,
            color="0.65", alpha=0.7, label="grade C/D (archive only)")
    if baseline is not None:
        ax.axvline(baseline, color="#D55E00", linestyle="--", linewidth=1.6,
                   label=f"Ristau vs USGS ({baseline:.0f}\N{DEGREE SIGN}: "
                         "difference between the two reference catalogues)")
    ax.set_xlabel("minimum rotation angle vs Ristau NZ CMT (deg)")
    ax.set_ylabel("events")
    ax.legend(fontsize=8)
    fig.suptitle("mechanism agreement, split by quality grade")
    fig.tight_layout()
    fig.savefig(out_dir / "comparison_rotation.jpg", dpi=150)
    plt.close(fig)

    print(f"\nwrote {out_dir}/comparison.csv and comparison_mw/"
          "comparison_depth/comparison_rotation .jpg")


if __name__ == "__main__":
    main()
