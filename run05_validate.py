"""Compare archived automated solutions against (1) the published NZ
regional CMT solutions (Ristau catalogue, GeoNet/data repository) and
(2) the Global CMT catalogue (Ekstrom et al., globalcmt.org; monthly +
quick NDK feeds), with GeoNet as the source of every original hypocentre.

    pixi run python run05_validate.py

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
    for path in sorted(config.EVENTS_DIR.glob("*/solution.json")):
        sol = json.loads(path.read_text())
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

    # inter-reference baseline: how far apart the REFERENCES are from
    # each other on co-matched events (Ristau vs USGS) — the floor any
    # catalogue could reach; drawn on the rotation figure
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
              f" median min rotation {baseline:.0f} deg — the floor any"
              " catalogue could reach")

    colors = {"NZ_CMT_Ristau": "#0072B2", "GlobalCMT": "#E69F00",
              "USGS_NEIC": "#CC79A7"}

    # ---- figure 1: magnitude -------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4))
    for ref, sub in df.groupby("reference"):
        c = colors.get(ref, "black")
        axes[0].scatter(sub.ref_Mw, sub.our_Mw, c=c, s=25, label=ref)
        axes[1].scatter(sub.ref_Mw, sub.dMw, c=c, s=25)
    axes[0].plot([3.5, 6.5], [3.5, 6.5], "-", color="0.7", zorder=0)
    axes[0].set_xlabel("published Mw")
    axes[0].set_ylabel("auto Mw")
    axes[0].legend(fontsize=8)
    axes[1].axhline(0, color="0.7")
    axes[1].set_xlabel("published Mw")
    axes[1].set_ylabel("auto - published Mw")
    fig.suptitle("moment magnitude vs published catalogues")
    fig.tight_layout()
    fig.savefig(out_dir / "comparison_mw.jpg", dpi=150)
    plt.close(fig)

    # ---- figure 2: depth ------------------------------------------------
    fig, ax = plt.subplots(figsize=(5.4, 5))
    for ref, sub in df.groupby("reference"):
        ax.scatter(sub.ref_depth, sub.our_depth,
                   c=colors.get(ref, "black"), s=25, label=ref)
    lim = max(df.ref_depth.max(), df.our_depth.max()) + 5
    ax.plot([0, lim], [0, lim], "-", color="0.7", zorder=0)
    ax.set_xlabel("published depth (km)")
    ax.set_ylabel("auto depth (km)")
    ax.legend(fontsize=8)
    fig.suptitle("centroid depth vs published catalogues")
    fig.tight_layout()
    fig.savefig(out_dir / "comparison_depth.jpg", dpi=150)
    plt.close(fig)

    # ---- figure 3: mechanism rotation, split by grade -------------------
    # the all-grades histogram alone is misleading (C/D dominates the
    # count); the grade split + the inter-agency floor is the story
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    rist = df[df.reference == "NZ_CMT_Ristau"]
    bins = np.arange(0, 121, 10)
    ax.hist(rist[rist.grade.isin(["A", "B"])].rotation_angle_deg, bins=bins,
            color="#0072B2", alpha=0.85, label="grade A/B (published tier)")
    ax.hist(rist[rist.grade.isin(["C", "D"])].rotation_angle_deg, bins=bins,
            color="0.65", alpha=0.7, label="grade C/D (archive only)")
    if baseline is not None:
        ax.axvline(baseline, color="#D55E00", linestyle="--", linewidth=1.6,
                   label=f"Ristau-vs-USGS floor ({baseline:.0f}\N{DEGREE SIGN}: "
                         "how far apart the references themselves are)")
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
