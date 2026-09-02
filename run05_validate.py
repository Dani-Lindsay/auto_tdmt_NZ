"""Compare archived automated solutions against (1) the published NZ
regional CMT solutions (Ristau catalogue, GeoNet/data repository) and
(2) the Global CMT catalogue (Ekstrom et al., globalcmt.org; monthly +
quick NDK feeds), with GeoNet as the source of every original hypocentre.

    pixi run python run05_validate.py

Writes <OUTPUT_BASE>/validation/comparison.csv and comparison.jpg plus a
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


def dc_tensor(strike: float, dip: float, rake: float) -> np.ndarray:
    """Unit-moment DC tensor in NED (Aki & Richards 1980 eqs 4.84-4.89)."""
    p, d, r = np.radians([strike, dip, rake])
    sp, cp = np.sin(p), np.cos(p)
    sd, cd = np.sin(d), np.cos(d)
    sr, cr = np.sin(r), np.cos(r)
    s2p, c2p = np.sin(2 * p), np.cos(2 * p)
    s2d, c2d = np.sin(2 * d), np.cos(2 * d)
    mnn = -(sd * cr * s2p + s2d * sr * sp**2)
    mee = sd * cr * s2p - s2d * sr * cp**2
    mdd = -(mnn + mee)
    mne = sd * cr * c2p + 0.5 * s2d * sr * s2p
    mnd = -(cd * cr * cp + c2d * sr * sp)
    med = -(cd * cr * sp - c2d * sr * cp)
    return np.array([[mnn, mne, mnd], [mne, mee, med], [mnd, med, mdd]])


def tensor_angle_deg(a, b) -> float:
    """Angle between two DC tensors (plane-choice independent)."""
    m1, m2 = dc_tensor(*a), dc_tensor(*b)
    cos = np.sum(m1 * m2) / (np.linalg.norm(m1) * np.linalg.norm(m2))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def main() -> None:
    gcmt = load_gcmt_2026()
    print(f"Global CMT: {len(gcmt) if gcmt else 0} events loaded")
    cmt = load_geonet_cmt()
    id_col = "PublicID" if "PublicID" in cmt.columns else "EVENT_ID"
    cmt = cmt.set_index(cmt[id_col].astype(str))
    print(f"reference catalogue: {len(cmt)} solutions "
          f"(columns: {list(cmt.columns)[:12]}...)")

    rows = []
    for path in sorted(config.EVENTS_DIR.glob("*/solution.json")):
        sol = json.loads(path.read_text())
        ev = sol["event"]
        pid = ev["public_id"]
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
                "tensor_angle_deg": round(tensor_angle_deg(
                    (p1["strike"], p1["dip"], p1["rake"]),
                    (float(ref["strike1"]), float(ref["dip1"]),
                     float(ref["rake1"]))), 1),
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
                "tensor_angle_deg": round(tensor_angle_deg(
                    (p1["strike"], p1["dip"], p1["rake"]),
                    (np1.strike, np1.dip, np1.rake)), 1),
            })

    assert rows, ("no common events with the reference catalogue yet — "
                  "rerun after it next updates")
    df = pd.DataFrame(rows).sort_values("PublicID")
    out_dir = config.OUTPUT_BASE / "validation"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "comparison.csv", index=False)

    print(df.to_string(index=False))
    for ref, sub in df.groupby("reference"):
        print(f"\n[{ref}] n = {len(sub)}")
        print(f"  Mw:    mean dMw {sub.dMw.mean():+.2f}, |dMw| median "
              f"{sub.dMw.abs().median():.2f}")
        print(f"  depth: mean dZ {sub.dDepth.mean():+.1f} km, |dZ| median "
              f"{sub.dDepth.abs().median():.1f} km")
        print(f"  mechanism: median tensor angle "
              f"{sub.tensor_angle_deg.median():.0f} deg")

    colors = {"NZ_CMT_Ristau": "#0173B2", "GlobalCMT": "#DE8F05"}
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    for ref, sub in df.groupby("reference"):
        c = colors.get(ref, "black")
        axes[0, 0].scatter(sub.ref_Mw, sub.our_Mw, c=c, s=25, label=ref)
        axes[0, 1].scatter(sub.ref_depth, sub.our_depth, c=c, s=25)
        axes[1, 0].hist(sub.tensor_angle_deg, bins=np.arange(0, 121, 10),
                        color=c, alpha=0.6, label=ref)
        axes[1, 1].scatter(sub.ref_Mw, sub.dMw, c=c, s=25)
    axes[0, 0].plot([3.5, 6.5], [3.5, 6.5], "-", color="0.7", zorder=0)
    axes[0, 0].set_xlabel("published Mw")
    axes[0, 0].set_ylabel("auto Mw")
    axes[0, 0].legend(fontsize=8)
    lim = max(df.ref_depth.max(), df.our_depth.max()) + 5
    axes[0, 1].plot([0, lim], [0, lim], "-", color="0.7", zorder=0)
    axes[0, 1].set_xlabel("published depth (km)")
    axes[0, 1].set_ylabel("auto depth (km)")
    axes[1, 0].set_xlabel("DC tensor angle (deg)")
    axes[1, 0].set_ylabel("events")
    axes[1, 1].axhline(0, color="0.7")
    axes[1, 1].set_xlabel("published Mw")
    axes[1, 1].set_ylabel("auto - published Mw")
    fig.suptitle(f"auto_tdmt_NZ vs published NZ CMT solutions (n={len(df)})")
    fig.tight_layout()
    fig.savefig(out_dir / "comparison.jpg", dpi=150)
    print(f"\nwrote {out_dir}/comparison.csv and comparison.jpg")


if __name__ == "__main__":
    main()
