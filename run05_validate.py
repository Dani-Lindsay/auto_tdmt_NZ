"""Compare archived automated solutions against the published NZ regional
CMT solutions (Ristau catalogue, GeoNet/data repository).

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
    cmt = load_geonet_cmt()
    id_col = "PublicID" if "PublicID" in cmt.columns else "EVENT_ID"
    cmt = cmt.set_index(cmt[id_col].astype(str))
    print(f"reference catalogue: {len(cmt)} solutions "
          f"(columns: {list(cmt.columns)[:12]}...)")

    rows = []
    for path in sorted(config.EVENTS_DIR.glob("*/solution.json")):
        sol = json.loads(path.read_text())
        pid = sol["event"]["public_id"]
        if pid not in cmt.index:
            continue
        ref = cmt.loc[pid]
        pref = sol["preferred"]
        p1 = pref["plane1"]
        angle = tensor_angle_deg(
            (p1["strike"], p1["dip"], p1["rake"]),
            (float(ref["strike1"]), float(ref["dip1"]), float(ref["rake1"])),
        )
        rows.append({
            "PublicID": pid,
            "our_Mw": round(pref["mw"], 2),
            "ref_Mw": float(ref["Mw"]),
            "dMw": round(pref["mw"] - float(ref["Mw"]), 2),
            "our_depth": pref["depth_km"],
            "ref_depth": float(ref["CD"]),
            "dDepth": round(pref["depth_km"] - float(ref["CD"]), 1),
            "tensor_angle_deg": round(angle, 1),
            "our_VR": round(pref["vr"], 1),
            "ref_VR": float(ref["VR"]) if "VR" in ref else np.nan,
            "our_DC": round(pref["pdc"], 0),
            "grade": sol["quality"].get("grade", "?"),
            "n_stations": sol["quality"]["n_stations_used"],
        })

    assert rows, ("no common events with the reference catalogue yet — "
                  "rerun after it next updates")
    df = pd.DataFrame(rows).sort_values("PublicID")
    out_dir = config.OUTPUT_BASE / "validation"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "comparison.csv", index=False)

    print(df.to_string(index=False))
    print(f"\nn = {len(df)}")
    print(f"Mw:    mean dMw {df.dMw.mean():+.2f}, |dMw| median "
          f"{df.dMw.abs().median():.2f}")
    print(f"depth: mean dZ {df.dDepth.mean():+.1f} km, |dZ| median "
          f"{df.dDepth.abs().median():.1f} km")
    print(f"mechanism: median tensor angle "
          f"{df.tensor_angle_deg.median():.0f} deg")

    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    ax = axes[0, 0]
    ax.plot([3.5, 6.5], [3.5, 6.5], "-", color="0.7")
    ax.scatter(df.ref_Mw, df.our_Mw, c="black", s=25)
    ax.set_xlabel("published Mw")
    ax.set_ylabel("auto Mw")
    ax = axes[0, 1]
    lim = max(df.ref_depth.max(), df.our_depth.max()) + 5
    ax.plot([0, lim], [0, lim], "-", color="0.7")
    ax.scatter(df.ref_depth, df.our_depth, c="black", s=25)
    ax.set_xlabel("published depth (km)")
    ax.set_ylabel("auto depth (km)")
    ax = axes[1, 0]
    ax.hist(df.tensor_angle_deg, bins=np.arange(0, 121, 10), color="0.4")
    ax.set_xlabel("DC tensor angle (deg)")
    ax.set_ylabel("events")
    ax = axes[1, 1]
    ax.scatter(df.ref_Mw, df.dMw, c="black", s=25)
    ax.axhline(0, color="0.7")
    ax.set_xlabel("published Mw")
    ax.set_ylabel("auto - published Mw")
    fig.suptitle(f"auto_tdmt_NZ vs published NZ CMT solutions (n={len(df)})")
    fig.tight_layout()
    fig.savefig(out_dir / "comparison.jpg", dpi=150)
    print(f"\nwrote {out_dir}/comparison.csv and comparison.jpg")


if __name__ == "__main__":
    main()
