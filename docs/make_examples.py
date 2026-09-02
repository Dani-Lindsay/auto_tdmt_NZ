"""Build docs/examples/ from archived Fiordland solutions: copy the outward
figures and write a results table into docs/examples/README.md.

    pixi run python docs/make_examples.py 2026p660242 2026p660160 ...
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402

OUT = Path(__file__).resolve().parent / "examples"


def main(event_ids: list[str]) -> None:
    OUT.mkdir(exist_ok=True)
    rows = []
    for pid in event_ids:
        ev_dir = config.EVENTS_DIR / pid
        sol = json.loads((ev_dir / "solution.json").read_text())
        pref = sol["preferred"]
        q = sol["quality"]
        d = sol["publish_decision"]
        p1 = pref["plane1"]
        rows.append(
            f"| {pid} | M{sol['event']['prelim_mag']:.1f} / "
            f"{sol['event']['depth_km']:g} km | "
            f"**{pref['mw']:.2f}** / {pref['depth_km']:g} km | "
            f"{pref['vr']:.0f}% | {pref['pdc']:.0f}/{pref['pclvd']:.0f} | "
            f"{p1['strike']:.0f}/{p1['dip']:.0f}/{p1['rake']:.0f} | "
            f"{q['n_stations_used']} | "
            f"{sol['chosen_band'].replace('band_', '')} | "
            f"{sol['forward_model']['peak_abs_cm']:.2f} cm | "
            f"{'YES' if d['publish'] else 'no — ' + d['reasons'][0][:60]} |"
        )
        for name in ("outward_figure.jpg", "share_figure.jpg"):
            src = ev_dir / name
            if src.exists():
                shutil.copy(src, OUT / f"{pid}_{name}")
        depth_figs = sorted((ev_dir / sol["chosen_band"]).glob("depth.*.jpg"))
        if depth_figs:
            shutil.copy(depth_figs[0], OUT / f"{pid}_depth_search.jpg")

    header = (
        "# Worked examples — Fiordland sequence, 2026-09-02\n\n"
        "Four events near Milford Sound at the southern termination of the\n"
        "Alpine Fault, processed automatically with the Ristau (2008) South\n"
        "Island model. GeoNet initial depths were 5 km placeholders for all\n"
        "four. Full archives (solution.json with complete provenance, all\n"
        "figures, per-band results) are under `events/<publicID>/`.\n\n"
        "| publicID | GeoNet prelim M / depth | our Mw / depth | VR | "
        "%DC/%CLVD | plane1 s/d/r | stations | band | pred. peak disp | "
        "email gate |\n"
        "|---|---|---|---|---|---|---|---|---|---|\n"
    )
    footer = (
        "\n`<publicID>_outward_figure.jpg` is the composite the email\n"
        "carries: mttime waveform fits with the Deviatoric = DC + CLVD\n"
        "decomposition, then station map + Okada-predicted displacement +\n"
        "NISAR pass table.\n"
    )
    (OUT / "README.md").write_text(header + "\n".join(rows) + "\n" + footer)
    print(f"wrote {OUT / 'README.md'} and {len(list(OUT.glob('*.jpg')))} figures")


if __name__ == "__main__":
    assert len(sys.argv) > 1, "pass publicIDs"
    main(sys.argv[1:])
