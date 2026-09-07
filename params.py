"""The parameter file, loaded.

    import params
    P = params.load()           # auto_tdmt.cfg next to this file (or $AUTO_TDMT_CFG)
    P.invert.depthWindowKm      # -> 30.0
    P.station.periodsMediumS    # -> [(10.0, 50.0), (20.0, 50.0)]
    P.as_dict()                 # every resolved value, for solution.json provenance

Format (MintPy style): one `auto_tdmt.<task>.<name> = <value>  # comment`
per line; `auto` means the default in DEFAULTS below. The trailing comment
in the cfg documents the default and its citation; THIS table is what the
code uses, so the two must agree (tests/test_rules.py checks every cfg
key exists here and vice versa).

Types are inferred from the default: a float stays a float, an int an
int, a string a string, a list of numbers a list, and "pairs" (a list of
2-lists, e.g. filter periods) is written in the cfg as `10 50, 20 50`.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from types import SimpleNamespace

CFG_PATH = Path(os.environ.get(
    "AUTO_TDMT_CFG", Path(__file__).resolve().parent / "auto_tdmt.cfg"))

# ---------------------------------------------------------------------------
# Defaults: the value "auto" resolves to. Keep the cfg comments in step.
# ---------------------------------------------------------------------------
DEFAULTS: dict[str, dict[str, object]] = {
    "geonet": {
        "minPrelimMag": 3.7,
        "maxDepthKm": 50.0,
        "placeholderDepthsKm": [5.0, 12.0, 33.0],
        "bbox": [-50.5, -33.0, 164.0, 182.5],
    },
    "station": {
        "minDistKm": 20.0,
        "maxDistKm": [120.0, 180.0, 250.0, 300.0],
        "maxStations": 16,
        "sectors": 8,
        "snrMin": 2.0,
        "snrWindowS": 200.0,
        "amplitudeMaxRatio": 3.0,
        "magBreaks": [4.5, 5.5],
        "periodsSmallS": [[10.0, 50.0]],
        "periodsMediumS": [[10.0, 50.0], [20.0, 50.0]],
        "periodsLargeS": [[20.0, 100.0], [30.0, 100.0]],
        "filterCorners": 3,
        "groupVelKms": 2.5,
        "windowTailS": [30.0, 40.0, 60.0],
        "windowMinS": 60.0,
        "maxWindowS": 200.0,
        "windowExtendMaxS": 60.0,
    },
    "invert": {
        "depthWindowKm": 30.0,
        "depthUncPct": 10.0,
        "degree": 5,
        "weight": "distance",
        "components": "ZRT",
        "correlate": 1,
        "stationVRDrop": 0.0,
        "stationVRFloor": 25.0,
        "maxTimeShiftS": 8.0,
        "maxRounds": 8,
        "maxDropPerRound": 3,
        "minStations": 3,
        "jackknife": 1,
        # INGV grade table, columns = N stations 3 | 4 | 5-8 | >8
        "gradeA_VR": [999.0, 60.0, 60.0, 50.0],
        "gradeB_VR": [70.0, 40.0, 40.0, 30.0],
        "gradeC_VR": [20.0, 20.0, 15.0, 15.0],
        "dcMinPublish": 60.0,
        "jackknifeMaxDeg": 25.0,
    },
    "forward": {
        "shearModulusPa": 30e9,
        "poissonRatio": 0.25,
        "wcLength": [-2.44, 0.59],
        "wcWidth": [-1.01, 0.32],
        "gridHalfwidthKm": "auto",
        "gridPoints": 120,
    },
    "publish": {
        "minGrade": "B",
        "minMw": 5.0,
        "minDisplacementM": 0.01,
        "maxPerDay": 3,
        "aftershockRadiusKm": 75.0,
        "aftershockWindowDays": 14.0,
        "aftershockMwMargin": 0.5,
    },
}

_LINE = re.compile(r"^\s*auto_tdmt\.(\w+)\.(\w+)\s*=\s*([^#]*?)\s*(#.*)?$")


def _coerce(text: str, default):
    """Parse a cfg value using the default's type as the template."""
    if text.lower() == "auto":
        return default
    if isinstance(default, list) and default and isinstance(default[0], list):
        # pairs: "10 50, 20 50"
        return [[float(x) for x in chunk.split()] for chunk in text.split(",")]
    if isinstance(default, list):
        return [type(default[0])(x) for x in text.replace(",", " ").split()]
    if isinstance(default, bool):
        return text.lower() in ("1", "true", "yes")
    if isinstance(default, int):
        return int(float(text))
    if isinstance(default, float):
        return float(text)
    return text  # string (weight, components, minGrade, "auto" halfwidth)


class Params(SimpleNamespace):
    """Nested namespace with the helpers the tasks need."""

    def as_dict(self) -> dict:
        return {task: dict(vars(ns)) for task, ns in vars(self).items()
                if isinstance(ns, SimpleNamespace)}

    # --- derived helpers (magnitude-binned lookups) -------------------------
    def mag_bin(self, prelim_mag: float) -> int:
        """0 small, 1 medium, 2 large (station.magBreaks)."""
        lo, hi = self.station.magBreaks
        return 0 if prelim_mag < lo else (1 if prelim_mag < hi else 2)

    def bands_hz(self, prelim_mag: float) -> list[tuple[float, float]]:
        """Ordered filter-band preference for this magnitude, in Hz."""
        periods = (self.station.periodsSmallS, self.station.periodsMediumS,
                   self.station.periodsLargeS)[self.mag_bin(prelim_mag)]
        return [(1.0 / long_s, 1.0 / short_s) for short_s, long_s in periods]

    def max_dist_km(self, prelim_mag: float) -> float:
        d = self.station.maxDistKm
        if prelim_mag < 4.0:
            return d[0]
        if prelim_mag < 4.5:
            return d[1]
        if prelim_mag < 5.0:
            return d[2]
        return d[3]

    def min_dist_km(self, prelim_mag: float) -> float:
        return 10.0 if prelim_mag < self.station.magBreaks[0] \
            else self.station.minDistKm

    def window_tail_s(self, prelim_mag: float) -> float:
        return self.station.windowTailS[self.mag_bin(prelim_mag)]

    def grade_vr_column(self, n_stations: int) -> int:
        """Column of the INGV table for this station count."""
        if n_stations <= 3:
            return 0
        if n_stations == 4:
            return 1
        if n_stations <= 8:
            return 2
        return 3


def load(path: Path | None = None) -> Params:
    """Read the cfg; unknown keys are an error (a typo must not silently
    fall back to a default), missing keys take the default."""
    path = Path(path) if path else CFG_PATH
    values = {task: dict(defaults) for task, defaults in DEFAULTS.items()}
    if path.exists():
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            m = _LINE.match(line)
            if not m:
                continue
            task, key, raw, _ = m.groups()
            assert task in DEFAULTS and key in DEFAULTS[task], (
                f"{path.name}:{lineno}: unknown parameter "
                f"auto_tdmt.{task}.{key}")
            values[task][key] = _coerce(raw, DEFAULTS[task][key])
    P = Params(**{task: SimpleNamespace(**v) for task, v in values.items()})
    P.source = str(path) if path.exists() else "defaults (no cfg found)"
    return P


def cfg_keys(path: Path | None = None) -> set[tuple[str, str]]:
    """(task, key) pairs present in the cfg file — for the consistency test."""
    path = Path(path) if path else CFG_PATH
    return {(m.group(1), m.group(2)) for m in
            (_LINE.match(l) for l in path.read_text().splitlines()) if m}


if __name__ == "__main__":
    import json
    P = load()
    print(f"# {P.source}")
    print(json.dumps(P.as_dict(), indent=2))
