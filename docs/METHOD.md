# auto_tdmt_NZ — Method and Implementation

Documentation of the automated regional moment tensor processing flow, for
review. Danielle Lindsay (danielle.lindsay@earthsciences.nz), September 2026.

Companion to the [README](../README.md) (setup/operations). This document
covers the science choices, module by module, and worked examples from the
2026-09-02 Fiordland sequence.

## 1. Motivation and scope

This is a personal, external project. Its purpose is rapid InSAR response:
within minutes of a New Zealand earthquake, decide whether it has likely
produced measurable surface displacement and whether action is needed
(acquisition tasking, processing preparation, checking the next NISAR
pass). A moment tensor is the prerequisite for that displacement forecast,
so the pipeline produces PRELIMINARY automated solutions
minutes-to-hours after an event using the Dreger-style time-domain method
with 1-D Green's functions and the published Ristau (2008) NZ velocity
models; solutions are therefore directly comparable to the published NZ
regional CMT solutions
([GeoNet/data moment-tensor](https://github.com/GeoNet/data/tree/main/moment-tensor)),
which serve as the validation reference. The MT solutions are a useful
by-product for other scientists; the displacement field is the point.
The project makes no representation about any organisation's internal or
operational systems.

The inversion engine is **mttime** (Chiang, LLNL;
github.com/LLNL/mttime) — the Python implementation of Dreger's TDMT used
operationally in California — driven end-to-end by this repository. We do
not reimplement any inversion mathematics.

## 2. Processing flow

```
GeoNet quake API poll (cron, 10 min) ........................ run01_watch.py
  processing floor: prelim M >= 3.7, NZ bbox, not "deleted"     trigger.py
     |
event processing ............................................ run02_process.py
  1. event metadata (quake API)                                 geonet.py
  2. station inventory + waveforms (FDSN NRT/archive)           waveforms.py
  3. pre-processing -> SAC displacement (cm), ZRT, 1 sps        waveforms.py
  4. Green's functions from precomputed CPS library             greens.py
  5. mttime deviatoric inversion + depth search                 invert.py
     station selection by FIT (the funnel, 3.1), grid-edge
     depth guard, ordered filter-band preference
  6. quality gates -> letter grade, or "no coherent solution"   invert.py
  7. Okada forward model of predicted surface displacement      okada_forward.py
  8. NISAR pass timing at the epicentre (NASA CMR)              nisar_dates.py
  9. figures (mttime waveform fits + cartopy maps)              figure.py
 10. publication gate -> email                                  trigger.py, publish.py
```

Every solution is archived (`events/<publicID>/`) whether or not it is
published, with full provenance in `solution.json`.

## 3. Data selection and pre-processing (`waveforms.py`)

### 3.1 Station selection — the funnel

The goal is the best data for the job. The previous scheme tried to
reach it by filtering candidates out before the inversion; an audit of
693 archived events showed that failing badly — 61% of all station
exclusions were decided by the solution itself, and 69% of those vetoes
were issued by reference solutions whose own variance reduction was
below 20 (junk vetoing good data). The far-field 3x-depth rule alone
removed the closest station from every event it touched, and 94% of
those events ended grade C or D.

Selection v4 (2026-09) therefore deletes nothing usable. The data
decides, in a funnel:

**Stage A — the pool** (`waveforms.py`). Hard rejection only for
genuinely unusable records: no waveform, no response, gaps, the wrong
sample count, a dead channel (peak-to-noise below 1.2), or a
broken-response amplitude outlier. That last screen is ONE-SIDED, only
catching stations far ABOVE the network median: a station far below it
may simply lie near a nodal plane, where small amplitude is real
information about the mechanism rather than evidence of a bad station.
Everything else enters the pool carrying demotion TAGS — `near_field`
(inside 3x the source depth; the CPS Green's functions are
complete-wavefield, so the fit decides), `weak_signal` (peak/noise below
5), `cluster_surplus` (beyond two stations within 25 km, e.g. the
Ruapehu ring) — which shift the burden of proof without excluding
anybody.

**Stage B — the funnel** (`invert.py`, `invert_with_rejection`).

1. **Survey, everyone in, no time shifts.** The whole pool is inverted
   over a coarse depth grid with mttime's cross-correlation disabled.
   This matters: mttime's shift search is unbounded (a 60 s window can
   slide past 100 s), so a noise trace can always find a chance
   alignment and earn undeserved variance reduction. With shifts off,
   the ranking is honest by construction. Each station's evidence is the
   MEDIAN of its own VR across the contiguous VR plateau — chance
   alignment is depth-specific, real coherence is not.
2. **Keep the majority.** The best ten stations that fit at all
   (own VR >= 10), plus the best station in any azimuth sector not yet
   represented, up to twelve.
3. **Re-search, then build a clean core.** Shifts are on from here, and
   bounded: a station whose solved shift exceeds 8 s is rejected, its
   fit no longer being evidence. (Archive calibration: stations used in
   grade A/B solutions sit at |zcor| <= 9 s at the 95th percentile,
   while grade-D stations reach 0.88 of their whole travel time.) The
   core is the best three to six by own VR, forced to contain two
   stations at least 90 degrees apart — the minimum geometry that can
   resolve a mechanism — and it gets its own full depth search. That
   solution is the reference every other station is judged against, so
   it must be clean.
4. **Earn your seat back.** Every station outside the core — tagged and
   pruned ones alike — is added at the core depth and kept if the core
   solution predicts its waveform (own VR >= 30, relaxed to 20 for a
   sparse core, or 10 if it fills an empty azimuth sector) without
   costing more than 3 joint VR points. A `cluster_surplus` station must
   additionally improve the joint fit. %DC never enters the decision: a
   noise station can inflate it by dragging the tensor toward a generic
   mechanism.
5. **Final search and the anti-fitting cull.** A station whose own VR is
   negative fits worse than silence and only steers the tensor; it is
   removed unconditionally, with no joint-gain test and no protection
   for being alone in its sector, and the search is repeated.
6. **Jackknife.** Leave-one-station-out at the preferred depth
   quantifies how much any single station moves the answer, and feeds
   the grade.

There is deliberately no greedy "drop whatever raises the joint VR"
pass. It was evicting stations fitting at 54-65% to gain two points —
variance-reduction vanity paid for in azimuth coverage.

When nothing coheres — survey and majority VR both below 20, no viable
core, or a final VR below 20 — the event is archived as
**no coherent solution** rather than reporting a mechanism fitted to
noise (§6).

Every pool station's whole history is recorded in `solution.json`:
peak/noise, tags, its own VR at each pass, the solved time shift, the
admission verdict and reason. The per-band all-station waveform figure
shows every candidate's record with those numbers printed beside it, so
any exclusion can be judged against the data that produced it.

### 3.2 Pre-processing

- **Far-field guard**: a station inside 3x the source depth is TAGGED
  `near_field`, not excluded (§3.1). The Green's functions include the
  near-field terms, so the fit is better evidence than the heuristic.
- **Waveform windows**: origin−150 s to origin+230 s downloaded; final
  cut origin−30 s to origin+200 s.
- **Inversion record length**: distance-adaptive at every magnitude —
  each station'''s fitted window is 30 s pre-origin + distance/2.5 km/s +
  a magnitude-dependent tail (30 s below M4.5, 40 s to M5.5, 60 s
  above), clamped to 60–150 s. A close station'''s train is over quickly
  regardless of event size; fitting the empty tail only taxes VR. This
  is a deterministic kinematic cut, not an amplitude-based one, for
  reproducibility.
- **Response removal**: to displacement with pre-filter
  (0.004, 0.007, 10, 20) Hz, then rotation to ZNE and NE→RT along the
  great-circle back-azimuth.
- **Peak-to-noise measurement**: per station, median over the three
  components of peak|signal| / RMS(pre-event noise), signal window = the
  distance-adaptive inversion window above, noise window −115 to −2 s.
  This replaced an RMS-ratio SNR gate (threshold after Ristau 2008)
  that could not order stations the way visual inspection does: an
  emergent band-limited packet is a localised peak, and RMS over a
  200 s record dilutes it with empty tail. Dropped stations are
  recorded with their peak/noise in `solution.json`.
- **Filtering**: zero-phase 3-corner Butterworth, band from the menu in
  §5; 5% taper; decimate/resample to 1 sps; convert m→cm (TDMT
  convention). Green's functions receive the identical filter.
- Data are never padded, interpolated across gaps, or substituted: any
  deficiency drops the station with a recorded reason (fail loud).

The chain follows the mttime example notebooks (Chiang) and matches the
preparation described in Ristau (2008) §"Preparation of the observed
waveforms".

## 4. Green's functions (`greens.py`)

- **Velocity models** (`models/*.d`, model96 format): Ristau (2008), SRL
  79(3), Table 1 — North Island (8 layers) and South Island (5 layers,
  Moho 39 km) models with Qp=400/Qs=200, transcribed verbatim from the
  paper (doi:10.1785/gssrl.79.3.400). Model chosen per event by a
  North/South Island split at Cook Strait (`config.model_for_event`);
- **Computation**: CPS 3.30 (Herrmann) wavenumber integration, run once
  locally per model:
  `hprep96 -EQEX | hspec96 | hpulse96 -D -i | f96tosac -B`,
  producing the ten Herrmann fundamental-source time series
  (ZDD RDD ZDS RDS TDS ZSS RSS TSS ZEX REX) that mttime's
  `green="herrmann"` mode expects (displacement, cm, 1e20 dyne-cm source).
- **Library grid**: distances 10–500 km at 5 km spacing; source depths
  1.0–5.0 km at 0.5 km (fine where shallow depth discrimination
  happens), 6–10 km at 1 km, 12–30 km at 2 km, and 34–58 km at 4 km
  (Fiordland subduction events exceed crustal depths); dt = 1 s,
  npts = 256, vred = 0 (traces start at origin time). Runtime lookup
  takes the nearest grid distance (≤2.5 km error, absorbed by the
  per-station time-shift search).
- Libraries are stored unfiltered and band-passed per event to match the
  data exactly. A manifest records the model file's SHA-256, grid, and
  build time. CI never runs CPS — it downloads the library tarball.

## 5. Inversion (`invert.py`)

- **Engine**: mttime deviatoric inversion (degree 5), ZRT components,
  inverse-distance weighting, per-station cross-correlation time shifts
  (zcor; the pandas<3 pin is required for this — mttime issue #15).
- **Depth search**: the FULL library grid, always — the solution is kept
  fully independent of GeoNet's hypocentral depth, which is recorded (with
  its uncertainty) for comparison only. The catalogue reports the final
  pick alongside the VR-max and DC-max depths and the plateau width.
- **Filter-band menu** (BSL practice: a small menu of period bands, longer
  periods for larger events; the pipeline tries each and picks by the rule
  below):
  | preliminary M | candidate bands |
  |---|---|
  | < 4.5 | 10–50 s only (no coherent energy above ~20 s period; longer-period trials only ever fit noise and inflate Mw — other bands remain testable via `run02 --band`) |
  | 4.5–5.5 | 10–50 s, then 20–50 s (20–100 s pruned: won 1/19 events in this bin, elsewhere only fit noise) |
  | ≥ 5.5 | 20–100 s, then 30–100 s |

  The menu is an **ordered preference at every magnitude**: the first
  band that produces a solution passing its gates wins. VR must not
  arbitrate across bands — a longer-period band is smoother and posts a
  higher VR even when it is fitting filtered noise (2026p033598: 20–50 s
  scored VR 61 and inflated Mw by 0.22 over the visibly signal-fitting
  10–50 s at VR 40), and at the largest magnitudes the longest band is
  visibly over-smoothed. A band is escalated when the inverted Mw
  overshoots the preliminary magnitude by ≥ 0.6, which means the event is
  bigger than the menu assumed.
- **Station selection**: the funnel (§3.1). Stations whose data the
  final solution cannot predict are removed rather than being allowed to
  dilute the %DC (cf. dropping persistently low-VR stations, Ristau 2008;
  Dreger & Helmberger 1993); every rejection is recorded with its reason
  and the numbers behind it.
- **Preferred solution rule** (VR first, then DC as a tie-break):
  candidates are the depths on the contiguous VR plateau around the
  reference (a bimodal VR curve must not let a disconnected lobe steal
  the pick), and among them the highest %DC wins. The DC tie-break window
  is deliberately narrow — **2 VR points** — because DC may break a
  near-tie but must not buy a real loss of fit: at the old 5-point window
  the pick landed on the plateau EDGE (2026p091845 took 24 km at VR 60.5
  / DC 98 over 18–22 km at VR 65 / DC 88–93). Recomputing the pick from
  all 343 archived depth searches that have a Ristau reference depth,
  the median depth error improves from 9.0 km to 8.0 km when the window
  tightens, and DC adds nothing beyond that (a pure VR maximum also
  scores 8.0). The wider 5-point plateau is still reported as
  `Plateau_km`: how well the depth is resolved.
  Rationale for using DC at all: VR is a weak discriminator with depth
  (it often climbs monotonically) while spurious CLVD grows where the
  depth or model is wrong; pure VR-max picks produced solutions with
  implausible ~70% CLVD. Both the VR-max and %DC-max depths are recorded
  and a disagreement flag is set.
- **Grid-edge guard**: a VR maximum sitting on the first or last depth of
  the library with a single-point plateau is an artifact — the smoothest
  Green's functions at the ends of the grid absorb noise — so the best
  interior local maximum within 5 VR points is preferred instead
  (2026p508890 rode a 58 km edge at VR 22.7 over the physical 8 km peak
  at VR 18.3 with DC 88–96; 520779, 348732 and 300334 the same).

## 6. Quality gates and publication (`invert.py`, `trigger.py`)

Every solution carries a letter grade built from EVIDENCE, not from how
much data went in. Station count and azimuthal gap are deliberately not
thresholds: three well-fitting stations spanning 90 degrees make a good
solution (standard BSL practice; Ristau's own catalogue has a median of
7 stations with quartiles 4–11, and uses as few as 1–3), while ten
stations carrying a passenger do not.

| Grade | VR | %DC | min own VR | jackknife rotation | depth |
|---|---|---|---|---|---|
| A | ≥ 70 | ≥ 60 | ≥ 40 | ≤ 15° (required) | interior, within 8 km of GeoNet |
| B | ≥ 60 | ≥ 60 | ≥ 25 | ≤ 25° or not possible | interior, within 8 km of GeoNet |
| C | ≥ 50 | — | ≥ 10 | — | — |
| D | anything below C, or no two stations ≥ 90° apart | | | | |

Grade B is exactly the BSL publishability rule (VR ≥ 60 and DC ≥ 60)
plus the evidence checks. Thresholds come from the archive's own A/B
statistics: minimum own-station VR sits at the 25th percentile of 41
(A) and 35 (B), and jackknife rotation at the 90th percentile of 12°
(A). **min own VR** is the worst-fitting station in the solution — one
passenger the mechanism cannot explain is reason to distrust the whole
answer. **Jackknife rotation** is the largest mechanism change when any
one station is removed (§ Validation metric); a solution that depends on
a single station is not a solution.

**Depth plausibility.** The depth search is never bounded by GeoNet, but
the result is judged: a centroid depth more than 8 km from a real
(non-placeholder) GeoNet hypocentre is flagged and capped at grade C, so
it cannot be published as though it were fine. Calibration: in the
Ristau catalogue, centroid depths sit a median 4 km from the GeoNet
hypocentre, 57% within 5 km and 82% within 10 km — so a large
disagreement is claiming something a careful analyst catalogue rarely
does. GeoNet's fixed placeholder depths (5/12/33 km) are not
measurements and are never used to judge us.

**No coherent solution.** When nothing coheres — the survey and majority
inversions both below VR 20, no viable core, or a final VR below 20 —
the event is archived with `"status": "no_coherent_solution"`, grade X:
the full station ledger and every pass's evidence, but no mechanism,
magnitude or depth. Such events never publish and are excluded from
validation statistics. Publishing a mechanism fitted to noise would be
worse than admitting the network could not constrain the event.

Publication gate, applied to **our inverted Mw** (never GeoNet's mixed
ML-type preliminary magnitudes): publish if Mw ≥ 5.0 OR the Okada-predicted
peak surface displacement ≥ 1 cm (the InSAR-interesting case). Anti-spam:
max 3 emails/day; within 75 km and 14 days of an already-published event,
a smaller event must be within 0.5 Mw of it (or exceed the Mw gate) to
publish. Everything processed is archived regardless.

## 7. Deformation forward model (`okada_forward.py`)

For each nodal plane (the MT cannot distinguish them): rectangular uniform-
slip Okada (1992) dislocation via okada4py; length/width from Wells &
Coppersmith (1994, all-type regressions); slip from M0 = μLW·s with
μ = 30 GPa; centroid at the inverted depth (plane pushed down if its top
would breach the surface); displacement evaluated on a 121×121 km grid.
Calibration anchors: Mw 5.5 at 10 km ≈ 1.6 cm peak |u|; Mw 4.0 at 5 km ≈
0.4 mm (hence the displacement gate passes only shallow M≳4.7 events).

## 8. Outputs and figures (`figure.py`, `diagnostics.py`)

Per event: `solution.json` (all numbers + provenance), `draft_email.txt`,
mttime waveform-fit and depth-search figures, and three outward figures for
the email: (1) mttime's waveform-fit page, untouched (full Deviatoric =
DC + CLVD decomposition, per-station fits, VR/%DC); (2) the modelling
figure — Mercator map panels of the station geometry with the moment
tensor beachball drawn from the actual tensor elements (obspy/mopad
full-MT rendering, so any CLVD is represented honestly rather than
collapsed to the closest DC) and the Okada-predicted E/N/U displacement
over the modelled area (cmcrameri 'vik', modelled nodal plane outlined),
with the NISAR pass table and provenance line; (3) a depth-sensitivity
summary (VR, %DC/%CLVD, Mw vs depth) that also displays the
selection-rule window. All plotting is matplotlib (one stack; maps
via cartopy with ocean/land/coastline features and gridlines).

`--debug` additionally saves stage-by-stage QC figures per band: raw
counts, response-removed displacement, filtered ZRT record sections, and a
station-geometry map.

Every `solution.json` records: velocity model + GF version, mttime/obspy
versions, filter band, all depth-search rows, stations used (with SNR,
distance, azimuth) and dropped (with reasons), the band-search table,
quality gates/warnings, forward-model parameters, NISAR passes, and the
publication decision with reasons.

## 9. Worked examples: Fiordland sequence, 2026-09-02 (Milford Sound)

Four events at the southern termination of the Alpine Fault / Fiordland
subduction corner, all with GeoNet fixed 5 km placeholder depths. Processed
with the South Island model. See `docs/examples/` for the figures; numbers
below are from the final run (SNR ≥ 2.0).

*(Table completed from the archived solution.json files — see
docs/examples/README.md.)*

Key observations:
1. **Velocity model dominates solution quality.** With GIL7 (California)
   the M5.6 mainshock inverted to 25% DC; with the Ristau South Island
   model, 81% DC at the same VR — the spurious CLVD was model error.
2. **Depth recovery from placeholders.** Both M5.5+ events prefer ~8 km
   over the 5 km fixed depth, with clean single-peaked VR(depth) curves.
3. **ML–Mw offset.** GeoNet preliminary magnitudes exceed our Mw by
   0.4–0.6 units across the sequence — the recurring NZ ML–Mw discrepancy
   discussed by Ristau (2008) — which is why all gating uses our Mw.
4. **Gates catch the marginal cases.** The M4.8 aftershock lost all but
   2 stations to the SNR gate in its first-choice band and was blocked by
   the 3-station minimum rather than emailing a poorly-constrained
   solution.

## 10. Known limitations / review questions

- 1-D models; no offshore-specific model yet (Ristau 2008 anticipates
  offshore-region issues north of North Island and at Puysegur).
- North Island model transcribed but not yet exercised on a real event.
- The North/South model split is a crude Cook Strait rule.
- Depth grid starts at 2 km; very shallow events sit at the grid edge
  (flagged, not blocked).
- No isotropic term (degree 5); volcanic/geothermal events with real ISO
  components will be forced deviatoric.
- Single-pass station rejection; no jackknife uncertainty yet (the EPS207
  jackknife scheme is the natural next addition).
- Validation against the published NZ CMT solutions is implemented as a
  comparison hook but awaits their solutions for overlapping events
  (their CSV updates ~monthly; note their MT units are 1e20 dyne-cm).

### Verification aids

Nodal-plane and fault-geometry conventions were cross-checked numerically
(conjugate-plane far-field equivalence, rake/side acid tests — now anchor
tests in the suite) and visually against the USGS finite-fault event pages
(black-outlined plane, red up-dip edge convention, adopted here) and the
interactive Focal Mechanism Explorer at https://eq.comoglu.com/bb/.

## Validation metric

Mechanism agreement with reference catalogues (`run05_validate.py`) and
jackknife stability are measured as the **minimum rotation angle**
between two double couples (J. Townend, pers. comm. 2026-08-20):
each mechanism is expressed as a rotation matrix of its principal axes
with respect to geographic coordinates (Walsh et al. 2009, eqs 1-3),
and the angle is arccos((tr(R1^T R2) - 1)/2), minimised over the
double-couple symmetry group so the result is independent of which
nodal plane parameterises either mechanism (cf. Kagan 1991). The
approach follows Townend et al. (2012), supplement eq. 1
(`docs/Townend_etal_2012_supplement.pdf`).

## References

- Walsh, D., Arnold, R., Townend, J. (2009). A Bayesian approach to
  determining and parameterising earthquake focal mechanisms. GJI 176,
  235-255.
- Townend, J., et al. (2012). Three-dimensional variations in present-day
  tectonic stress along the Australia-Pacific plate boundary in New
  Zealand. EPSL, doi:10.1016/j.epsl.2012.08.003 (supplement in docs/).
- Kagan, Y. Y. (1991). 3-D rotation of double-couple earthquake sources.
  GJI 106, 709-716.
- Chiang, A. MTtime: Time Domain Moment Tensor Inversion in Python.
  LLNL-CODE-814839. github.com/LLNL/mttime
- Dreger, D., & Helmberger, D. (1993). Determination of source parameters
  at regional distances with three-component sparse network data. JGR 98.
- Herrmann, R. B. (2013). Computer Programs in Seismology. SRL 84.
- Okada, Y. (1992). Internal deformation due to shear and tensile faults
  in a half-space. BSSA 82.
- Ristau, J. (2008). Implementation of routine regional moment tensor
  analysis in New Zealand. SRL 79(3), 400–415. doi:10.1785/gssrl.79.3.400
- Wells, D. L., & Coppersmith, K. J. (1994). New empirical relationships
  among magnitude, rupture length, rupture width, rupture area, and
  surface displacement. BSSA 84.
