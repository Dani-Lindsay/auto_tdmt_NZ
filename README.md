# auto_tdmt_NZ

**A personal tool for rapid InSAR response.** When a New Zealand earthquake
happens, the question I need answered within minutes is: *has this event
likely produced measurable surface displacement, and do I need to act as an
InSAR scientist* (task acquisitions, prepare processing, look at the next
NISAR/Sentinel-1 pass)? Answering that requires a moment tensor, so this
pipeline watches the GeoNet quake feed, runs a Dreger-style time-domain
moment tensor inversion ([LLNL mttime](https://github.com/LLNL/mttime))
with CPS Green's functions and the Ristau (2008) NZ velocity models,
forward-models the predicted surface displacement for both nodal planes
(Okada), reports NISAR acquisition timing over the epicentre, and emails
the result to a small list. The moment tensor solutions are a useful
by-product for other scientists; the displacement field is the point.

Event detection and every original hypocentre (time, location, preliminary
magnitude, initial depth) come from GeoNet; this project adds the moment
tensor, the revised centroid depth, and the displacement forecast on top
of that origin, and records both the GeoNet values and the revisions in
[`events/catalogue.csv`](events/catalogue.csv).

This is a personal, external project by Danielle Lindsay, not an
operational product of any agency, and it makes no representation about
any organisation's internal systems.

**All solutions are PRELIMINARY, deviatoric-only, and produced without
human review** — do not interpret mechanisms in volcanic/geothermal
settings from these solutions.

<p align="center"><img src="events/solutions_map.jpg" width="480" alt="All automated moment tensor solutions to date: beachballs sized by Mw (solid = grade A/B, washed = C/D) over the NZ Active Faults Database"></p>

## How it works

```
GeoNet quake API (poll, 10 min cron)
  -> processing floor (prelim M >= 3.7, NZ bbox)          run01_watch.py
  -> waveforms: GeoNet NRT FDSN, NZ broadbands <= 400 km   waveforms.py
     response removal -> ZRT -> bandpass -> 1 sps SAC
  -> Green's functions: precomputed CPS library            greens.py
     (Ristau 2008 North/South Island models, 10-500 km,
      depths 2-58 km, 10 Herrmann fundamental sources)
  -> mttime deviatoric inversion, depth search,            invert.py
     station selection by FIT (the funnel: survey the whole
     pool, keep the majority, build a clean core, then let
     every station earn its seat back), grid-edge depth
     guard, ordered filter-band preference
  -> quality gates (fit, no passengers, stability, depth,  invert.py
     azimuth pair) -> letter grade, or "no coherent
     solution" when nothing fits
  -> Okada forward model (both nodal planes,               okada_forward.py
     Wells & Coppersmith dimensions) -> predicted peak
     surface displacement
  -> NISAR last/next pass at epicentre (CMR)               nisar_dates.py
  -> figures (matplotlib/cartopy + mttime fits)            figure.py
  -> publication gates (our Mw >= 5.0 OR predicted         trigger.py
     displacement >= 1 cm; aftershock throttle; daily cap)
  -> email to the list (SMTP secrets)                      publish.py
```

Solutions, figures, and provenance are committed to `events/` — the repo is
the public archive — and [`events/catalogue.csv`](events/catalogue.csv) (column reference:
[`events/CATALOGUE_README.md`](events/CATALOGUE_README.md)) is
regenerated from the archived solutions after every event: one row per
solution with origin, nodal planes, Mw/Mo, centroid depth, %DC/%CLVD, VR and
MT elements (1e20 dyne-cm, matching the GeoNet CMT catalogue conventions). Every `solution.json` records the velocity model, GF
version, package versions, stations used/dropped (with reasons), the full
depth-search table and the band search.

## Students: the Human Review catalogue

Alongside the automated archive lives a **human-reviewed catalogue**
([events_human/](events_human/README.md)) built gradually by students —
final-year undergraduates and masters students learning regional moment
tensors by re-examining the automated solutions.

Start here: **[human_review.ipynb](human_review.ipynb)** — a Jupyter
notebook (runs in the browser or in VS Code) that walks you from
installation through the whole workflow: browse the automated
catalogue, watch a seismogram travel from raw counts to
inversion-ready displacement, run the inversion yourself with full
manual control (stations, time shifts, depth, record length, filter
band, velocity model), decide whether you can beat the machine, and
submit your reviewed solution by Pull Request.

Before opening the notebook, read [docs/METHOD.md](docs/METHOD.md)
(how the pipeline works) and
[docs/REVIEW_LEARNINGS.md](docs/REVIEW_LEARNINGS.md) (the reviewer
watch-list of known failure modes — read it twice). Multiple reviews of
the same event by different authors are welcome: each is its own row in
the human catalogue, and disagreement between reviewers is itself
useful information. The automated archive `events/` is read-only for
reviewers; PRs may only touch `events_human/`.

> The student workflow was *Made by Claude* and lightly tested (one
> full end-to-end run) — expect to troubleshoot, and treat that as
> part of the course.

## Local setup (macOS / Linux)

```sh
pixi install
pixi run test
# one-time Green's function library build (needs CPS):
#   brew install gcc && download+build CPS NP330 (see docs), then
pixi run python greens.py --build
# process one event:
pixi run python run02_process.py --event 2026p660242 --debug
```

`--debug` writes stage-by-stage troubleshooting figures (raw counts,
displacement, filtered ZRT record sections, station map) under
`<event>/<band>/diagnostics/`.

Outputs default to `~/work/proj_tdmt_NZ`; override with `AUTO_TDMT_OUTPUT`,
`AUTO_TDMT_EVENTS`, `AUTO_TDMT_GF` (CI points these into the checkout).

## CI (GitHub Actions)

- `watch.yml` — cron every 10 min: poll, process new events, email passing
  solutions, commit results.
- `process.yml` — manual reprocess of one publicID (with optional debug).
- `publish.yml` — manual email of one processed event (`force` to override
  gates).

The GF libraries are attached to the `gf-latest` release as
`gf_library.tar.zst` (containing `gf_cache/<model>/v1/...`) and cached in CI;
CPS never runs in CI. Rebuild + re-upload after any velocity model change:

```sh
pixi run python greens.py --build
cd ~/work/proj_tdmt_NZ && cp -r gf_library gf_cache \
  && tar --zstd -cf gf_library.tar.zst gf_cache && rm -r gf_cache
gh release create gf-latest gf_library.tar.zst --notes "GF libraries"
```

Secrets required for email: `SMTP_HOST`, `SMTP_PORT`, `SMTP_USER`,
`SMTP_PASS`, `MAIL_FROM`, `MAIL_TO` (ideally a single Google Group address so
subscriber addresses never live in this public repo).

## Velocity models

`models/nz_{north,south}_ristau2008.d` — Ristau (2008), SRL 79(3) Table 1,
doi:10.1785/gssrl.79.3.400 — so solutions are directly comparable to the
published NZ regional CMT solutions
([GeoNet/data moment-tensor](https://github.com/GeoNet/data/tree/main/moment-tensor);
their MT elements are in 1e20 dyne-cm).

## Quality grades — how a solution is rated

A grade is a statement about EVIDENCE, not about how much data went in.
Station count and azimuthal gap are deliberately not thresholds: three
well-fitting stations spanning 90 degrees make a good solution (standard
BSL practice), while ten stations carrying a passenger do not. Each
grade is the first row whose conditions all hold, top down:

| Grade | Fit (VR) | %DC | Every station fits | Mechanism stability | Depth | Meaning |
|---|---|---|---|---|---|---|
| **A** | ≥ 70% | ≥ 60 | min own VR ≥ 40 | jackknife rotation ≤ 15° (required) | interior, agrees with GeoNet | Publishable as-is |
| **B** | ≥ 60% | ≥ 60 | min own VR ≥ 25 | rotation ≤ 25° or jackknife not possible | interior, agrees with GeoNet | Publishable |
| **C** | ≥ 50% | — | min own VR ≥ 10 | — | — | Indicative only — archived, never emailed |
| **D** | below any C bar, or no two stations ≥ 90° apart | | | | | Archive/diagnostic material |
| **X** | no coherent solution — see below | | | | | No mechanism is reported at all |

Grade B is exactly the BSL publishability rule (VR ≥ 60 **and** DC ≥ 60)
plus the evidence checks. The individual criteria mean:

- **min own VR** — the worst-fitting station in the solution. A single
  "passenger" that the mechanism does not explain is a reason to
  distrust the whole answer, however good the total.
- **jackknife rotation** — the largest mechanism change when any one
  station is removed ([minimum rotation angle](validation/README.md)).
  A solution that depends on one station is not a solution.
- **depth** — the pick must not be a single-point maximum at the edge of
  the Green's-function depth grid (a known artifact), and the centroid
  depth must sit within 8 km of a real GeoNet hypocentre. The depth
  search itself is never bounded by GeoNet: it always covers the full
  grid, and the check only flags and downgrades the result.
- **90° azimuth pair** — the minimum geometry that can resolve a
  mechanism at all.

Only A/B solutions can pass the publication gate. The grades are
validated against independent catalogues
([validation/](validation/README.md)).

## No coherent solution

Some earthquakes cannot be inverted from this network: too small, too
far offshore, or arriving inside the long-period coda of a larger event.
Rather than publish a mechanism fitted to noise, the pipeline archives
those events with `"status": "no_coherent_solution"` — the full station
ledger and every pass's evidence, but no mechanism, magnitude or depth.
They appear in the catalogue with grade `X` and empty solution columns,
and are excluded from validation statistics.

## Operating thresholds

| Stage | Rule |
|---|---|
| **Triggers** | GeoNet preliminary magnitude >= 4.0, inside the NZ box (33-50.5 S, 164 E-177.5 W), event type "earthquake", depth <= 30 km (GeoNet fixed placeholder depths 5/12/33 km are exempt: true depth unknown, the depth search decides) |
| **Station selection** | NZ broadbands (HH? preferred over BH?), near-field exclusion 10 km (<M4.5) / 20 km, magnitude-scaled radius 120-300 km with a one-shot +100 km extension when the pool is thin. NOTHING USABLE IS PRE-FILTERED OUT: only no-data, dead channels (peak/noise < 1.2) and broken-response amplitude outliers (> 8x the network median, high side only) are removed. Near-field, weak-signal and cluster-surplus stations are TAGGED and must earn a seat by fit. Then the funnel: pass 1 inverts the whole pool with NO time shifts (so nothing can slide into a chance alignment) and ranks each station by the median of its own VR across the depth plateau; the best ~10 survive, plus any station filling an empty azimuth sector; pass 2 re-searches with shifts bounded at 8 s; the core (3-6 stations, forced to span >= 90 deg) gets its own depth search as the clean reference; every other station is then added at the core depth and kept if its own VR reaches 30 (20 for a sparse core, 10 if it fills an empty sector) without costing more than 3 joint VR points; finally a full search and an unconditional cull of any station fitting worse than silence |
| **Depth search** | 1-5 km at 0.5 km, to 10 km at 1 km, to 30 km at 2 km, to 58 km at 4 km; ALWAYS the full grid — the search is never bounded by GeoNet. A single-point VR maximum at the edge of the grid is rejected as an artifact in favour of the best interior local maximum, and the resulting depth is flagged (and capped at grade C) if it sits more than 8 km from a real GeoNet hypocentre |
| **Filter bands** | < M4.5: 10-50 s only; M4.5-5.5: 10-50 s then 20-50 s; >= M5.5: 20-100 s then 30-100 s. The menu is an ORDERED PREFERENCE at every magnitude — the first band that passes its gates wins, because VR is not comparable across bands (a longer period is smoother and scores higher even when fitting noise). A band is escalated when the inverted Mw overshoots the preliminary magnitude by >= 0.6 |
| **Rated** | see the grade rubric above: evidence (VR, %DC, worst station's own VR, jackknife stability, depth plausibility, a 90 deg azimuth pair), not station count or azimuthal gap |
| **No coherent solution** | when nothing coheres (survey and majority VR both < 20, no core, or a final VR < 20) the event is archived with `status: no_coherent_solution` and grade X — station ledger kept, no mechanism reported |
| **Publishes** | grade A or B AND (our Mw >= 5.0 OR Okada-predicted peak displacement >= 1 cm); max 3 emails/day; aftershock throttle (within 75 km/14 d of a published event, must be within 0.5 Mw of it or above the Mw gate) |
| **Preferred solution** | the depth whose VR is highest, with %DC breaking near-ties only (within 2 VR points); the wider 5-point plateau is reported as `Plateau_km` — how well the depth is resolved |

## Data sources

- GeoNet (Earth Sciences New Zealand) quake API + FDSN (NRT + archive):
  event detection, all original hypocentres, waveforms and station
  metadata. CC BY 3.0 NZ. Polling is polite: one request per 10-minute
  cron tick, gzip, descriptive User-Agent.
- NASA CMR for NISAR GSLC granule timing.
- NZ Active Faults Database (GNS Science) and GeoNet delta GNSS marks for
  map context.

## References

Software:

- Chiang, A. — **MTtime**, Time Domain Moment Tensor Inversion in Python,
  LLNL-CODE-814839, https://github.com/LLNL/mttime. Methodology after
  Dreger & Helmberger (1993), Dreger (2003) and Minson & Dreger (2008).
- Herrmann, R. B. (2013). Computer Programs in Seismology: an evolving
  tool for instruction and research. *Seism. Res. Lett.* 84, 1081-1088.
  https://rbherrmann.github.io/ComputerProgramsSeismology/
- Beyreuther, M., R. Barsch, L. Krischer, T. Megies, Y. Behr &
  J. Wassermann (2010). ObsPy: a Python toolbox for seismology.
  *Seism. Res. Lett.* 81(3), 530-533. https://github.com/obspy/obspy
- Jolivet, R. — **okada4py**, Python/C implementation of Okada (1992),
  https://github.com/jolivetr/okada4py
- Crameri, F. — Scientific colour maps,
  https://www.fabiocrameri.ch/colourmaps/ (via cmcrameri).

Method:

- Dreger, D. S., & D. V. Helmberger (1993). Determination of source
  parameters at regional distances with three-component sparse network
  data. *J. Geophys. Res.* 98, 8107-8125.
- Dreger, D. S. (2003). TDMT_INV: Time Domain Seismic Moment Tensor
  INVersion. *International Handbook of Earthquake and Engineering
  Seismology* 81B, 1627.
- Minson, S. E., & D. S. Dreger (2008). Stable inversions for complete
  moment tensors. *Geophys. J. Int.* 174, 585-592.
- Ristau, J. (2008). Implementation of routine regional moment tensor
  analysis in New Zealand. *Seism. Res. Lett.* 79(3), 400-415.
  doi:10.1785/gssrl.79.3.400 (velocity models, Table 1).
- Okada, Y. (1985). Surface deformation due to shear and tensile faults
  in a half-space. *Bull. Seism. Soc. Am.* 75(4), 1135-1154.
- Okada, Y. (1992). Internal deformation due to shear and tensile faults
  in a half-space. *Bull. Seism. Soc. Am.* 82(2), 1018-1040.
- Wells, D. L., & K. J. Coppersmith (1994). New empirical relationships
  among magnitude, rupture length, rupture width, rupture area, and
  surface displacement. *Bull. Seism. Soc. Am.* 84(4), 974-1002.
- Aki, K., & P. G. Richards (1980). *Quantitative Seismology*. W.H.
  Freeman (double-couple tensor construction used in validation).

This workflow was compiled with Claude (Anthropic) assistance under the
direction of Danielle Lindsay; the science stands on the shoulders of the
authors above — cite them, not this repository, for the methods.
