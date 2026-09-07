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
the result to a small list.

Event detection and every original hypocentre (time, location, preliminary
magnitude, initial depth) come from GeoNet; this project adds the moment
tensor, the revised centroid depth, and the displacement model on top, and 
records both the GeoNet values and the revisions in
[`catalogue.csv`](catalogue.csv). 

**All solutions are PRELIMINARY, deviatoric-only, and produced without
human review** — do not interpret mechanisms in volcanic/geothermal
settings from these solutions.

<p align="center"><img src="events/solutions_map.jpg" width="480" alt="All automated moment tensor solutions to date: beachballs sized by Mw (solid = grade A/B, washed = C/D) over the NZ Active Faults Database"></p>

## The five tasks

All code is in [`src/`](src/). The pipeline is five tasks, each a
module (or two) whose docstring opens with its task number:

```
 1  geonet    poll the quake API, apply the processing floor      geonet.py  trigger.py
 2  stations  select stations, download, process to SAC           waveforms.py
 3  invert    bound the depth grid, run mttime, the selection      invert.py  (greens.py stages the GFs)
              loop, jackknife, grade
 4  forward   Okada displacement for both nodal planes            okada_forward.py
 5  publish   figures, catalogue tables, publish decision, email  figure.py catalogue.py trigger.py publish.py
```

Each task has a walkthrough notebook — [`docs/task_1_geonet.ipynb`](docs/task_1_geonet.ipynb)
… [`docs/task_5_publish.ipynb`](docs/task_5_publish.ipynb) — that
imports the module and calls its functions on a real event, so you can
run and inspect any single step (`pixi run notebooks` rebuilds them).

## Entry points

The launchers carry no task number: they are how the tasks get run, not
tasks themselves. Each has a `pixi run` alias and a `--help`.

```
 watch.py           cron    task 1, then process_event and publish_event for every new event
 process_event.py   by hand one event, tasks 1 → 5:  --event <publicID> [--debug] [--band]
 publish_event.py   by hand email one processed event (task 5's last step; --force skips the gates)
 backsweep.py       by hand process_event over a date range; never emails
 validate.py        by hand the whole archive against the Ristau NZ CMT and GCMT catalogues
 sync_repo.py       by hand local archive → repo, rebuild the tables and validation, commit, push
 replay_labels.py   by hand acceptance test: replay the reviewer's station labels
```

`validate.py` is catalogue-level, not per event: `sync_repo.py` runs it
after every sync so the validation report is refreshed whenever the
archive is published.

## The rules and where they come from

| Task | Rule | Default | Source |
|---|---|---|---|
| 1 | processing floor | prelim M ≥ 3.7, depth ≤ 50 km (placeholders exempt) | Herrmann, Benz & Ammon (2011) reach Mw 3.7 in the same band |
| 2 | distance window | 120–300 km by magnitude; near limit 10–20 km | Gisola (Triantafyllis et al. 2022) |
| 2 | signal quality | per-component SNR = RMS(signal)/RMS(noise) over the inverted window; station usable when its best component ≥ 2 | BMKG (Halauwet et al. 2024); best-component rule so nodal components do not veto a station |
| 2 | broken responses | reject peak×distance > 3× network median, high side only | Duputel et al. (2012), one-sided so nodal stations survive |
| 2 | azimuth balance | rank by SNR, fill 8 sectors round-robin to 16 stations | Gisola, SCARDEC, INGV; relaxed for one-sided offshore geometry |
| 2 | filter bands | 10–50 s (M<4.5) · 10–50 then 20–50 s · 20–100 then 30–100 s (M≥5.5), ordered preference | Clinton et al. (2006); INGV |
| 2 | record window | 30 s before origin to distance/2.5 km/s + tail, envelope-extended ≤ 60 s, cap 200 s | Herrmann's group-velocity cut; the 2026 figure review |
| 3 | depth grid | GeoNet depth ± 30 km (full library for placeholders) | F-net ±30, Gisola ±31, W-phase/SCARDEC ±50 km |
| 3 | preferred depth | maximum VR (mttime's own rule); uncertainty = range within 10% of max | every system reviewed; Vallée et al. (2011), Bernardi et al. (2004) |
| 3 | station selection | invert all; drop \|shift\| > 8 s at once and own VR < 25 three per round, worst first; re-invert until stable | Clinton, Hauksson & Solanki (2006); NEIC SynDepth for the shift cap; gempa's fixed floor |
| 3 | stability | leave-one-out jackknife, minimum rotation angle | Fukuyama et al. (1998); Townend et al. (2012) |
| 3 | grade | INGV table (VR bar falls with station count) + DC ≥ 60 + rotation ≤ 25° for A/B | INGV TDMT quality legend; BSL; scisola |
| 5 | publish | A/B and (Mw ≥ 5.0 or predicted displacement ≥ 1 cm); 3/day; aftershock throttle | this project |

Every number in that table is one line in `auto_tdmt.cfg`.

## Quality grades

| N stations | D | C | B | A |
|---|---|---|---|---|
| 3 | VR < 20 | 20–70 | ≥ 70 | — |
| 4 | < 20 | 20–40 | 40–60 | ≥ 60 |
| 5–8 | < 15 | 15–40 | 40–60 | ≥ 60 |
| > 8 | < 15 | 15–30 | 30–50 | ≥ 50 |

A and B additionally require **%DC ≥ 60** (the BSL publishability rule)
and a jackknife rotation ≤ 25° when the jackknife is possible. Grade A
is unreachable with three stations by construction: the bar falls as
the station count rises so that dropping stations can never buy a
better grade (Triantafyllis et al. 2016 show two stations reaching
VR 0.9 with a condition number above 10). Only A/B are emailed.

**No coherent solution.** When fewer than three stations survive the
loop, the event is recorded in `events/NOSOL/` with
`"status": "no_coherent_solution"` — the full station ledger, no
mechanism — rather than a number fitted to noise; F-net likewise does
not publish below its station floor (Fukuyama et al. 1998). Such events
are not rows of `catalogue.csv`; their stage and reason are listed in
`not_published.csv` and printed by the human-review notebook.

## Outputs

- `events/NOSOL/<publicID>.json` — every event attempted without a
  coherent solution, in one directory, with its station map and
  all-station waveform figure beside it.
- `events/<publicID>_<date>_Mw…/solution.json` — everything: origin,
  depth search, preferred solution, both tensors, every station used or
  not with its SNR, own VR, time shift and reason, the jackknife, the
  quality block, the forward model, and `provenance.params` (the
  resolved parameter file the solution was made with).
- [`catalogue.csv`](catalogue.csv) — one row per event
  (column reference: [`events/CATALOGUE_README.md`](events/CATALOGUE_README.md)).
- [`not_published.csv`](not_published.csv) — every event with no
  solution: attempted but no coherent solution (stage, reason, best VR),
  or seen but not attempted (below the processing floor: too deep, too
  small, outside the box).
- [`station_ledger.jsonl`](station_ledger.jsonl) — one JSON line per
  event holding a per-station dictionary of numbers (distance, azimuth,
  SNR per component, amplitude ratio, own VR, time shift, used or drop
  class): the raw material for a station quality index.
  [`station_performance.csv`](station_performance.csv) is its
  per-station aggregate.
- Figures per event: stations + displacement field, depth sensitivity,
  mttime waveform fits, and the all-station waveform figure showing
  every candidate's record with its SNR, own VR and time shift.

## Students: the Human Review catalogue

Alongside the automated archive lives a **human-reviewed catalogue**
([events_human/](events_human/README.md)) built gradually by students
re-examining the automated solutions. Start with
**[human_review.ipynb](human_review.ipynb)**: installation, browse the
catalogue, watch a seismogram travel from raw counts to inversion-ready
displacement, run the inversion with full manual control, and submit
your reviewed solution by Pull Request. The human catalogue carries the
same columns as the automated one (including the full moment tensor)
plus the review fields. Read [docs/METHOD.md](docs/METHOD.md) and
[docs/REVIEW_LEARNINGS.md](docs/REVIEW_LEARNINGS.md) first.

## Local setup (macOS / Linux)

```sh
pixi install
pixi run test                 # anchor tests, including one that checks
                              # auto_tdmt.cfg and params.py agree
pixi run get-gfs              # one-time Green's function download (~1 GB)
pixi run process --event 2026p669681 --debug    # = python src/process_event.py
pixi run params               # print every resolved parameter
```

Outputs default to `~/work/proj_tdmt_NZ`; override with `AUTO_TDMT_OUTPUT`,
`AUTO_TDMT_EVENTS`, `AUTO_TDMT_GF`, and `AUTO_TDMT_CFG` for an alternative
parameter file. `--debug` keeps the staged Green's functions and SAC
data and writes stage-by-stage figures under `<event>/<band>/diagnostics/`.

## CI (GitHub Actions)

- `watch.yml` — twice an hour: poll, process every new event above the
  floor, largest first, inside a 300-minute budget (the rest wait for
  the next run), email passing solutions, commit results. This is not a
  real-time tool: GitHub queues scheduled workflows at low priority, so
  a new event is typically picked up within one to a few hours; the
  state file dedupes, so frequent polling never repeats work, and
  nothing this pipeline is for is lost by the delay (NISAR and
  Sentinel-1 passes are days apart). Run it by hand
  (`workflow_dispatch`) for anything urgent — that path skips the queue.
- `process.yml` — manual reprocess of one publicID.
- `publish.yml` — manual email of one processed event.
- `human_catalogue.yml` — rebuilds the human catalogue on PR merge.

The GF libraries are attached to the `gf-latest` release as
`gf_library.tar.zst` and cached in CI; CPS never runs in CI. Rebuild and
re-upload after any velocity-model change:

```sh
pixi run python src/greens.py --build
cd ~/work/proj_tdmt_NZ && cp -r gf_library gf_cache \
  && tar --zstd -cf gf_library.tar.zst gf_cache && rm -r gf_cache
gh release create gf-latest gf_library.tar.zst --notes "GF libraries"
```

Secrets required for email: `SMTP_HOST`, `SMTP_PORT`, `SMTP_USER`,
`SMTP_PASS`, `MAIL_FROM`, `MAIL_TO`.

## Velocity models

`models/nz_{north,south}_ristau2008.d` — Ristau (2008), SRL 79(3) Table 1,
doi:10.1785/gssrl.79.3.400 — the velocity models published for NZ
regional CMT analysis, so solutions can be compared with the published NZ
solutions ([GeoNet/data moment-tensor](https://github.com/GeoNet/data/tree/main/moment-tensor);
MT elements in both catalogues are in 1e20 dyne-cm).

## Validation

[validation/README.md](validation/README.md) compares the catalogue
against the Ristau NZ CMT, USGS and GCMT references on magnitude, depth
and mechanism rotation, and says which version of the selection its
numbers describe. `Selection` and `Code` columns in the catalogue make a
mixed-vintage catalogue self-describing; `tools/regrade_archive.py`
re-grades the archive under the current table without reprocessing.
The four CSV tables live at the repository root so they are the first
thing on the page. `pixi run validate` regenerates the comparison by
hand; `pixi run sync` regenerates it as part of publishing the archive.

## For the next maintainer

1. Change a number in `auto_tdmt.cfg`; `pixi run test` tells you if you
   broke the contract between the cfg and the code.
2. Run one task at a time with its notebook or its `--event` CLI.
3. After a year: `station_performance.csv` says which stations
   are consistently dropped and by which rule; that is the evidence for
   the next change. `docs/METHOD.md` §10 lists what to add next
   (distance-dependent bandwidth, MouseTrap, the condition number, a
   distance-regressed time-shift residual), each with its citation.
4. `git tag selection-v3` and `selection-v4` are the two earlier
   designs; `docs/REVIEW_LEARNINGS.md` records why they were replaced.

## Data sources

- GeoNet (Earth Sciences New Zealand) quake API + FDSN (NRT + archive):
  event detection, all original hypocentres, waveforms and station
  metadata. CC BY 3.0 NZ. Polling uses one request per run, gzip
  encoding and a descriptive User-Agent.
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
  doi:10.1785/0220110096
- Beyreuther, M., et al. (2010). ObsPy: a Python toolbox for seismology.
  *Seism. Res. Lett.* 81(3), 530-533. https://github.com/obspy/obspy
- Jolivet, R. — **okada4py**, https://github.com/jolivetr/okada4py
- Crameri, F. — Scientific colour maps, https://www.fabiocrameri.ch/colourmaps/

Operational systems reviewed for the station-selection and grading
rules, and what each one contributed (the full review, with DOIs and
what could not be verified, is
[docs/lit_review/operational_systems.md](docs/lit_review/operational_systems.md)):

| System | Reference | What we took from it |
|---|---|---|
| **Berkeley / SCSN** (the TDMT lineage this code descends from) | Dreger & Helmberger (1993) *JGR* 98, 8107-8125, doi:10.1029/93JB00023; Dreger (2003) *IASPEI Handbook* 81B; **Clinton, Hauksson & Solanki (2006)** *BSSA* 96(5), 1689-1705, doi:10.1785/0120050241 | the selection loop (invert, drop the worst-fitting stations, re-invert), the magnitude-scaled filter-band menu, the quality-B own-VR floor of 25 |
| **GeoNet / Ristau** (the reference catalogue) | Ristau (2008) *SRL* 79(3), 400-415, doi:10.1785/gssrl.79.3.400; Ristau (2013) *BSSA* 103(4), 2520-2533, doi:10.1785/0120120339 | the velocity models; the validation baseline (NS, DC, VR) |
| **Herrmann / CPS** (our Green's functions) | Herrmann (2013) *SRL* 84, 1081-1088; Herrmann, Benz & Ammon (2011) *BSSA* 101(6), 2609-2625, doi:10.1785/0120110095 | the 0.02-0.10 Hz band and the Mw 3.7 floor; per-station time shifts and VR as published diagnostics; the group-velocity record window |
| **USGS NEIC** | regional MT: Herrmann, Benz & Ammon (2011), above. Teleseismic W-phase and SynDepth — Duputel, Rivera, Kanamori & Hayes (2012) *GJI* 189(2), 1125-1147, doi:10.1111/j.1365-246X.2012.05419.x; Yeck et al. (2025) *SRL* 96(6), doi:10.1785/0220240372 — cited for design only, not as a method | the two-stage architecture (reject what needs no forward model, then prune by fit); the amplitude-ratio screen (used one-sided here); the 8 s time-shift cap |
| **INGV Italy** | Scognamiglio, Tinti & Michelini (2009) *BSSA* 99(4), 2223-2242, doi:10.1785/0120080104; INGV TDMT quality legend | the **grade table** (the VR bar falls as station count rises); 8 azimuth sectors; higher frequencies for small events |
| **ISOLA family** (Gisola, scisola, Bayesian ISOLA, BMKG) | Triantafyllis et al. (2022) *SRL* 93(2A), 957-966; Triantafyllis, Sokos, Ilias & Zahradník (2016) *SRL* 87(1), 157-163; Vackář et al. (2017) *GJI* 210(2), 693-705, doi:10.1093/gji/ggx158; Halauwet et al. (2024) *GJI* 239(2), 1000-1020, doi:10.1093/gji/ggae309; Zahradník & Sokos (2018) doi:10.1007/978-3-319-77359-9_1 | the magnitude-scaled distance window; the sector-balanced station cap of 16; the per-component RMS SNR with threshold 2; the ±30 km depth window; the warning that VR alone is gameable by dropping stations |
| **NIED F-net** (Japan) | Fukuyama, Ishida, Dreger & Kawai (1998) *Zisin* 51(1), 149-156, doi:10.4294/zisin1948.51.1_149; Kubo et al. (2002) *Tectonophysics* 356, 23-48 | depth searched within ±30 km of the hypocentre; distance weighting; the jackknife as the detector of a bad station; "no solution below the floor" |
| **gempa / GEOFON** | docs.gempa.de scautomt | the fixed per-station fit floor that pruning cannot breach below the minimum station count |
| **SCARDEC, SED, AutoBATS** | Vallée et al. (2011) *GJI* 184(1), 338-358; Bernardi et al. (2004) *GJI* 157(2), 703-716; Jian et al. (2018) *BSSA* 108, doi:10.1785/0120170231 | depth uncertainty as the range within 10% of the best misfit; best-SNR-per-azimuth-bin selection |
| **Mechanism comparison** | Townend et al. (2012) supplement; Walsh, Arnold & Townend (2009) *GJI*; Kagan (1991) | the minimum rotation angle used in validation and the jackknife |
| **Forward model** | Okada (1992) *BSSA* 82(2), 1018-1040; Wells & Coppersmith (1994) *BSSA* 84(4), 974-1002; Aki & Richards (1980) | the displacement forecast and its fault dimensions |

This workflow was compiled with Claude (Anthropic) assistance under the
direction of Danielle Lindsay; the science stands on the shoulders of the
authors above — cite them, not this repository, for the methods.
