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

## How it works

```
GeoNet quake API (poll, 10 min cron)
  -> processing floor (prelim M >= 4.0, NZ bbox)          run01_watch.py
  -> waveforms: GeoNet NRT FDSN, NZ broadbands <= 400 km   waveforms.py
     response removal -> ZRT -> bandpass -> 1 sps SAC
  -> Green's functions: precomputed CPS library            greens.py
     (Ristau 2008 North/South Island models, 10-500 km,
      depths 2-58 km, 10 Herrmann fundamental sources)
  -> mttime deviatoric inversion, depth search,            invert.py
     BSL-style filter-band menu, low-VR station rejection,
     preferred pick = max %DC within 5 VR points of max
  -> quality gates (stations, VR, azimuthal gap, depth)    invert.py
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
the public archive — and [`events/catalogue.csv`](events/catalogue.csv) is
regenerated from the archived solutions after every event: one row per
solution with origin, nodal planes, Mw/Mo, centroid depth, %DC/%CLVD, VR and
MT elements (1e20 dyne-cm, matching the GeoNet CMT catalogue conventions). Every `solution.json` records the velocity model, GF
version, package versions, stations used/dropped (with reasons), the full
depth-search table and the band search.

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

## Data sources

- GeoNet quake API + FDSN (NRT + archive), CC BY 3.0 NZ. Polling is polite:
  one request per 10-minute cron tick, gzip, descriptive User-Agent.
- NASA CMR for NISAR GSLC granule timing.
