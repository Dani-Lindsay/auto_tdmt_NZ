# Literature review: station selection and quality control in operational moment-tensor systems

Before settling the station-selection and grading rules of this
pipeline, we reviewed how established operational regional
moment-tensor services make the same decisions. The aim was to adopt
published practice wherever it exists, cite it, and reserve invention
for the few places where the literature is silent.

The review asked six questions of each system:

1. What is excluded **before** the inversion, and what is judged by fit afterwards?
2. How is signal quality measured (SNR definition; per station or per component)?
3. How are per-station time shifts bounded or vetted?
4. What distance range and azimuthal coverage are required?
5. How is source depth chosen, and is the search tied to the hypocentre?
6. What quality metrics gate publication?

Files:

- [`operational_systems.md`](operational_systems.md) — the review across
  twelve systems (Berkeley/SCSN, GeoNet/Ristau, Herrmann/CPS, USGS
  W-phase, INGV, ISOLA/Gisola/scisola/BMKG, NIED F-net, GFZ/gempa, SED,
  SCARDEC, AutoBATS), a comparison table, and for each recommendation
  its status in the current pipeline (adopted, adapted, or deferred).
- [`neic.md`](neic.md) — a closer look at the USGS NEIC stack, including
  the operational SynDepth parameters.

Conventions: **[verified]** means read in the primary source, code or
configuration file; **[secondary]** means reported only via a citing
paper; **[not verified]** means it could not be confirmed and is
flagged rather than assumed. No number or DOI in these documents is
inferred; where a source could not be obtained, that is stated.

The rules that resulted are listed with their sources in the README
("The rules and where they come from") and in `auto_tdmt.cfg`, where
each parameter's default carries its citation.
