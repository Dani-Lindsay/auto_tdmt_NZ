# events_human — the human-reviewed catalogue

> **Made by Claude. Lightly tested.** The workflow that populates this
> directory (`human_review.ipynb` in the repository root) was generated
> as a starting point for students; one full end-to-end execution has
> succeeded, but expect rough edges.

This directory mirrors the automated archive `events/`, but every entry
here has been **reviewed by a person**: a student or researcher stepped
through the automated solution, re-ran the inversion, tweaked station
selection / filter band / depth range where justified, and recorded
their reasoning. It grows gradually as events receive human eyes.

- `catalogue_human.csv` — the SHARED master catalogue that many
  reviewers append to, one row per reviewed event. Its columns
  mirror the automated `events/catalogue.csv` (so the two can be
  compared directly) plus reviewer columns: the reviewer, the decision (`accept_automated` / `revised` /
  `reject`), and what changed relative to the automated answer.
- `<eventID>_Mw..._..._.../` — the reviewed solution: `solution.json`
  (with a `human_review` block and the automated reference), the
  reviewer's waveform-fit figure, and copies of the automated figures
  (`auto_*.jpg`) for side-by-side comparison.

**Rules**
1. The automated archive `events/` is never modified by review work.
   The two catalogues coexist: `events/` is what the machine said,
   `events_human/` is what a person concluded.
2. Contributions arrive by Pull Request touching **only** this
   directory (see the notebook, step 8). PRs touching anything else are
   closed unmerged.
3. Every solution names its reviewer and states its reasoning. "The
   automated solution is correct" and "no defensible solution exists"
   are both valid, valuable review outcomes.
4. Merge conflicts in `catalogue_human.csv` (two reviews merged close
   together) are normal: `git pull --rebase`, re-run the notebook's
   staging cell (it re-appends your row idempotently), commit again.

**Getting started**: open `human_review.ipynb` in the repository root
and follow it from installation to your first Pull Request. Read
`docs/METHOD.md` and `docs/REVIEW_LEARNINGS.md` first.

Data: GeoNet (CC BY 3.0 NZ). Method credits as in the repository README
— cite the original authors (Ristau 2008; Dreger; Chiang/mttime;
Herrmann/CPS), not this repository.
