"""Sync the local event archive into the repo and rebuild the public
products (catalogue, station ledger, overview map, validation), then
commit and push. This is the DELIBERATE publish step: test sweeps
(run02/run04) write only to the local archive, and nothing reaches
GitHub until this is run — so experimental rule iterations do not churn
the public archive or bloat git history with throwaway figures.

    pixi run python run06_sync_repo.py --message "why this state is worth publishing"
    pixi run python run06_sync_repo.py --no-push   # local commit only
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import warnings

import config

GIT_AUTHOR = ["-c", "user.name=dani-lindsay",
              "-c", "user.email=danielle.lindsay@earthsciences.nz"]


def sync_events() -> int:
    repo_events = config.REPO_DIR / "events"
    assert config.EVENTS_DIR.resolve() != repo_events.resolve(), \
        "EVENTS_DIR already is the repo events/ (CI layout) — nothing to sync"
    n = 0
    for p in sorted(config.EVENTS_DIR.glob("*/solution.json")):
        pid = json.loads(p.read_text())["event"]["public_id"]
        for old in repo_events.glob(f"{pid}*"):
            if old.is_dir():
                shutil.rmtree(old)
        dst = repo_events / p.parent.name
        dst.mkdir()
        for f in [p, p.parent / "draft_email.txt", *p.parent.glob("*_*.jpg")]:
            if f.exists():
                shutil.copy(f, dst / f.name)
        n += 1
    return n


def main(message: str, push: bool) -> None:
    n = sync_events()
    print(f"synced {n} events into repo events/")

    import catalogue
    import figure
    import station_performance
    warnings.filterwarnings("ignore")
    repo_events = config.REPO_DIR / "events"
    catalogue.build_catalogue(repo_events)
    station_performance.build_station_performance(repo_events)
    figure.make_overview_map(repo_events, repo_events / "solutions_map.jpg")
    (config.REPO_DIR / "validation").mkdir(exist_ok=True)
    report = config.REPO_DIR / "validation" / "validation_report.txt"
    with open(report, "w") as f:
        subprocess.run(
            ["pixi", "run", "python", "run05_validate.py"],
            cwd=config.REPO_DIR, stdout=f, stderr=subprocess.STDOUT,
            check=False)
    print(f"validation report: {report}")

    def git(*args, check=True):
        return subprocess.run(["git", *args], cwd=config.REPO_DIR,
                              check=check)

    git("add", "-A")
    r = subprocess.run(
        ["git", *GIT_AUTHOR, "commit", "-qm",
         f"{message}\n\nCo-Authored-By: Claude Fable 5 "
         "<noreply@anthropic.com>"],
        cwd=config.REPO_DIR, check=False)
    if r.returncode != 0:
        print("nothing to commit")
    elif push:
        git("pull", "-q", "--rebase", "origin", "main")
        git("push", "-q", "origin", "main")
        print("pushed")
    else:
        print("committed locally (--no-push)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--message", default="sync event archive to repo",
                    help="commit message (say why this state is publishable)")
    ap.add_argument("--no-push", action="store_true",
                    help="commit locally without pushing")
    args = ap.parse_args()
    main(args.message, push=not args.no_push)
