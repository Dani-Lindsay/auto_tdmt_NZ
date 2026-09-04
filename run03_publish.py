"""Email one processed event's solution to the distribution list.

    pixi run python run03_publish.py --event 2026p660242
    pixi run python run03_publish.py --event 2026p660242 --force  # skip gates
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import config
import publish


def publish_event(public_id: str, state: dict | None = None,
                  force: bool = False) -> None:
    event_dir = config.find_event_dir(public_id)
    assert event_dir is not None, f"no archive for {public_id} — run run02 first"
    sol_path = event_dir / "solution.json"
    assert sol_path.exists(), f"no solution for {public_id} — run run02 first"
    solution = json.loads(sol_path.read_text())

    decision = solution["publish_decision"]
    if not decision["publish"] and not force:
        raise AssertionError(
            f"{public_id} did not pass publication gates: "
            f"{decision['reasons']} (use --force to override)"
        )

    draft = (event_dir / "draft_email.txt").read_text()
    subject, _, body = draft.partition("\n\n")

    attachments = sorted(event_dir.glob(f"{public_id}_waveform_fits_*.jpg"))
    attachments += [
        event_dir / f"{public_id}_stations_displacement_field.jpg",
        event_dir / f"{public_id}_depth_sensitivity.jpg",
    ]

    publish.send_email(subject, body.rstrip("\n"), attachments)

    own_state = state is None
    if own_state:
        from run01_watch import load_state
        state = load_state()
    assert config.is_solved(solution), \
        f"{public_id} has no coherent solution — nothing to publish"
    state["published"].append({
        "public_id": public_id,
        "mw": solution["preferred"]["mw"],
        "latitude": solution["event"]["latitude"],
        "longitude": solution["event"]["longitude"],
        "published_utc": datetime.now(timezone.utc).isoformat(),
    })
    if own_state:
        from run01_watch import save_state
        save_state(state)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--event", required=True)
    ap.add_argument("--force", action="store_true",
                    help="send even if publication gates failed")
    args = ap.parse_args()
    publish_event(args.event, force=args.force)
