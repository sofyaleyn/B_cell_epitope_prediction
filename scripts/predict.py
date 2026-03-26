"""Stage 1 — run any or all predictors for a target.

Saves raw results to outputs/<target>/:
    bepipred_raw.csv
    discotope_raw.csv
    graphbepi_raw.csv

Usage:
    uv run python scripts/predict.py --target ERCC1 --predictors all
    uv run python scripts/predict.py --target ERCC1 --predictors graphbepi
    uv run python scripts/predict.py --target ERCC1 --predictors bepipred discotope
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import yaml

from epitope_pipeline.predictors import bepipred, discotope, graphbepi

PREDICTORS = {
    "bepipred":  bepipred.run,
    "discotope": discotope.run,
    "graphbepi": graphbepi.run,
}


def main() -> None:
    parser = argparse.ArgumentParser(description="Run epitope predictors for a target")
    parser.add_argument("--target", required=True, metavar="NAME",
                        help="Target name from config/targets.yaml")
    parser.add_argument("--predictors", nargs="+",
                        choices=[*PREDICTORS, "all"], default=["all"],
                        metavar="NAME",
                        help=f"Predictors to run: all | {' | '.join(PREDICTORS)}")
    args = parser.parse_args()

    config = yaml.safe_load((REPO_ROOT / "config" / "targets.yaml").read_text())
    if args.target not in config["targets"]:
        raise SystemExit(
            f"Unknown target '{args.target}'. Available: {list(config['targets'])}"
        )

    target_config = config["targets"][args.target]

    # Resolve relative paths to absolute and inject runtime dirs
    for key in ("bepipred_dir", "discotope_dir", "pdb_dir"):
        if key in target_config:
            target_config[key] = REPO_ROOT / target_config[key]

    out_dir = REPO_ROOT / "outputs" / args.target
    out_dir.mkdir(parents=True, exist_ok=True)
    target_config["out_dir"] = out_dir

    to_run = list(PREDICTORS) if "all" in args.predictors else args.predictors

    for name in to_run:
        print(f"[{args.target}] Running {name}...")
        try:
            PREDICTORS[name](target_config)
            print(f"[{args.target}] {name} done.")
        except NotImplementedError:
            print(f"[{args.target}] {name} — not yet implemented, skipping.")
        except Exception as e:
            print(f"[{args.target}] {name} — ERROR: {e}")
            raise


if __name__ == "__main__":
    main()
