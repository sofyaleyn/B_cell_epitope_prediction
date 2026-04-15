"""Stage 1 — run any or all predictors for a target.

Saves raw results to outputs/<target>/:
    bepipred_raw.csv
    discotope_raw.csv
    graphbepi_raw.csv

Usage:
    uv run python scripts/predict.py --target <TARGET> --predictors all
    uv run python scripts/predict.py --target <TARGET> --predictors graphbepi
    uv run python scripts/predict.py --target <TARGET> --predictors bepipred discotope
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import requests

from epitope_pipeline.config import resolve_target_config
from epitope_pipeline.predictors import bepipred, discotope, graphbepi

# Exceptions that indicate a web service is down, not a code bug.
_NETWORK_ERRORS = (
    requests.exceptions.ConnectionError,
    requests.exceptions.Timeout,
    requests.exceptions.HTTPError,
    TimeoutError,
    OSError,
)

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

    target_config = resolve_target_config(args.target, REPO_ROOT)

    to_run = list(PREDICTORS) if "all" in args.predictors else args.predictors

    for name in to_run:
        print(f"[{args.target}] Running {name}...")
        try:
            PREDICTORS[name](target_config)
            print(f"[{args.target}] {name} done.")
        except NotImplementedError:
            print(f"[{args.target}] {name} — not yet implemented, skipping.")
        except _NETWORK_ERRORS as e:
            print(f"\n[{args.target}] {name} — web service unavailable: {e}")
            print(f"[{args.target}] Pipeline stopped. Fix connectivity and re-run predict.py.")
            sys.exit(1)
        except Exception as e:
            print(f"[{args.target}] {name} — ERROR: {e}")
            raise


if __name__ == "__main__":
    main()
