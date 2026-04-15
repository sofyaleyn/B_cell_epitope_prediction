"""Run the full epitope prediction pipeline for one or more targets.

Executes all three stages in sequence:
    1. predict   — run BepiPred, DiscoTope, GraphBepi
    2. combine   — merge results into combined_scores.csv
    3. report    — plots, tables, HTML, notebook

Usage:
    uv run python scripts/run_pipeline.py --target ERCC1
    uv run python scripts/run_pipeline.py --target ERCC1 --target E3
    uv run python scripts/run_pipeline.py                    # all targets in config
    uv run python scripts/run_pipeline.py --skip-notebook
    uv run python scripts/run_pipeline.py --predictors graphbepi
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import requests
import yaml

from epitope_pipeline.config import resolve_target_config
from epitope_pipeline.integration.combine import load_and_combine
from epitope_pipeline.reporting.plots import plot_epitope_scores
from epitope_pipeline.reporting.tables import export_epitope_table
from epitope_pipeline.reporting.report import export_html
from epitope_pipeline.reporting.notebook import execute_template
from epitope_pipeline.predictors import bepipred, discotope, graphbepi

PREDICTORS = {
    "bepipred":  bepipred.run,
    "discotope": discotope.run,
    "graphbepi": graphbepi.run,
}

_NETWORK_ERRORS = (
    requests.exceptions.ConnectionError,
    requests.exceptions.Timeout,
    requests.exceptions.HTTPError,
    TimeoutError,
    OSError,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run full B-cell epitope prediction pipeline")
    parser.add_argument(
        "--target", action="append", dest="targets", metavar="NAME",
        help="Target name from config/targets.yaml (repeat for multiple; default: all)",
    )
    parser.add_argument(
        "--predictors", nargs="+", choices=[*PREDICTORS, "all"], default=["all"],
        metavar="NAME",
        help=f"Predictors to run: all | {' | '.join(PREDICTORS)}",
    )
    parser.add_argument(
        "--skip-notebook", action="store_true",
        help="Skip notebook execution",
    )
    args = parser.parse_args()

    # Load yaml only to enumerate all targets for the default "run all" case
    config = yaml.safe_load((REPO_ROOT / "config" / "targets.yaml").read_text())
    all_targets = list(config["targets"].keys())
    selected = args.targets if args.targets else all_targets

    to_run = list(PREDICTORS) if "all" in args.predictors else args.predictors

    for name in selected:
        _run_target(name, to_run, args.skip_notebook)


def _run_target(name: str, predictors: list[str], skip_notebook: bool) -> None:
    cfg = resolve_target_config(name, REPO_ROOT)
    target_dir = REPO_ROOT / "data" / name
    out_dir    = cfg["out_dir"]

    print(f"[{name}] Stage 1 — predictors: {predictors}")
    for predictor in predictors:
        print(f"[{name}]   running {predictor}...")
        try:
            PREDICTORS[predictor](cfg)
        except NotImplementedError:
            print(f"[{name}]   {predictor} — not yet implemented, skipping.")
        except _NETWORK_ERRORS as e:
            print(f"\n[{name}]   {predictor} — web service unavailable: {e}")
            print(f"[{name}]   Pipeline stopped. Stages 2–3 skipped. Fix connectivity and re-run.")
            sys.exit(1)
        except Exception as e:
            print(f"[{name}]   {predictor} — ERROR: {e}")
            raise

    print(f"[{name}] Stage 2 — combining results...")
    df = load_and_combine(target_dir, out_dir=out_dir)
    df.to_csv(out_dir / "combined_scores.csv", index=False)
    for col in ["is_epitope_bepipred", "is_epitope_discotope", "is_epitope_graphbepi", "is_epitope_AND"]:
        if col in df.columns:
            n = df[col].sum()
            print(f"  {col}: {n} / {len(df)} ({100*n/len(df):.1f}%)")

    print(f"[{name}] Stage 3 — generating report...")
    plot_epitope_scores(df, name, out_dir / "B_cell_epitope_combined.png")
    export_epitope_table(df, out_dir / "epitope_table.csv")
    export_html(df, name, out_dir / "summary_report.html")
    if not skip_notebook:
        execute_template(
            template=REPO_ROOT / "notebooks" / "analysis.ipynb",
            output=out_dir / "analysis.ipynb",
            target=name,
        )

    print(f"[{name}] Done → {out_dir}")


if __name__ == "__main__":
    main()
