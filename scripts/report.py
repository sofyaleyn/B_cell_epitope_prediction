"""Stage 3 — generate plots, tables, and execute the analysis notebook.

Reads outputs/<target>/combined_scores.csv and writes:
    bepipred_profile.png
    discotope_profile.png
    B_cell_epitope_combined.png
    epitope_table.csv
    summary_report.html
    analysis.ipynb          (executed notebook)

Usage:
    uv run python scripts/report.py --target <TARGET>
    uv run python scripts/report.py --target <TARGET> --skip-notebook
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import pandas as pd

from epitope_pipeline.reporting.plots import plot_epitope_scores
from epitope_pipeline.reporting.tables import export_epitope_table
from epitope_pipeline.reporting.report import export_html
from epitope_pipeline.reporting.notebook import execute_template


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate epitope prediction report")
    parser.add_argument("--target", required=True, metavar="NAME",
                        help="Target name matching a folder under outputs/")
    parser.add_argument("--skip-notebook", action="store_true",
                        help="Skip notebook execution (plots and tables only)")
    args = parser.parse_args()

    out_dir = REPO_ROOT / "outputs" / args.target
    csv_path = out_dir / "combined_scores.csv"
    if not csv_path.exists():
        raise SystemExit(
            f"combined_scores.csv not found in {out_dir}\n"
            f"Run stage 2 first:  uv run python scripts/combine.py --target {args.target}"
        )

    df = pd.read_csv(csv_path)
    print(f"[{args.target}] Generating plots...")
    plot_epitope_scores(df, args.target, out_dir / "B_cell_epitope_combined.png")

    print(f"[{args.target}] Exporting tables...")
    export_epitope_table(df, out_dir / "epitope_table.csv")

    print(f"[{args.target}] Exporting HTML report...")
    export_html(df, args.target, out_dir / "summary_report.html")

    if not args.skip_notebook:
        print(f"[{args.target}] Executing analysis notebook...")
        execute_template(
            template=REPO_ROOT / "notebooks" / "analysis.ipynb",
            output=out_dir / "analysis.ipynb",
            target=args.target,
        )

    print(f"[{args.target}] Done → {out_dir}")


if __name__ == "__main__":
    main()
