"""Stage 3 — generate plots, tables, and execute the analysis notebook.

Reads outputs/<target>/combined_scores.csv and writes:
    bepipred_profile.png
    discotope_profile.png
    B_cell_epitope_combined.png
    B_cell_epitope_combined_with_interface.png  (only with --include-interface)
    epitope_table.csv
    summary_report.html
    analysis.ipynb          (executed notebook)

Usage:
    uv run python scripts/report.py --target <TARGET>
    uv run python scripts/report.py --target <TARGET> --skip-notebook
    uv run python scripts/report.py --target <TARGET> --include-interface
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import pandas as pd
import yaml

from epitope_pipeline.reporting.plots import plot_epitope_scores, plot_per_structure
from epitope_pipeline.reporting.tables import export_epitope_table
from epitope_pipeline.reporting.report import export_html
from epitope_pipeline.reporting.notebook import execute_template


def _load_interface(target: str) -> tuple[pd.DataFrame | None, Path | None]:
    """Load antigen interface CSV for target if configured. Returns (df, path) or (None, None)."""
    config_path = REPO_ROOT / "config" / "targets.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)

    target_cfg = config.get("targets", {}).get(target, {})
    rel_path = target_cfg.get("antigen_interface")
    if not rel_path:
        return None, None

    interface_path = REPO_ROOT / rel_path
    if not interface_path.exists():
        print(f"[{target}] WARNING: antigen_interface path not found: {interface_path} — skipping interface panel")
        return None, None

    df = pd.read_csv(interface_path).rename(columns={"position": "res_id", "amino_acid": "residue"})
    return df, interface_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate epitope prediction report")
    parser.add_argument("--target", required=True, metavar="NAME",
                        help="Target name matching a folder under outputs/")
    parser.add_argument("--skip-notebook", action="store_true",
                        help="Skip notebook execution (plots and tables only)")
    parser.add_argument("--include-interface", action="store_true",
                        help="Overlay antigen-antibody interface contacts on the combined plot "
                             "(requires antigen_interface field in config/targets.yaml)")
    args = parser.parse_args()

    out_dir = REPO_ROOT / "outputs" / args.target
    csv_path = out_dir / "combined_scores.csv"
    if not csv_path.exists():
        raise SystemExit(
            f"combined_scores.csv not found in {out_dir}\n"
            f"Run stage 2 first:  uv run python scripts/combine.py --target {args.target}"
        )

    df = pd.read_csv(csv_path)

    # Resolve interface data if requested
    interface_df = None
    interface_path = None
    if args.include_interface:
        interface_df, interface_path = _load_interface(args.target)
        if interface_df is None:
            print(f"[{args.target}] No interface data available — generating plain combined plot only")

    print(f"[{args.target}] Generating plots...")
    combined_plot_path = out_dir / "B_cell_epitope_combined.png"
    # Always generate the plain combined plot
    plot_epitope_scores(df, args.target, combined_plot_path)
    # Generate interface-annotated plot if interface data available
    if interface_df is not None:
        plot_epitope_scores(df, args.target, combined_plot_path, interface_df=interface_df)

    # Per-structure DiscoTope plot (auto-discovers *Discotope* dir inside out_dir)
    discotope_per_struct_path = None
    discotope_multi_dirs = sorted(d for d in out_dir.iterdir()
                                  if d.is_dir() and "discotope" in d.name.lower())
    if discotope_multi_dirs:
        dt_csvs = sorted(discotope_multi_dirs[0].glob("*.csv"))
        discotope_per_struct_path = out_dir / "discotope_per_structure.png"
        if not plot_per_structure(
            csv_paths=dt_csvs,
            target_name=args.target,
            output_path=discotope_per_struct_path,
            score_col="calibrated_score",
            threshold=0.5,
            ylabel="Calibrated Score",
            title_prefix="DiscoTope 3.0: Per-Structure Predictions",
        ):
            discotope_per_struct_path = None

    # Per-structure GraphBepi plot (graphbepi_*.csv files in out_dir, excluding raw)
    graphbepi_per_struct_path = None
    gb_csvs = [p for p in out_dir.glob("graphbepi_*.csv") if p.name != "graphbepi_raw.csv"]
    if gb_csvs:
        graphbepi_per_struct_path = out_dir / "graphbepi_per_structure.png"
        if not plot_per_structure(
            csv_paths=gb_csvs,
            target_name=args.target,
            output_path=graphbepi_per_struct_path,
            score_col="score",
            threshold=0.1763,
            ylabel="Score",
            title_prefix="GraphBepi: Per-Structure Predictions",
        ):
            graphbepi_per_struct_path = None

    print(f"[{args.target}] Exporting tables...")
    export_epitope_table(df, out_dir / "epitope_table.csv")

    print(f"[{args.target}] Exporting HTML report...")
    export_html(
        df, args.target, out_dir / "summary_report.html",
        discotope_per_struct_img=discotope_per_struct_path,
        graphbepi_per_struct_img=graphbepi_per_struct_path,
    )

    if not args.skip_notebook:
        print(f"[{args.target}] Executing analysis notebook...")
        execute_template(
            template=REPO_ROOT / "notebooks" / "analysis.ipynb",
            output=out_dir / "analysis.ipynb",
            target=args.target,
            INTERFACE_CSV=str(interface_path) if interface_path is not None else "",
        )

    print(f"[{args.target}] Done → {out_dir}")


if __name__ == "__main__":
    main()
