"""CLI to run the full epitope prediction pipeline for one or more targets.

Usage:
    uv run python scripts/run_pipeline.py --target her2
    uv run python scripts/run_pipeline.py --target her2 --target egfr
    uv run python scripts/run_pipeline.py          # runs all targets in config
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).parents[1]
CONFIG_PATH = REPO_ROOT / "config" / "targets.yaml"
OUTPUTS_DIR = REPO_ROOT / "outputs"


def main() -> None:
    parser = argparse.ArgumentParser(description="Run B-cell epitope prediction pipeline")
    parser.add_argument(
        "--target",
        action="append",
        dest="targets",
        metavar="NAME",
        help="Target name from config/targets.yaml (repeat for multiple targets; default: all)",
    )
    args = parser.parse_args()

    config = yaml.safe_load(CONFIG_PATH.read_text())
    all_targets = {t["name"]: t for t in config["targets"]}

    selected = args.targets if args.targets else list(all_targets.keys())
    for name in selected:
        if name not in all_targets:
            raise SystemExit(f"Unknown target '{name}'. Available: {list(all_targets)}")

    for name in selected:
        _run_target(all_targets[name])


def _run_target(target: dict) -> None:
    from epitope_pipeline.predictors import bepipred, discotope, graphbepi
    from epitope_pipeline.integration.combine import combine
    from epitope_pipeline.reporting.plots import plot_epitope_scores
    from epitope_pipeline.reporting.report import export_html
    from Bio import SeqIO
    from Bio.PDB import PDBParser

    name = target["name"]
    out_dir = OUTPUTS_DIR / name
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{name}] Running predictors...")

    results = {}
    for chain in target["chains"]:
        pdb_path = REPO_ROOT / target["pdb"]
        chain_pdb_path = REPO_ROOT / target["chain_pdbs"][chain]

        sequence = _extract_sequence(chain_pdb_path, chain)

        results["bepipred"] = bepipred.run(
            sequence=sequence,
            chain=chain,
            cache_path=out_dir / "bepipred_raw.json",
        )
        results["discotope"] = discotope.run(
            pdb_path=chain_pdb_path,
            chain=chain,
            cache_path=out_dir / "discotope_raw.json",
        )
        results["graphbepi"] = graphbepi.run(
            pdb_path=chain_pdb_path,
            chain=chain,
            cache_path=out_dir / "graphbepi_raw.csv",
        )

    print(f"[{name}] Combining results...")
    combined = combine(results)
    combined.to_csv(out_dir / "combined_scores.csv", index=False)

    print(f"[{name}] Generating report...")
    plot_epitope_scores(combined, name, out_dir / "epitope_plot.png")
    export_html(combined, name, out_dir / "summary_report.html")

    print(f"[{name}] Done → {out_dir}")


def _extract_sequence(pdb_path: Path, chain: str) -> str:
    from Bio.PDB import PDBParser
    from Bio.SeqUtils import seq1

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    residues = [
        r for r in structure[0][chain].get_residues()
        if r.get_id()[0] == " "  # exclude HETATM
    ]
    return "".join(seq1(r.get_resname()) for r in residues)


if __name__ == "__main__":
    main()
