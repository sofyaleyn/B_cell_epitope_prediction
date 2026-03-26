"""GraphBepi local runner and output parser.

Runs graphbepi/test.py on a single-chain PDB (CPU only, M4 compatible).

Pipeline entry point: run(target_config)
    Reads pdb_dir and chains from target_config, runs GraphBepi on each
    chain, saves graphbepi_raw.csv to out_dir.

Output columns: res_id (int), chain (str), score (float), is_epitope (bool)
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd


THRESHOLD = 0.5
GRAPHBEPI_DIR = Path(__file__).parents[4] / "graphbepi"


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def run(target_config: dict) -> None:
    """Run GraphBepi for all chains and save to out_dir/graphbepi_raw.csv.

    Reads target_config keys:
        pdb_dir (Path): folder containing per-chain .pdb files
        chains  (list): chain identifiers to process
        out_dir (Path): destination for graphbepi_raw.csv
    """
    pdb_dir = Path(target_config["pdb_dir"])
    chains  = target_config.get("chains", ["A"])
    out_dir = Path(target_config["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    cache_path = out_dir / "graphbepi_raw.csv"

    if cache_path.exists():
        return

    pdb_files = list(pdb_dir.glob("*.pdb"))
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in {pdb_dir}")

    # Run on first chain's PDB; extend when multi-chain support is needed
    chain = chains[0]
    _run_graphbepi(pdb_files[0], chain, cache_path)


# ---------------------------------------------------------------------------
# Local runner
# ---------------------------------------------------------------------------

def _run_graphbepi(pdb_path: Path, chain: str, cache_path: Path) -> None:
    """Execute graphbepi/test.py and write output to cache_path."""
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(GRAPHBEPI_DIR / "test.py"),
        "--pdb", str(pdb_path),
        "--chain", chain,
        "--output", str(cache_path),
        # Force CPU — no --gpu flag; GraphBepi defaults to CPU when omitted
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=GRAPHBEPI_DIR)
    if result.returncode != 0:
        raise RuntimeError(f"GraphBepi failed:\n{result.stderr}")


# ---------------------------------------------------------------------------
# Output parser (column mapping TBD once real output is confirmed)
# ---------------------------------------------------------------------------

def _parse(csv_path: Path, chain: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    # TODO: map GraphBepi column names to standard schema once output confirmed
    raise NotImplementedError("GraphBepi output parser not yet implemented")
