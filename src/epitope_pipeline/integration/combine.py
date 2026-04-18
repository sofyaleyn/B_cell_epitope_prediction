"""Load and combine B-cell epitope predictions for a target.

Main entry point: load_and_combine(target_dir, out_dir) — auto-discovers all
predictor result files, aligns predictions to the full antigen sequence,
and returns a single merged DataFrame.

Expected target_dir layout:
    <target_dir>/
    ├── *.fasta                   # antigen sequence (first ERCC1-like record used)
    ├── *bepipred*/
    │   └── raw_output.csv
    ├── *discotope*/
    │   └── <pdb_stem>_discotope3.csv   (one per structure)
    └── pdb/
        └── <pdb_stem>.pdb

GraphBepi results are read from out_dir/graphbepi_raw.csv (written by predict.py).

Output columns:
    res_id (int), residue (str), bepipred_score (float),
    discotope_score (float), average_rsa (float), graphbepi_score (float),
    is_epitope_bepipred (bool), is_epitope_discotope (bool),
    is_epitope_graphbepi (bool), is_epitope_AND (bool)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from epitope_pipeline.predictors import bepipred as _bepipred
from epitope_pipeline.predictors import discotope as _discotope
from epitope_pipeline.predictors.graphbepi import THRESHOLD as GRAPHBEPI_THRESHOLD

BEPIPRED_THRESHOLD  = 0.50   # BepiPred-3.0 recommended default
DISCOTOPE_THRESHOLD = 0.70   # DiscoTope-3.0 calibrated score


def load_and_combine(target_dir: Path, out_dir: Path | None = None) -> pd.DataFrame:
    """Discover, parse, and merge all predictor results for a target.

    Args:
        target_dir: root directory for one target (e.g. data/ERCC1/).
        out_dir:    outputs directory for the target (e.g. outputs/ERCC1/).
                    If provided and graphbepi_raw.csv exists there, GraphBepi
                    scores are merged in. If None, GraphBepi is skipped.

    Returns:
        DataFrame (one row per residue in the antigen sequence).
        Scores are 0.0 where a predictor had no coverage for that residue.
    """
    target_dir = Path(target_dir)

    sequence = _load_sequence(target_dir)
    base = pd.DataFrame({
        "res_id":  range(1, len(sequence) + 1),
        "residue": list(sequence),
    })

    pdb_dir = target_dir / "pdb"

    # BepiPred — use cached raw CSV from out_dir if available (written by predict.py),
    # otherwise fall back to parsing from data/ directory.
    bp_path = Path(out_dir) / "bepipred_raw.csv" if out_dir else None
    if bp_path and bp_path.exists():
        bp = pd.read_csv(bp_path)[["res_id", "bepipred_score"]]
    else:
        bepipred_dir = _find_subdir(target_dir, "*bepipred*")
        bp = _bepipred.parse_results_dir(bepipred_dir)[["res_id", "bepipred_score"]]

    # DiscoTope — same: prefer cached raw CSV so combine.py works regardless of
    # where the raw DiscoTope results live (data/ vs outputs/).
    dt_path = Path(out_dir) / "discotope_raw.csv" if out_dir else None
    if dt_path and dt_path.exists():
        _dt_all = pd.read_csv(dt_path)
        _dt_cols = ["res_id", "discotope_score", "average_rsa"]
        if "discotope_score_std" in _dt_all.columns:
            _dt_cols.append("discotope_score_std")
        dt = _dt_all[_dt_cols]
    else:
        discotope_dir = _find_subdir(target_dir, "*discotope*")
        dt = _discotope.parse_results_dir(discotope_dir, pdb_dir)[
            ["res_id", "discotope_score", "discotope_score_std", "average_rsa"]
        ]

    df = base.merge(dt, on="res_id", how="left")
    df = df.merge(bp, on="res_id", how="left")

    df["bepipred_score"]  = df["bepipred_score"].fillna(0.0)
    df["discotope_score"] = df["discotope_score"].fillna(0.0)

    # GraphBepi — optional, read from out_dir if available
    gb_path = Path(out_dir) / "graphbepi_raw.csv" if out_dir else None
    if gb_path and gb_path.exists():
        _gb_all = pd.read_csv(gb_path)
        _gb_cols = {"score": "graphbepi_score"}
        if "score_std" in _gb_all.columns:
            _gb_cols["score_std"] = "graphbepi_score_std"
        gb = _gb_all.rename(columns=_gb_cols)[list(_gb_cols.values()) + ["res_id"]]
        df = df.merge(gb, on="res_id", how="left")
    else:
        df["graphbepi_score"] = float("nan")
    df["graphbepi_score"] = df["graphbepi_score"].fillna(0.0)

    df["is_epitope_bepipred"]  = df["bepipred_score"]  >= BEPIPRED_THRESHOLD
    df["is_epitope_discotope"] = df["discotope_score"] >= DISCOTOPE_THRESHOLD
    df["is_epitope_graphbepi"] = df["graphbepi_score"] >= GRAPHBEPI_THRESHOLD
    df["is_epitope_AND"] = (
        df["is_epitope_bepipred"]
        | df["is_epitope_discotope"]
        | df["is_epitope_graphbepi"]
    )

    dt_std = ["discotope_score_std"] if "discotope_score_std" in df.columns else []
    gb_std = ["graphbepi_score_std"] if "graphbepi_score_std" in df.columns else []
    return df[["res_id", "residue",
               "bepipred_score",
               "discotope_score", *dt_std, "average_rsa",
               "graphbepi_score", *gb_std,
               "is_epitope_bepipred", "is_epitope_discotope", "is_epitope_graphbepi",
               "is_epitope_AND"]]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_sequence(target_dir: Path) -> str:
    """Return the antigen sequence from FASTA files in target_dir.

    Searches all .fasta/.fa files for a record whose header contains the
    target directory name (case-insensitive). If none match across any file,
    falls back to the first record of the first file.
    """
    fasta_files = sorted(
        list(target_dir.glob("*.fasta")) + list(target_dir.glob("*.fa"))
    )
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA file found in {target_dir}")

    target_name = target_dir.name.upper()

    # Search all files for a record whose header contains the target name.
    for fasta_path in fasta_files:
        sequence = ""
        capture = False
        with open(fasta_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if capture:
                        break   # finished the matching record
                    if target_name in line.upper():
                        capture = True
                elif capture:
                    sequence += line
        if sequence:
            return sequence.upper()

    # Fall back: first record of first file.
    sequence = ""
    with open(fasta_files[0]) as f:
        capture = False
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if capture:
                    break
                capture = True
            elif capture:
                sequence += line

    if not sequence:
        raise ValueError(f"Could not extract sequence from {fasta_files[0]}")
    return sequence.upper()


def _find_subdir(target_dir: Path, pattern: str) -> Path:
    matches = [p for p in target_dir.iterdir() if p.is_dir() and p.match(pattern)]
    if not matches:
        raise FileNotFoundError(
            f"No subdirectory matching '{pattern}' found in {target_dir}"
        )
    return matches[0]
