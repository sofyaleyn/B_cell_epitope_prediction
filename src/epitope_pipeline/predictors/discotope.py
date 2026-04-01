"""DiscoTope web server scraper and local result parser.

Web API (run): submits a single-chain PDB to the DiscoTope-3.0 web server,
polls for results, caches raw JSON to outputs/<target>/discotope_raw.json.

Local parser (parse_results_dir): auto-discovers CSVs in a *discotope*
results directory whose stems match PDB files in the companion pdb/ dir,
averages calibrated_score across structures, and returns a standardised
DataFrame aligned to the full target sequence.

Output columns: res_id (int), residue (str), discotope_score (float),
                average_rsa (float), n_structures (int)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import requests


THRESHOLD = 0.5


# ---------------------------------------------------------------------------
# Local result parser
# ---------------------------------------------------------------------------

def parse_results_dir(results_dir: Path, pdb_dir: Path) -> pd.DataFrame:
    """Parse DiscoTope-3.0 CSVs for all structures and average across them.

    Discovers CSVs in results_dir whose stem matches a PDB file stem in
    pdb_dir (e.g. 2A1I_A_discotope3.csv matches 2A1I_A_discotope3.pdb).
    Structures with different residue ranges are outer-merged on res_id so
    every covered position is included.

    Args:
        results_dir: directory containing *_discotope3.csv files
                     (typically data/<target>/*discotope*/).
        pdb_dir:     directory containing the corresponding .pdb files
                     (typically data/<target>/pdb/).

    Returns:
        DataFrame with columns: res_id, residue, discotope_score,
        average_rsa, n_structures.
        discotope_score is the mean calibrated_score across all structures.
        Residue identity is taken from the first structure that covers each
        position.
    """
    pdb_stems = {p.stem.lower() for p in pdb_dir.glob("*.pdb")}
    csv_paths = [
        p for p in results_dir.glob("*.csv")
        if any(p.stem.lower().startswith(stem) for stem in pdb_stems)
    ]
    if not csv_paths:
        raise FileNotFoundError(
            f"No DiscoTope CSVs matching PDB stems in {pdb_dir} found under {results_dir}"
        )

    frames = [pd.read_csv(p) for p in csv_paths]
    struct_names = [p.stem for p in csv_paths]

    # Outer-merge all structures on res_id
    merged = frames[0][["res_id", "residue", "calibrated_score", "rsa"]].rename(
        columns={
            "calibrated_score": f"score_{struct_names[0]}",
            "rsa": f"rsa_{struct_names[0]}",
            "residue": f"residue_{struct_names[0]}",
        }
    )

    for name, df in zip(struct_names[1:], frames[1:]):
        right = df[["res_id", "residue", "calibrated_score", "rsa"]].rename(
            columns={
                "calibrated_score": f"score_{name}",
                "rsa": f"rsa_{name}",
                "residue": f"residue_{name}",
            }
        )
        merged = merged.merge(right, on="res_id", how="outer")

    score_cols = [c for c in merged.columns if c.startswith("score_")]
    rsa_cols   = [c for c in merged.columns if c.startswith("rsa_")]
    res_cols   = [c for c in merged.columns if c.startswith("residue_")]

    merged["discotope_score"] = merged[score_cols].mean(axis=1)
    merged["average_rsa"]     = merged[rsa_cols].mean(axis=1)
    merged["n_structures"]    = merged[score_cols].notna().sum(axis=1)
    # take residue identity from first structure that has it
    merged["residue"] = merged[res_cols[0]]
    for col in res_cols[1:]:
        merged["residue"] = merged["residue"].fillna(merged[col])

    return (
        merged[["res_id", "residue", "discotope_score", "average_rsa", "n_structures"]]
        .sort_values("res_id")
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def run(target_config: dict) -> None:
    """Parse local DiscoTope results and save to out_dir/discotope_raw.csv.

    Reads target_config keys:
        discotope_dir (Path): folder containing *_discotope3.csv files
        pdb_dir       (Path): folder containing matching .pdb files
        out_dir       (Path): destination for discotope_raw.csv
    """
    results_dir = Path(target_config["discotope_dir"])
    pdb_dir     = Path(target_config["pdb_dir"])
    out_dir     = Path(target_config["out_dir"])
    df = parse_results_dir(results_dir, pdb_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "discotope_raw.csv", index=False)


# ---------------------------------------------------------------------------
# Web scraper (not yet implemented)
# ---------------------------------------------------------------------------

def _run_web(
    pdb_path: Path,
    chain: str,
    cache_path: Path,
) -> pd.DataFrame:
    """Submit a PDB to the DiscoTope-3.0 web server and return results.

    Loads from cache if available, otherwise queries the web server.
    """
    if cache_path.exists():
        raw = json.loads(cache_path.read_text())
    else:
        raw = _query_server(pdb_path, chain)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(raw, indent=2))

    return _parse(raw, chain)


def _query_server(pdb_path: Path, chain: str) -> dict:
    raise NotImplementedError("DiscoTope web server query not yet implemented")


def _parse(raw: dict, chain: str) -> pd.DataFrame:
    raise NotImplementedError("DiscoTope result parser not yet implemented")
