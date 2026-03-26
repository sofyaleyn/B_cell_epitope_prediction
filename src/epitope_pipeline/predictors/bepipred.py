"""BepiPred web server scraper and local result parser.

Web API (run): submits a FASTA sequence to the BepiPred-3.0 web server,
polls for results, caches raw JSON to outputs/<target>/bepipred_raw.json.

Local parser (parse_results_dir): reads a pre-downloaded raw_output.csv
from a *bepipred* results directory and returns a standardised DataFrame.

Output columns: res_id (int), residue (str), bepipred_score (float)
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import requests


THRESHOLD = 0.5


# ---------------------------------------------------------------------------
# Local result parser
# ---------------------------------------------------------------------------

def parse_results_dir(results_dir: Path) -> pd.DataFrame:
    """Parse a BepiPred-3.0 raw_output.csv into a standardised DataFrame.

    Args:
        results_dir: directory containing raw_output.csv
                     (typically data/<target>/*bepipred*/).

    Returns:
        DataFrame with columns: res_id, residue, bepipred_score.
        res_id is 1-based, matching the order of residues in the CSV.
        Only the first accession is returned (the target protein).
    """
    csv_path = results_dir / "raw_output.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"BepiPred raw_output.csv not found in {results_dir}")

    df = pd.read_csv(csv_path)
    first_accession = df["Accession"].iloc[0]
    df = df[df["Accession"] == first_accession].copy()

    df = df.rename(columns={
        "Residue": "residue",
        "BepiPred-3.0 score": "bepipred_score",
    })
    df["res_id"] = range(1, len(df) + 1)

    return df[["res_id", "residue", "bepipred_score"]].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def run(target_config: dict) -> None:
    """Parse local BepiPred results and save to out_dir/bepipred_raw.csv.

    Reads target_config keys:
        bepipred_dir (Path): folder containing raw_output.csv
        out_dir      (Path): destination for bepipred_raw.csv
    """
    results_dir = Path(target_config["bepipred_dir"])
    out_dir     = Path(target_config["out_dir"])
    df = parse_results_dir(results_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "bepipred_raw.csv", index=False)


# ---------------------------------------------------------------------------
# Web scraper (not yet implemented)
# ---------------------------------------------------------------------------

def _run_web(
    sequence: str,
    chain: str,
    cache_path: Path,
) -> pd.DataFrame:
    """Submit sequence to the BepiPred-3.0 web server and return results.

    Loads from cache if available, otherwise queries the web server.
    """
    if cache_path.exists():
        raw = json.loads(cache_path.read_text())
    else:
        raw = _query_server(sequence)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(raw, indent=2))

    return _parse(raw, chain)


def _query_server(sequence: str) -> dict:
    raise NotImplementedError("BepiPred web server query not yet implemented")


def _parse(raw: dict, chain: str) -> pd.DataFrame:
    raise NotImplementedError("BepiPred result parser not yet implemented")
