"""GraphBepi local runner and output parser.

Runs graphbepi/test.py on cleaned single-chain PDB files.
Pipeline entry point: run(target_config)

Workflow per chain:
  1. Detect chains present in PDB; warn if multi-chain.
  2. Write ATOM-only, single-chain PDB to out_dir/clean_pdb/.
  3. Call test.py (-p -i <clean.pdb> -o <graphbepi_out/> --gpu 0 --threshold ...).
  4. Parse <pdb_stem>.csv output → standard schema.
  5. Save per-chain CSV (graphbepi_<chain>.csv) + combined graphbepi_raw.csv.

Output columns: res_id (int), chain (str), residue (str), score (float), is_epitope (bool)
"""

from __future__ import annotations

import subprocess
import sys
import warnings
from pathlib import Path

import pandas as pd
from Bio import PDB as biopdb


THRESHOLD = 0.1763          # GraphBepi model calibrated default
GRAPHBEPI_DIR = Path(__file__).parents[3] / "graphbepi"


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def run(target_config: dict) -> None:
    """Run GraphBepi on every PDB in pdb_dir; save per-PDB and combined CSVs.

    Iterates all .pdb files in pdb_dir and processes whichever chain is inside
    each file (auto-detected), ignoring the `chains` config field. This handles
    multi-structure folders where the antigen chain letter differs per PDB.

    Reads target_config keys:
        pdb_dir (Path): folder containing .pdb files for the target
        out_dir (Path): destination for output CSVs
    """
    pdb_dir = Path(target_config["pdb_dir"])
    out_dir = Path(target_config["out_dir"])

    cache_path = out_dir / "graphbepi_raw.csv"
    if cache_path.exists():
        return

    pdb_files = sorted(pdb_dir.glob("*.pdb"))
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in {pdb_dir}")

    chains_filter: set[str] | None = (
        set(target_config["chains"]) if target_config.get("chains") else None
    )

    clean_pdb_dir     = out_dir / "clean_pdb"
    graphbepi_out_dir = out_dir / "graphbepi_out"
    clean_pdb_dir.mkdir(parents=True, exist_ok=True)
    graphbepi_out_dir.mkdir(parents=True, exist_ok=True)

    all_dfs = []
    for pdb_path in pdb_files:
        detected = _detect_chains(pdb_path)
        if not detected:
            warnings.warn(f"{pdb_path.name}: no ATOM records found, skipping.", stacklevel=2)
            continue

        to_process = sorted(detected & chains_filter if chains_filter else detected)
        if not to_process:
            warnings.warn(
                f"{pdb_path.name}: none of the configured chains {chains_filter} "
                f"found in {detected}, skipping.",
                stacklevel=2,
            )
            continue
        if len(detected) > 1 and not chains_filter:
            warnings.warn(
                f"{pdb_path.name} contains multiple chains {detected}; "
                f"processing each separately.",
                stacklevel=2,
            )

        for chain_id in to_process:
            pdb_stem   = pdb_path.stem
            clean_name = pdb_stem if pdb_stem.endswith(f"_{chain_id}") else f"{pdb_stem}_{chain_id}"
            clean_pdb  = clean_pdb_dir / f"{clean_name}.pdb"
            _clean_pdb(pdb_path, chain_id, clean_pdb)

            _run_graphbepi(clean_pdb, graphbepi_out_dir)

            raw_csv = graphbepi_out_dir / f"{clean_name}.csv"
            res_ids = _pdb_res_ids(clean_pdb, chain_id)
            df = _parse(raw_csv, chain_id, res_ids)

            per_chain_path = out_dir / f"graphbepi_{chain_id}.csv"
            df.to_csv(per_chain_path, index=False)

            all_dfs.append(df)

    if not all_dfs:
        raise RuntimeError("GraphBepi produced no output for any configured chain.")

    combined = pd.concat(all_dfs, ignore_index=True)
    combined.to_csv(cache_path, index=False)


# ---------------------------------------------------------------------------
# PDB helpers
# ---------------------------------------------------------------------------

def _pdb_res_ids(pdb_path: Path, chain_id: str) -> list[int]:
    """Return PDB residue sequence numbers for each unique residue in chain_id.

    Reads CA atoms in file order; one entry per residue. These numbers match
    the res_id values used by DiscoTope, allowing the outputs to be merged.
    """
    seen: list[int] = []
    seen_set: set[tuple[int, str]] = set()
    with open(pdb_path) as fh:
        for line in fh:
            if line[:4] != "ATOM":
                continue
            if line[21] != chain_id:
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            res_seq = int(line[22:26].strip())
            ins_code = line[26].strip()
            key = (res_seq, ins_code)
            if key not in seen_set:
                seen_set.add(key)
                seen.append(res_seq)
    return seen


def _detect_chains(pdb_path: Path) -> set[str]:
    """Return the set of chain IDs present in ATOM records."""
    chains: set[str] = set()
    with open(pdb_path) as fh:
        for line in fh:
            if line[:6] in ("ATOM  ", "ATOM"):
                if len(line) > 21:
                    chains.add(line[21])
    return chains


def _clean_pdb(pdb_path: Path, chain_id: str, out_path: Path) -> None:
    """Write a normalised PDB with only ATOM records for chain_id.

    Uses biopython PDBIO to rewrite the file with correct column alignment,
    which fixes B-factor overflow and other formatting issues that mkdssp 4.x
    rejects.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    class _ChainATOMSelect(biopdb.Select):
        def accept_model(self, model):
            return model.id == 0  # first model only (handles NMR multi-model PDBs)
        def accept_chain(self, chain):
            return chain.id == chain_id
        def accept_atom(self, atom):
            return not atom.is_disordered() or atom.get_altloc() == "A"

    parser = biopdb.PDBParser(QUIET=True)
    structure = parser.get_structure("target", str(pdb_path))
    io = biopdb.PDBIO()
    io.set_structure(structure)
    io.save(str(out_path), select=_ChainATOMSelect())


# ---------------------------------------------------------------------------
# Local runner
# ---------------------------------------------------------------------------

def _run_graphbepi(pdb_path: Path, graphbepi_out_dir: Path) -> None:
    """Execute graphbepi/test.py for a cleaned single-chain PDB."""
    cmd = [
        sys.executable,
        str(GRAPHBEPI_DIR / "test.py"),
        "-p",                               # PDB mode (boolean flag)
        "-i", str(pdb_path),                # input PDB
        "-o", str(graphbepi_out_dir),       # output directory
        "--gpu", "0",                       # 0 → MPS on Apple Silicon, cuda:0 on Linux
        "--threshold", str(THRESHOLD),
    ]
    result = subprocess.run(cmd, cwd=GRAPHBEPI_DIR)
    if result.returncode != 0:
        raise RuntimeError(
            f"GraphBepi failed for {pdb_path.name} (see output above)"
        )


# ---------------------------------------------------------------------------
# Output parser
# ---------------------------------------------------------------------------

def _parse(csv_path: Path, chain: str, res_ids: list[int]) -> pd.DataFrame:
    """Parse test.py output CSV into standard pipeline schema.

    Input columns (from graphbepi/output/7s2r_A.csv):
        resn, score, is epitope
    Output columns:
        res_id (PDB residue number), chain, residue, score, is_epitope

    res_ids must be the PDB residue sequence numbers for each row, in order,
    as returned by _pdb_res_ids(). This ensures res_id matches DiscoTope numbering.
    """
    df = pd.read_csv(csv_path)
    df = df.rename(columns={"resn": "residue", "is epitope": "is_epitope"})
    if len(df) != len(res_ids):
        raise ValueError(
            f"{csv_path.name}: {len(df)} rows but {len(res_ids)} residues in PDB"
        )
    df["res_id"]     = res_ids
    df["chain"]      = chain
    df["is_epitope"] = df["is_epitope"].astype(bool)
    return df[["res_id", "chain", "residue", "score", "is_epitope"]]
