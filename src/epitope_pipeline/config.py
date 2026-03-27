"""Target config resolution for the epitope prediction pipeline.

resolve_target_config(target_name, repo_root) is the single entry point used
by all pipeline scripts. It:
  1. Looks up the target in config/targets.yaml (explicit config wins).
  2. If not found, auto-detects from data/<target_name>/ and writes the entry
     back to targets.yaml for future runs.
"""

from __future__ import annotations

from pathlib import Path

import yaml


CONFIG_PATH_REL = Path("config") / "targets.yaml"


def resolve_target_config(target_name: str, repo_root: Path) -> dict:
    """Return target config dict with all paths resolved to absolute.

    If the target is not in targets.yaml, auto-detects from data/<target_name>/
    and registers it in targets.yaml before returning.

    Args:
        target_name: name of the target (must match a folder under data/).
        repo_root:   absolute path to the repository root.

    Returns:
        Dict with keys: chains, pdb_dir, and optionally bepipred_dir,
        discotope_dir. All path values are absolute Path objects.
        Also includes out_dir (outputs/<target_name>/).

    Raises:
        FileNotFoundError: if data/<target_name>/ does not exist or has no PDB.
        SystemExit: if target is not in yaml and auto-detection also fails.
    """
    config_path = repo_root / CONFIG_PATH_REL
    config = yaml.safe_load(config_path.read_text())

    if target_name in config["targets"]:
        return _resolve_paths(config["targets"][target_name], target_name, repo_root)

    # Auto-detect from data/<target_name>/
    data_dir = repo_root / "data" / target_name
    if not data_dir.exists():
        raise SystemExit(
            f"Unknown target '{target_name}' and data directory not found: {data_dir}\n"
            f"Create data/{target_name}/ with a FASTA file and pdb/ subfolder."
        )

    detected = _auto_detect(target_name, data_dir, repo_root)
    _write_to_yaml(target_name, detected, config, config_path)
    print(f"[{target_name}] Auto-detected config registered in {CONFIG_PATH_REL}")

    return _resolve_paths(detected, target_name, repo_root)


# ---------------------------------------------------------------------------
# Auto-detection
# ---------------------------------------------------------------------------

def _auto_detect(target_name: str, data_dir: Path, repo_root: Path) -> dict:
    """Infer target config from the data folder structure."""
    # pdb_dir: prefer data/<target>/pdb/, fall back to data/<target>/ itself
    pdb_subdir = data_dir / "pdb"
    if pdb_subdir.is_dir() and list(pdb_subdir.glob("*.pdb")):
        pdb_dir = pdb_subdir
    else:
        pdbs = list(data_dir.glob("*.pdb"))
        if not pdbs:
            raise FileNotFoundError(
                f"No .pdb files found under {data_dir} or {pdb_subdir}. "
                f"Place PDB files in data/{target_name}/pdb/."
            )
        pdb_dir = data_dir

    # chains: read ATOM records from first PDB
    first_pdb = sorted(pdb_dir.glob("*.pdb"))[0]
    chains = _detect_chains(first_pdb)
    if not chains:
        raise FileNotFoundError(f"No ATOM records found in {first_pdb}.")

    entry: dict = {
        "chains": chains,
        "pdb_dir": str(pdb_dir.relative_to(repo_root)),
    }

    # optional dirs
    bepipred_dir = _find_optional_subdir(data_dir, "*bepipred*")
    if bepipred_dir:
        entry["bepipred_dir"] = str(bepipred_dir.relative_to(repo_root))

    discotope_dir = _find_optional_subdir(data_dir, "*discotope*")
    if discotope_dir:
        entry["discotope_dir"] = str(discotope_dir.relative_to(repo_root))

    return entry


def _detect_chains(pdb_path: Path) -> list[str]:
    chains: set[str] = set()
    with open(pdb_path) as fh:
        for line in fh:
            if line[:4] == "ATOM" and len(line) > 21:
                chains.add(line[21])
    return sorted(chains)


def _find_optional_subdir(data_dir: Path, pattern: str) -> Path | None:
    matches = [p for p in data_dir.iterdir() if p.is_dir() and p.match(pattern)]
    return matches[0] if matches else None


# ---------------------------------------------------------------------------
# YAML write-back
# ---------------------------------------------------------------------------

def _write_to_yaml(target_name: str, entry: dict, config: dict, config_path: Path) -> None:
    """Add the new target entry to targets.yaml, preserving existing content."""
    config["targets"][target_name] = entry
    with open(config_path, "w") as fh:
        yaml.dump(config, fh, default_flow_style=False, sort_keys=False, allow_unicode=True)


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

def _resolve_paths(entry: dict, target_name: str, repo_root: Path) -> dict:
    """Return a copy of entry with all path values resolved to absolute Paths."""
    cfg = dict(entry)
    for key in ("bepipred_dir", "discotope_dir", "pdb_dir"):
        if key in cfg:
            cfg[key] = repo_root / cfg[key]
    cfg["out_dir"] = repo_root / "outputs" / target_name
    cfg["out_dir"].mkdir(parents=True, exist_ok=True)
    return cfg
