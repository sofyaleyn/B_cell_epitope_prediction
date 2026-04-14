"""BepiPred web server scraper and local result parser.

Web mode (run, no bepipred_dir): submits a FASTA sequence to the BepiPred-3.0
web server, polls for results, caches raw CSV to outputs/<target>/bepipred_cache.csv,
and saves the normalised output to outputs/<target>/bepipred_raw.csv.

Local parser (parse_results_dir): reads a pre-downloaded raw_output.csv
from a *bepipred* results directory and returns a standardised DataFrame.

Output columns: res_id (int), residue (str), bepipred_score (float)
"""

from __future__ import annotations

import time
from io import StringIO
from pathlib import Path

import pandas as pd
import requests
from bs4 import BeautifulSoup


THRESHOLD = 0.5

BASE_URL = "https://services.healthtech.dtu.dk"
SUBMIT_URL = f"{BASE_URL}/cgi-bin/webface2.fcgi"


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
    """Parse BepiPred results (local or web) and save to out_dir/bepipred_raw.csv.

    Local mode (bepipred_dir present in config):
        Reads raw_output.csv from bepipred_dir.

    Web mode (no bepipred_dir):
        Discovers the first *.fasta in target_dir, submits to the BepiPred-3.0
        web server, caches raw CSV to out_dir/bepipred_cache.csv.
    """
    out_dir = Path(target_config["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    if target_config.get("bepipred_dir"):
        df = parse_results_dir(Path(target_config["bepipred_dir"]))
    else:
        target_dir = Path(target_config["target_dir"])
        fasta_files = sorted(target_dir.glob("*.fasta"))
        if not fasta_files:
            raise FileNotFoundError(
                f"No .fasta file found in {target_dir} — "
                "either add a FASTA or set bepipred_dir in targets.yaml"
            )
        fasta_text = fasta_files[0].read_text()
        cache_path = out_dir / "bepipred_cache.csv"
        df = _run_web(fasta_text, cache_path)

    df.to_csv(out_dir / "bepipred_raw.csv", index=False)


# ---------------------------------------------------------------------------
# Web scraper
# ---------------------------------------------------------------------------

def _run_web(fasta_text: str, cache_path: Path) -> pd.DataFrame:
    """Submit FASTA to BepiPred-3.0 and return normalised DataFrame.

    Loads from cache if available, otherwise queries the web server.
    """
    if cache_path.exists():
        print(f"[bepipred] Using cached results: {cache_path}")
        return _parse_csv_text(cache_path.read_text())

    print("[bepipred] Submitting to BepiPred-3.0 web server...")
    csv_text = _query_server(fasta_text)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(csv_text)
    print(f"[bepipred] Cached raw results → {cache_path}")
    return _parse_csv_text(csv_text)


def _query_server(fasta_text: str) -> str:
    """Submit FASTA to BepiPred-3.0, poll until done, return raw CSV text."""
    session = requests.Session()

    # Step 1: GET submission page to read hidden form fields and form action
    page = session.get(f"{BASE_URL}/services/BepiPred-3.0/", timeout=30)
    page.raise_for_status()
    soup = BeautifulSoup(page.text, "html.parser")

    form = soup.find("form")
    if form is None:
        raise RuntimeError("Could not find submission form on BepiPred-3.0 page")

    hidden = {
        inp["name"]: inp.get("value", "")
        for inp in form.find_all("input", type="hidden")
        if inp.get("name")
    }
    action = form.get("action") or SUBMIT_URL
    if not action.startswith("http"):
        action = BASE_URL + ("" if action.startswith("/") else "/") + action

    # Find the textarea name for the FASTA input
    textarea = form.find("textarea")
    fasta_field = textarea["name"] if textarea and textarea.get("name") else "SEQPASTE"

    # Step 2: POST FASTA as multipart/form-data (matches form enctype)
    # Field names confirmed from form inspection: thr_epi, roll_mean, top_epi
    fields = {
        **{k: (None, v) for k, v in hidden.items()},
        fasta_field: (None, fasta_text),
        "thr_epi": (None, "0.5"),
        "roll_mean": (None, "yes"),
    }
    resp = session.post(action, files=fields, timeout=60)
    resp.raise_for_status()

    # Step 3: Extract job ID — DTU redirects to ?jobid=<id> or uses meta-refresh
    job_url = resp.url
    if "jobid" not in job_url:
        soup2 = BeautifulSoup(resp.text, "html.parser")
        meta = soup2.find("meta", attrs={"http-equiv": "refresh"})
        if meta and "url=" in meta.get("content", "").lower():
            rel = meta["content"].split("url=")[-1].strip()
            job_url = rel if rel.startswith("http") else BASE_URL + "/" + rel.lstrip("/")

    if "jobid" not in job_url:
        # Print response snippet to help diagnose the actual server response
        snippet = resp.text[:3000].replace("\n", " ")
        raise RuntimeError(
            f"Could not extract job ID from BepiPred-3.0 response.\n"
            f"Response URL: {resp.url}\n"
            f"Response snippet:\n{snippet}\n\n"
            "Look for 'jobid' in the snippet above to find where it is embedded."
        )

    jobid = job_url.split("jobid=")[-1].split("&")[0]
    print(f"[bepipred] Job submitted: {jobid}")

    # Step 4: Poll until finished (max ~10 minutes)
    poll_url = f"{SUBMIT_URL}?jobid={jobid}&wait=20"
    for attempt in range(30):
        time.sleep(20)
        poll = session.get(poll_url, timeout=60)
        poll.raise_for_status()
        text_lower = poll.text.lower()
        if "raw_output.csv" in poll.text or "finished" in text_lower or "result" in text_lower:
            break
        print(f"[bepipred] Waiting for results (attempt {attempt + 1}/30)...")
    else:
        raise TimeoutError("BepiPred-3.0 job did not complete within 10 minutes")

    # Step 5: Download raw_output.csv
    csv_url = _find_csv_link(poll.text, jobid)
    csv_resp = session.get(csv_url, timeout=60)
    csv_resp.raise_for_status()
    return csv_resp.text


def _find_csv_link(html: str, jobid: str) -> str:
    """Extract the raw_output.csv download URL from the results page HTML."""
    soup = BeautifulSoup(html, "html.parser")
    for a in soup.find_all("a", href=True):
        if "raw_output.csv" in a["href"]:
            href = a["href"]
            return href if href.startswith("http") else BASE_URL + "/" + href.lstrip("/")
    # Fallback: DTU standard download path
    return f"{SUBMIT_URL}?jobid={jobid}&outfile=raw_output.csv"


def _parse_csv_text(csv_text: str) -> pd.DataFrame:
    """Parse raw_output.csv text into the normalised bepipred DataFrame."""
    df = pd.read_csv(StringIO(csv_text))
    first_accession = df["Accession"].iloc[0]
    df = df[df["Accession"] == first_accession].copy()
    df = df.rename(columns={
        "Residue": "residue",
        "BepiPred-3.0 score": "bepipred_score",
    })
    df["res_id"] = range(1, len(df) + 1)
    return df[["res_id", "residue", "bepipred_score"]].reset_index(drop=True)
