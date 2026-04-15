"""Summary report export.

export_html: Produces a human-readable HTML summary of predicted epitope
             regions per tool, consensus overlap, and (optionally) the
             antibody-antigen interface.  No raw data tables.

export_pdf:  Multi-page PDF via matplotlib PdfPages (plots + lightweight
             region summary).  Use --format pdf; note plot resolution is
             limited by the PNG DPI of the source images.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from jinja2 import Environment, BaseLoader


# ── HTML summary ──────────────────────────────────────────────────────────────

_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>{{ target_name }} — Epitope Prediction Report</title>
  <style>
    body { font-family: Georgia, serif; max-width: 860px; margin: 2.5rem auto; color: #222; line-height: 1.65; }
    h1   { font-size: 1.55rem; margin-bottom: 0.15rem; }
    h2   { font-size: 1.1rem; margin-top: 2rem; border-bottom: 1px solid #ccc; padding-bottom: 4px; }
    h3   { font-size: 0.98rem; margin-top: 1.2rem; margin-bottom: 0.25rem; }
    p    { margin-top: 0.35rem; }
    ul   { margin-top: 0.3rem; padding-left: 1.4rem; }
    li   { margin-bottom: 0.2rem; font-family: 'Courier New', monospace; font-size: 0.92rem; }
    .meta   { font-size: 0.83rem; color: #666; }
    .none   { color: #888; font-style: italic; }
    .box    { background: #f0f7ff; border-left: 3px solid #4a90d9;
              padding: 0.65rem 1rem; margin: 1rem 0; border-radius: 0 4px 4px 0; }
    .warn   { background: #fff8e1; border-left: 3px solid #e6a817;
              padding: 0.65rem 1rem; margin: 1rem 0; border-radius: 0 4px 4px 0; }
  </style>
</head>
<body>
  <h1>{{ target_name }} — B-cell Epitope Prediction Summary</h1>
  <p class="meta">{{ n_residues }} residues in full sequence</p>

  <!-- ── BepiPred ── -->
  <h2>BepiPred 3.0 <span class="meta">(sequence-based · threshold 0.50)</span></h2>
  {% if bp_regions %}
  <p>{{ bp_regions | length }} region{{ 's' if bp_regions | length != 1 else '' }},
     {{ bp_total }} residues:</p>
  <ul>{% for r in bp_regions %}<li>{{ r }}</li>{% endfor %}</ul>
  {% else %}
  <p class="none">No epitope regions predicted above threshold.</p>
  {% endif %}

  <!-- ── DiscoTope ── -->
  <h2>DiscoTope 3.0 <span class="meta">(structure-based · calibrated score threshold 0.90)</span></h2>
  {% if dt_regions %}
  <p>{{ dt_regions | length }} region{{ 's' if dt_regions | length != 1 else '' }},
     {{ dt_total }} residues:</p>
  <ul>{% for r in dt_regions %}<li>{{ r }}</li>{% endfor %}</ul>
  {% else %}
  <p class="none">No epitope regions predicted above threshold.</p>
  {% endif %}

  <!-- ── GraphBepi ── -->
  <h2>GraphBepi <span class="meta">(graph neural network · structure-based · threshold 0.1763)</span></h2>
  {% if gb_regions %}
  <p>{{ gb_regions | length }} region{{ 's' if gb_regions | length != 1 else '' }},
     {{ gb_total }} residues:</p>
  <ul>{% for r in gb_regions %}<li>{{ r }}</li>{% endfor %}</ul>
  {% else %}
  <p class="none">No epitope regions predicted above threshold.</p>
  {% endif %}

  <!-- ── Consensus ── -->
  <h2>Consensus — any structure-based predictor (DiscoTope ∪ GraphBepi)</h2>
  {% if and_regions %}
  <div class="box">
    <p><strong>{{ and_regions | length }} region{{ 's' if and_regions | length != 1 else '' }},
       {{ and_total }} residues</strong> flagged by at least one structure-based predictor:</p>
    <ul>{% for r in and_regions %}<li>{{ r }}</li>{% endfor %}</ul>
    {% if and_both_regions %}
    <p style="margin-top:0.8rem"><strong>Both DiscoTope and GraphBepi agree:</strong>
       {{ and_both_regions | length }} region{{ 's' if and_both_regions | length != 1 else '' }},
       {{ and_both_total }} residues:</p>
    <ul>{% for r in and_both_regions %}<li>{{ r }}</li>{% endfor %}</ul>
    {% endif %}
  </div>
  {% else %}
  <p class="none">No consensus regions found.</p>
  {% endif %}

  <!-- ── Interface ── -->
  {% if has_interface %}
  <h2>Antibody–Antigen Interface <span class="meta">(crystal structure contacts)</span></h2>

  <h3>Direct contacts — {{ iface_direct_total }} residue{{ 's' if iface_direct_total != 1 else '' }}</h3>
  {% if iface_direct_regions %}
  <ul>{% for r in iface_direct_regions %}<li>{{ r }}</li>{% endfor %}</ul>
  {% else %}<p class="none">None recorded.</p>{% endif %}

  <h3>Adjacent contacts — {{ iface_adj_total }} residue{{ 's' if iface_adj_total != 1 else '' }}</h3>
  {% if iface_adj_regions %}
  <ul>{% for r in iface_adj_regions %}<li>{{ r }}</li>{% endfor %}</ul>
  {% else %}<p class="none">None recorded.</p>{% endif %}

  {% if overlap_lines %}
  <h3>Overlap with predicted epitopes</h3>
  {% for line in overlap_lines %}<p>{{ line }}</p>{% endfor %}
  {% endif %}
  {% endif %}

</body>
</html>
"""


def export_html(
    df: pd.DataFrame,
    target_name: str,
    output_path: Path,
    interface_df: pd.DataFrame | None = None,
) -> None:
    """Render a text-summary HTML report and write it to output_path."""
    env = Environment(loader=BaseLoader())
    template = env.from_string(_TEMPLATE)

    bp_raw  = _get_regions(df, "is_epitope_bepipred")  if "is_epitope_bepipred"  in df.columns else []
    dt_raw  = _get_regions(df, "is_epitope_discotope") if "is_epitope_discotope" in df.columns else []
    gb_raw  = _get_regions(df, "is_epitope_graphbepi") if "is_epitope_graphbepi" in df.columns else []
    and_raw = _get_regions(df, "is_epitope_AND")       if "is_epitope_AND"       in df.columns else []

    # Residues flagged by BOTH DiscoTope and GraphBepi (true intersection)
    both_col = None
    if "is_epitope_discotope" in df.columns and "is_epitope_graphbepi" in df.columns:
        df = df.copy()
        df["_both"] = df["is_epitope_discotope"] & df["is_epitope_graphbepi"]
        both_raw = _get_regions(df, "_both")
    else:
        both_raw = []

    # Interface regions
    has_interface = interface_df is not None and not interface_df.empty
    iface_direct_raw = iface_adj_raw = []
    overlap_lines = []
    if has_interface:
        direct_df = interface_df[interface_df["contact_type"] == "direct"]
        adj_df    = interface_df[interface_df["contact_type"] == "adjacent"]
        iface_direct_raw = _get_regions_from_series(direct_df["res_id"])
        iface_adj_raw    = _get_regions_from_series(adj_df["res_id"])
        overlap_lines    = _overlap_summary(df, interface_df)

    ctx = dict(
        target_name=target_name,
        n_residues=len(df),
        bp_regions=_fmt(bp_raw),   bp_total=_total(bp_raw),
        dt_regions=_fmt(dt_raw),   dt_total=_total(dt_raw),
        gb_regions=_fmt(gb_raw),   gb_total=_total(gb_raw),
        and_regions=_fmt(and_raw), and_total=_total(and_raw),
        and_both_regions=_fmt(both_raw), and_both_total=_total(both_raw),
        has_interface=has_interface,
        iface_direct_regions=_fmt(iface_direct_raw),
        iface_direct_total=_total(iface_direct_raw),
        iface_adj_regions=_fmt(iface_adj_raw),
        iface_adj_total=_total(iface_adj_raw),
        overlap_lines=overlap_lines,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(template.render(**ctx), encoding="utf-8")


# ── helpers ───────────────────────────────────────────────────────────────────

def _get_regions(df: pd.DataFrame, flag_col: str) -> list[tuple[int, int, str]]:
    """Group consecutive flagged residues into (start_id, end_id, sequence) tuples."""
    flagged = df[df[flag_col]].sort_values("res_id")
    if flagged.empty:
        return []
    return _consecutive(flagged["res_id"].tolist(), flagged["residue"].tolist())


def _get_regions_from_series(res_ids: pd.Series) -> list[tuple[int, int, str]]:
    """Group consecutive res_id values (no sequence available) into region tuples."""
    ids = sorted(res_ids.dropna().astype(int).tolist())
    if not ids:
        return []
    return _consecutive(ids, ["?" for _ in ids])


def _consecutive(ids: list, aas: list) -> list[tuple[int, int, str]]:
    regions: list[tuple[int, int, str]] = []
    start = ids[0]
    prev  = ids[0]
    seq   = [aas[0]]
    for rid, aa in zip(ids[1:], aas[1:]):
        if rid > prev + 1:
            regions.append((start, prev, "".join(seq)))
            start, seq = rid, [aa]
        else:
            seq.append(aa)
        prev = rid
    regions.append((start, prev, "".join(seq)))
    return regions


def _fmt(regions: list[tuple[int, int, str]]) -> list[str]:
    out = []
    for start, end, seq in regions:
        n = end - start + 1
        if start == end:
            out.append(f"residue {start} ({seq})")
        else:
            out.append(f"residues {start}–{end}  ({n} aa)  {seq}")
    return out


def _total(regions: list[tuple[int, int, str]]) -> int:
    return sum(end - start + 1 for start, end, _ in regions)


def _overlap_summary(df: pd.DataFrame, interface_df: pd.DataFrame) -> list[str]:
    """Return prose lines describing how many interface residues fall in predicted epitope regions."""
    epitope_ids = set(df.loc[df.get("is_epitope_AND", pd.Series(False, index=df.index)), "res_id"].tolist())
    lines = []
    for contact_type in ("direct", "adjacent"):
        sub = interface_df[interface_df["contact_type"] == contact_type]
        total = len(sub)
        if total == 0:
            continue
        in_epitope = sub[sub["res_id"].isin(epitope_ids)]
        n = len(in_epitope)
        pct = round(100 * n / total)
        missed = sub[~sub["res_id"].isin(epitope_ids)]["res_id"].tolist()
        line = (
            f"{contact_type.capitalize()} contacts: "
            f"{n}/{total} residues ({pct}%) fall within predicted epitope regions."
        )
        if missed:
            line += f"  Missed: {', '.join(str(r) for r in sorted(missed))}."
        lines.append(line)
    return lines


# ── PDF export (unchanged) ────────────────────────────────────────────────────

_ROWS_PER_PAGE = 40
_PAGE_SIZE = (11, 8.5)


def export_pdf(
    df: pd.DataFrame,
    target_name: str,
    output_path: Path,
    discotope_per_struct_img: Path | None = None,
    graphbepi_per_struct_img: Path | None = None,
) -> None:
    """Render a multi-page PDF report and write it to output_path."""
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(output_path) as pdf:
        fig, ax = plt.subplots(figsize=_PAGE_SIZE)
        ax.axis("off")
        ax.text(0.5, 0.55, target_name, ha="center", va="center",
                fontsize=28, fontweight="bold", transform=ax.transAxes)
        ax.text(0.5, 0.45, "B-cell Epitope Prediction Report",
                ha="center", va="center", fontsize=18, transform=ax.transAxes)
        pdf.savefig(fig)
        plt.close(fig)

        combined_png = output_path.parent / "B_cell_epitope_combined.png"
        plot_pages = [
            (combined_png, "Combined predictor scores"),
            (discotope_per_struct_img, "DiscoTope 3.0 — per-structure predictions"),
            (graphbepi_per_struct_img, "GraphBepi — per-structure predictions"),
        ]
        for img_path, title in plot_pages:
            if img_path and Path(img_path).exists():
                img = mpimg.imread(str(img_path))
                fig, ax = plt.subplots(figsize=_PAGE_SIZE)
                ax.imshow(img)
                ax.axis("off")
                ax.set_title(title, fontsize=12, pad=6)
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
