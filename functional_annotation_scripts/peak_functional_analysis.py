import glob
import os
from collections import Counter
from io import StringIO

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# config

GFF_PATH = "genome/kieran_braker3_Sratti.gff3"
EMAPPER_PATH = "genome/strongyloides_ratti.emapper.annotations.tsv"
OUT_DIR = "functional_annotation"

PEAK_PREFIXES = [
    "FLK4", "FLK9", "FLK27", "FLK36", "FLK92",
    "FL2K4", "FL2K9", "FL2K27", "FL2K36", "FL2K92",
    "PFK4", "PFK9", "PFK27", "PFK36", "PFK92",
]

GENE_FLANK = 50

ALTERNATIVE = "greater"  # one-sided: testing for over-representation only
SIG_ALPHA_PLOT = 0.05

os.makedirs(OUT_DIR, exist_ok=True)

def get_gene_id(attrs):
    for part in attrs.split(";"):
        if part.startswith("ID="):
            return part.replace("ID=", "")
    return None


def load_background(path):
    with open(path) as f:
        lines = f.readlines()
    header_idx = next(i for i, l in enumerate(lines) if l.startswith("#query"))
    header = lines[header_idx].lstrip("#").strip().split("\t")
    data_lines = [l for l in lines[header_idx + 1:] if not l.startswith("#") and l.strip()]
    df = pd.read_csv(StringIO("".join(data_lines)), sep="\t", header=None, names=header)
    print(f"Loaded background: {len(df)} genes from {path}")
    return df


def explode_terms_pfam(df, gene_col="gene_id", term_col="PFAMs"):
    rows = []
    for gene, val in zip(df[gene_col], df[term_col]):
        if pd.isna(val) or val in ("-", ""):
            continue
        for t in val.split(","):
            t = t.strip()
            if t and t != "-":
                rows.append((gene, t))
    return pd.DataFrame(rows, columns=["gene", "term"])


def benjamini_hochberg(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    bh = ranked * n / (np.arange(n) + 1)
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    bh = np.clip(bh, 0, 1)
    out = np.empty(n)
    out[order] = bh
    return out


def run_fisher(a, group_total, c, bg_total):
    table = [[a, group_total - a], [c, bg_total - c]]
    _, p = fisher_exact(table, alternative=ALTERNATIVE)
    fold = (
        (a / group_total) / (c / bg_total)
        if c > 0 and group_total > 0
        else float("inf") if a > 0 else float("nan")
    )
    return (round(fold, 2) if np.isfinite(fold) else fold), p


def pfam_stats_for_geneset(gene_ids, genes_df, emapper, background_term_counts,
                            background_total, pfam_to_desc, dataset_label):
    subset = genes_df[genes_df["gene_id"].isin(gene_ids)].copy()
    annot = subset.merge(emapper, on="gene_id", how="left")
    n_annotated = annot["Description"].notna().sum()
    pct = round(n_annotated / len(annot) * 100, 1) if len(annot) else 0
    print(f"  {dataset_label}: {len(annot)} genes, {n_annotated} annotated ({pct}%)")

    pfam_counts = Counter()
    for val in annot["PFAMs"].dropna():
        for domain in val.split(","):
            domain = domain.strip()
            if domain and domain != "-":
                pfam_counts[domain] += 1

    if not pfam_counts:
        print(f"  {dataset_label}: no PFAM domains found -- skipping stats.")
        return annot, None

    group_total = len(annot)
    terms = list(pfam_counts.keys())
    folds, pvals = [], []
    for term in terms:
        a = pfam_counts[term]
        c = background_term_counts.get(term, 0)
        fold, p = run_fisher(a, group_total, c, background_total)
        folds.append(fold)
        pvals.append(p)
    padj = benjamini_hochberg(pvals)

    stats_df = pd.DataFrame({
        "dataset": dataset_label,
        "PFAM": terms,
        "domain_name": [pfam_to_desc.get(t, "") for t in terms],
        "n_genes": [pfam_counts[t] for t in terms],
        "group_total": group_total,
        "background_n": [background_term_counts.get(t, 0) for t in terms],
        "background_total": background_total,
        "fold_enrichment": folds,
        "pvalue": pvals,
        "padj_bh": padj,
    })
    return annot, stats_df


def make_dotplot_and_heatmap(pfam_all, dataset_order, out_prefix, title_stub):
    sig_domains = pfam_all.loc[pfam_all["padj_bh"] < SIG_ALPHA_PLOT, "PFAM"].unique()
    plot_df = pfam_all[pfam_all["PFAM"].isin(sig_domains)].copy()

    if plot_df.empty:
        print(f"No PFAM domains significant at padj_bh < {SIG_ALPHA_PLOT} for {title_stub} -- skipping plots.")
        return

    plot_df["label"] = plot_df.apply(
        lambda r: f"{r['PFAM']} ({r['domain_name']})" if r["domain_name"] else r["PFAM"], axis=1
    )
    domain_order = plot_df.groupby("label")["padj_bh"].min().sort_values().index.tolist()
    plot_df["neg_log10_padj"] = -np.log10(plot_df["padj_bh"].clip(lower=1e-300))

    fig_w = 1.6 * len(dataset_order) + 4
    fig_h = 0.35 * len(domain_order) + 2
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    x_pos = {d: i for i, d in enumerate(dataset_order)}
    y_pos = {d: i for i, d in enumerate(domain_order)}
    xs = plot_df["dataset"].map(x_pos)
    ys = plot_df["label"].map(y_pos)
    max_genes = plot_df["n_genes"].max()
    sizes = 20 + 180 * (plot_df["n_genes"] / max_genes)

    sc = ax.scatter(xs, ys, s=sizes, c=plot_df["neg_log10_padj"], cmap="viridis",
                     edgecolor="black", linewidth=0.4)
    ax.set_xticks(range(len(dataset_order)))
    ax.set_xticklabels(dataset_order, rotation=45, ha="right")
    ax.set_yticks(range(len(domain_order)))
    ax.set_yticklabels(domain_order, fontsize=8)
    ax.set_xlim(-0.5, len(dataset_order) - 0.5)
    ax.set_ylim(-0.5, len(domain_order) - 0.5)
    ax.set_title(f"PFAM enrichment: {title_stub} (padj_bh < {SIG_ALPHA_PLOT} in >=1 dataset)")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("-log10(padj_bh)")

    for n in sorted(set([plot_df["n_genes"].min(), plot_df["n_genes"].median(), max_genes])):
        ax.scatter([], [], s=20 + 180 * (n / max_genes), c="gray", edgecolor="black",
                   label=f"{int(n)} genes")
    ax.legend(scatterpoints=1, frameon=True, labelspacing=1.5, title="Gene count",
              loc="upper left", bbox_to_anchor=(1.3, 1))

    fig.tight_layout()
    dot_path = os.path.join(OUT_DIR, f"{out_prefix}_dotplot.png")
    fig.savefig(dot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {dot_path}")

    heat_pivot = plot_df.pivot_table(index="label", columns="dataset", values="fold_enrichment")
    heat_pivot = heat_pivot.reindex(index=domain_order, columns=dataset_order)

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("lightgray")

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(heat_pivot.values, cmap=cmap, aspect="auto")
    ax.set_xticks(range(len(dataset_order)))
    ax.set_xticklabels(dataset_order, rotation=45, ha="right")
    ax.set_yticks(range(len(domain_order)))
    ax.set_yticklabels(domain_order, fontsize=8)
    ax.set_title(f"PFAM fold enrichment: {title_stub} (padj_bh < {SIG_ALPHA_PLOT} in >=1 dataset; "
                 f"gray = domain absent from that dataset's gene set)")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Fold enrichment")

    sig_lookup = plot_df.set_index(["label", "dataset"])["padj_bh"]
    for yi, label in enumerate(domain_order):
        for xi, dataset in enumerate(dataset_order):
            p = sig_lookup.get((label, dataset))
            if p is not None and p < SIG_ALPHA_PLOT:
                ax.text(xi, yi, "*", ha="center", va="center", color="white",
                        fontsize=10, fontweight="bold")

    fig.tight_layout()
    heat_path = os.path.join(OUT_DIR, f"{out_prefix}_heatmap.png")
    fig.savefig(heat_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {heat_path}")


def find_seacr_file(seacr_dir, prefix):
    matches = glob.glob(f"{seacr_dir}/{prefix}_*_stringent.stringent.bed")
    if len(matches) != 1:
        raise ValueError(
            f"Expected 1 SEACR file for '{prefix}' in {seacr_dir}, found {len(matches)}: {matches}"
        )
    return matches[0]


def load_peaks_bed(path):
    return pd.read_csv(
        path, sep="\t", header=None,
        names=["chr", "start", "end", "auc", "max", "region"],
    )


print(f"Loading GFF3: {GFF_PATH}")
gff = pd.read_csv(
    GFF_PATH, sep="\t", comment="#", header=None,
    names=["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"],
)
genes_all = gff[gff["feature"] == "gene"].copy()
genes_all["gene_id"] = genes_all["attributes"].apply(get_gene_id)
print(f"Total genes in GFF (genome-wide): {len(genes_all)}")

genes_by_chrom = {c: sub for c, sub in genes_all.groupby("chrom")}


def genes_overlapping_peak(chrom, start, end, flank=GENE_FLANK):
    sub = genes_by_chrom.get(chrom)
    if sub is None:
        return []
    hits = sub[(sub["start"] - flank <= end) & (sub["end"] + flank >= start)]
    return hits["gene_id"].tolist()


print(f"Loading emapper background: {EMAPPER_PATH}")
emapper = load_background(EMAPPER_PATH)
emapper = emapper.rename(columns={"query": "gene_id"})

pfam_to_desc = {}
for _, row in emapper.iterrows():
    if pd.isna(row["PFAMs"]) or row["PFAMs"] == "-":
        continue
    preferred = row["Preferred_name"] if pd.notna(row["Preferred_name"]) and row["Preferred_name"] != "-" else ""
    desc = row["Description"] if pd.notna(row["Description"]) and row["Description"] != "-" else ""
    label = f"{preferred} | {desc}" if preferred else desc
    for domain in str(row["PFAMs"]).split(","):
        domain = domain.strip()
        if domain and domain not in pfam_to_desc:
            pfam_to_desc[domain] = label

bg_long = explode_terms_pfam(emapper)
background_term_counts = bg_long.groupby("term")["gene"].nunique()
background_total = emapper["gene_id"].nunique()
print(f"Background genome: {background_total} genes, {len(background_term_counts)} distinct PFAM domains")

def run_peaks_pfam(seacr_dir, out_prefix, title_stub):
    print("\n" + "=" * 60)
    print(f"PEAKS: {title_stub} ({seacr_dir})")
    print("=" * 60)

    pfam_rows = []
    gene_sheets = {}

    for prefix in PEAK_PREFIXES:
        print(f"\n=== {prefix} ===")
        try:
            path = find_seacr_file(seacr_dir, prefix)
        except ValueError as e:
            print(f"  WARNING: {e} -- skipping.")
            continue

        peaks = load_peaks_bed(path)
        print(f"  {len(peaks)} peaks")

        gene_ids = set()
        for _, row in peaks.iterrows():
            gene_ids.update(genes_overlapping_peak(row["chr"], row["start"], row["end"]))
        print(f"  {len(gene_ids)} unique genes within {GENE_FLANK}bp of a peak")

        annot, stats_df = pfam_stats_for_geneset(
            gene_ids, genes_all, emapper, background_term_counts, background_total, pfam_to_desc, prefix
        )
        gene_sheets[f"{prefix}_genes"] = annot
        if stats_df is not None:
            pfam_rows.append(stats_df)

    if pfam_rows:
        pfam_all = pd.concat(pfam_rows, ignore_index=True).sort_values(["dataset", "padj_bh"])
        out_path = os.path.join(OUT_DIR, f"{out_prefix}_PFAM_by_dataset.xlsx")
        with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
            pfam_all.to_excel(writer, sheet_name="PFAM_all", index=False)
            for name, gdf in gene_sheets.items():
                gdf.to_excel(writer, sheet_name=name[:31], index=False)
        print(f"\nSaved {out_path}")

        make_dotplot_and_heatmap(pfam_all, PEAK_PREFIXES, f"{out_prefix}_PFAM", title_stub)
    else:
        print(f"\nNo PFAM stats computed for any {title_stub} dataset.")


run_peaks_pfam("SEACR", "peaks_unfiltered", "Peaks (unfiltered)")
run_peaks_pfam("SEACR_filtered", "peaks_filtered", "Peaks (FE-filtered)")

print("\nDone.")
