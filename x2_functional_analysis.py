import os
from collections import Counter
from io import StringIO

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# config

GFF_PATH = "genome/kieran_braker3_Sratti.gff3"
EMAPPER_PATH = "genome/strongyloides_ratti.emapper.annotations.tsv"
OUT_DIR = "functional_annotation"
CHROM = "SRv7_ChrX2"

CONDITIONS = {
    "Arima_FLF": "fanc/Arima_FLF_eigenvector.bed",
    "Arima_PF": "fanc/Arima_PF_eigenvector.bed",
    "Dovetail_FLF": "fanc/Dovetail_FLF_eigenvector.bed",
    "Dovetail_PF": "fanc/Dovetail_PF_eigenvector.bed",
}

ALTERNATIVE = "greater"

os.makedirs(OUT_DIR, exist_ok=True)

# load

def load_eigenvector(path):
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "score", "strand"],
    )
    df["compartment"] = df["score"].apply(
        lambda x: "A" if x > 0 else ("B" if x < 0 else "NA")
    )
    return df


def get_gene_id(attrs):
    for part in attrs.split(";"):
        if part.startswith("ID="):
            return part.replace("ID=", "")
    return None


def in_A_compartment(start, end, a_bins):
    return ((a_bins["start"] < end) & (a_bins["end"] > start)).any()


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

print(f"Loading GFF3: {GFF_PATH}")
gff = pd.read_csv(
    GFF_PATH,
    sep="\t",
    comment="#",
    header=None,
    names=["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"],
)
genes_all = gff[(gff["feature"] == "gene") & (gff["chrom"] == CHROM)].copy()
genes_all["gene_id"] = genes_all["attributes"].apply(get_gene_id)
print(f"{CHROM} genes in GFF: {len(genes_all)}")

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

#main

all_pfam_rows = []
gene_sheets = {}

for condition, eigenvector_path in CONDITIONS.items():
    print(f"\n=== {condition} ===")
    balanced = load_eigenvector(eigenvector_path)
    a_bins = balanced[(balanced["chrom"] == CHROM) & (balanced["compartment"] == "A")]
    print(f"  {CHROM} A-compartment bins: {len(a_bins)}")

    genes = genes_all.copy()
    genes["in_A_compartment"] = genes.apply(
        lambda r: in_A_compartment(r["start"], r["end"], a_bins), axis=1
    )
    a_genes = genes[genes["in_A_compartment"]].copy()
    print(f"  Genes in {CHROM} A compartment: {len(a_genes)}")

    annot = a_genes.merge(emapper, on="gene_id", how="left")
    n_annotated = annot["Description"].notna().sum()
    pct = round(n_annotated / len(annot) * 100, 1) if len(annot) else 0
    print(f"  Annotated: {n_annotated} / {len(annot)} genes ({pct}%)")
    gene_sheets[f"{condition}_genes"] = annot

    pfam_counts = Counter()
    for val in annot["PFAMs"].dropna():
        for domain in val.split(","):
            domain = domain.strip()
            if domain and domain != "-":
                pfam_counts[domain] += 1

    if not pfam_counts:
        print("  No PFAM domains found in this gene set -- skipping stats.")
        continue

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

    cond_df = pd.DataFrame({
        "dataset": condition,
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
    all_pfam_rows.append(cond_df)

pfam_all = pd.concat(all_pfam_rows, ignore_index=True).sort_values(["dataset", "padj_bh"])

#output
out_path = os.path.join(OUT_DIR, "ChrX2_A_compartment_PFAM_by_dataset.xlsx")
with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
    pfam_all.to_excel(writer, sheet_name="PFAM_all", index=False)
    for name, df in gene_sheets.items():
        df.to_excel(writer, sheet_name=name[:31], index=False)

print(f"\nSaved {out_path}")


import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SIG_ALPHA_PLOT = 0.05

sig_domains = pfam_all.loc[pfam_all["padj_bh"] < SIG_ALPHA_PLOT, "PFAM"].unique()
plot_df = pfam_all[pfam_all["PFAM"].isin(sig_domains)].copy()

if plot_df.empty:
    print(f"\nNo PFAM domains significant at padj_bh < {SIG_ALPHA_PLOT} in any dataset -- skipping plots.")
else:
    plot_df["label"] = plot_df["PFAM"]
    dataset_order = list(CONDITIONS.keys())
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
    ax.set_title(f"PFAM enrichment (padj_bh < {SIG_ALPHA_PLOT} in >=1 dataset)")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("-log10(padj_bh)")

    for n in sorted(set([plot_df["n_genes"].min(), plot_df["n_genes"].median(), max_genes])):
        ax.scatter([], [], s=20 + 180 * (n / max_genes), c="gray", edgecolor="black",
                   label=f"{int(n)} genes")
    ax.legend(scatterpoints=1, frameon=True, labelspacing=1.5, title="Gene count",
              loc="upper left", bbox_to_anchor=(1.3, 1))

    fig.tight_layout()
    dot_path = os.path.join(OUT_DIR, "ChrX2_A_compartment_PFAM_dotplot.png")
    fig.savefig(dot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {dot_path}")

    heat_pivot = plot_df.pivot_table(index="label", columns="dataset", values="fold_enrichment")
    heat_pivot = heat_pivot.reindex(index=domain_order, columns=dataset_order)

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("lightgray")  # domain not found at all in that dataset's gene set

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(heat_pivot.values, cmap=cmap, aspect="auto")
    ax.set_xticks(range(len(dataset_order)))
    ax.set_xticklabels(dataset_order, rotation=45, ha="right")
    ax.set_yticks(range(len(domain_order)))
    ax.set_yticklabels(domain_order, fontsize=8)
    ax.set_title(f"PFAM fold enrichment (padj_bh < {SIG_ALPHA_PLOT} in >=1 dataset; "
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
    heat_path = os.path.join(OUT_DIR, "ChrX2_A_compartment_PFAM_heatmap.png")
    fig.savefig(heat_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {heat_path}")

BUBBLE_PADJ_MAX = 0.1
SIG_LINE = 0.05

qualifying = pfam_all[pfam_all["padj_bh"] < BUBBLE_PADJ_MAX].copy()

if qualifying.empty:
    print(f"\nNo PFAM domains with padj_bh < {BUBBLE_PADJ_MAX} in any dataset -- skipping bubble plot.")
else:
    best_idx = qualifying.groupby("PFAM")["padj_bh"].idxmin()
    best = qualifying.loc[best_idx].copy()

    n_sharing = qualifying.groupby("PFAM")["dataset"].nunique()
    best["n_datasets"] = best["PFAM"].map(n_sharing)
    best["label"] = best["PFAM"]

    fig, ax = plt.subplots(figsize=(18, 9))
    fig.patch.set_facecolor("#e6e6e6")
    ax.set_facecolor("white")

    max_share = best["n_datasets"].max()
    sizes = 60 + 260 * (best["n_datasets"] / max_share)

    ax.scatter(best["padj_bh"], best["fold_enrichment"], s=sizes,
               c="steelblue", alpha=0.55, edgecolor="black", linewidth=0.8, zorder=3)

    for _, r in best.iterrows():
        ax.annotate(r["label"], (r["padj_bh"], r["fold_enrichment"]),
                    fontsize=7, ha="center", va="center", zorder=4)

    ax.axvline(SIG_LINE, color="black", linestyle=":", linewidth=1.2, zorder=2)
    ax.text(SIG_LINE, ax.get_ylim()[1], "padj = 0.05 ", fontsize=8, va="bottom", ha="right")

    ax.set_xlim(BUBBLE_PADJ_MAX, 0)
    ax.set_xlabel("padj (BH-adjusted)")
    ax.set_ylabel("Fold enrichment")
    ax.set_title(f"PFAM enrichment across datasets (padj_bh < {BUBBLE_PADJ_MAX})",
                 fontsize=13, fontweight="bold")
    ax.text(0.5, 1.04, "Bubble size: number of datasets sharing the domain",
            transform=ax.transAxes, ha="center", fontsize=9)

    for n in sorted(best["n_datasets"].unique()):
        ax.scatter([], [], s=60 + 260 * (n / max_share), c="steelblue", alpha=0.55,
                   edgecolor="black", linewidth=0.8, label=f"{n} dataset{'s' if n != 1 else ''}")
    ax.legend(scatterpoints=1, frameon=True, labelspacing=1.3, title="Shared in",
              loc="upper left", bbox_to_anchor=(1.02, 1))

    fig.tight_layout()
    bubble_path = os.path.join(OUT_DIR, "ChrX2_A_compartment_PFAM_bubble.png")
    fig.savefig(bubble_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"Saved {bubble_path}")

print("Done.")
