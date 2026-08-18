import os
from collections import Counter
from io import StringIO

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

matplotlib.use("Agg")
import matplotlib.pyplot as plt

GFF_PATH = "genome/kieran_braker3_Sratti.gff3"
EMAPPER_PATH = "genome/strongyloides_ratti.emapper.annotations.tsv"
OUT_DIR = "functional_annotation"

CONTACT_CONDITIONS = ["ArimaFLF", "ArimaPF", "DovetailFLF", "DovetailPF"]
LOOP_CONDITIONS = ["DovetailFLF", "DovetailPF"]  # Arima dot files exist but are empty

CONTACT_PATH_TEMPLATE = "FitHiC/{cond}/{cond}_10kb.spline_pass2.res10000.significances.txt.gz"
LOOP_PATH_TEMPLATE = "cooltools/{cond}_dots.bedpe"

HALF_WINDOW = 5000
Q_THRESHOLD = 0.01
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
    cmap.set_bad("lightgray")  # domain not found at all in that dataset's gene set

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


print(f"Loading GFF3: {GFF_PATH}")
gff = pd.read_csv(
    GFF_PATH, sep="\t", comment="#", header=None,
    names=["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"],
)
genes_all = gff[gff["feature"] == "gene"].copy()
genes_all["gene_id"] = genes_all["attributes"].apply(get_gene_id)
print(f"Total genes in GFF (genome-wide): {len(genes_all)}")

genes_by_chrom = {c: sub for c, sub in genes_all.groupby("chrom")}


def genes_in_anchor(chrom, start, end):
    sub = genes_by_chrom.get(chrom)
    if sub is None:
        return []
    hits = sub[(sub["start"] <= end) & (sub["end"] >= start)]
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


print("\n" + "=" * 60)
print("CONTACTS")
print("=" * 60)

contact_pfam_rows = []
contact_gene_sheets = {}

for cond in CONTACT_CONDITIONS:
    path = CONTACT_PATH_TEMPLATE.format(cond=cond)
    print(f"\n=== {cond} ===")
    df = pd.read_csv(path, sep="\t")
    sig = df[df["q-value"] < Q_THRESHOLD].copy()
    print(f"  {len(sig)} significant contacts (q<{Q_THRESHOLD}) of {len(df)} total")

    gene_ids = set()
    for _, row in sig.iterrows():
        gene_ids.update(genes_in_anchor(row["chr1"], row["fragmentMid1"] - HALF_WINDOW,
                                         row["fragmentMid1"] + HALF_WINDOW))
        gene_ids.update(genes_in_anchor(row["chr2"], row["fragmentMid2"] - HALF_WINDOW,
                                         row["fragmentMid2"] + HALF_WINDOW))
    print(f"  {len(gene_ids)} unique anchor genes")

    annot, stats_df = pfam_stats_for_geneset(
        gene_ids, genes_all, emapper, background_term_counts, background_total, pfam_to_desc, cond
    )
    contact_gene_sheets[f"{cond}_genes"] = annot
    if stats_df is not None:
        contact_pfam_rows.append(stats_df)

if contact_pfam_rows:
    contact_pfam_all = pd.concat(contact_pfam_rows, ignore_index=True).sort_values(["dataset", "padj_bh"])
    contact_out_path = os.path.join(OUT_DIR, "contacts_PFAM_by_dataset.xlsx")
    with pd.ExcelWriter(contact_out_path, engine="openpyxl") as writer:
        contact_pfam_all.to_excel(writer, sheet_name="PFAM_all", index=False)
        for name, gdf in contact_gene_sheets.items():
            gdf.to_excel(writer, sheet_name=name[:31], index=False)
    print(f"\nSaved {contact_out_path}")

    make_dotplot_and_heatmap(contact_pfam_all, CONTACT_CONDITIONS, "contacts_PFAM", "Significant contacts")
else:
    print("\nNo PFAM stats computed for any contact dataset.")


print("\n" + "=" * 60)
print("LOOPS (Dovetail only -- Arima dot-call files exist but are empty)")
print("=" * 60)

loop_pfam_rows = []
loop_gene_sheets = {}

for cond in LOOP_CONDITIONS:
    path = LOOP_PATH_TEMPLATE.format(cond=cond)
    print(f"\n=== {cond} ===")
    loops = pd.read_csv(path, sep="\t")
    print(f"  {len(loops)} loops")

    gene_ids = set()
    for _, row in loops.iterrows():
        gene_ids.update(genes_in_anchor(row["chrom1"], row["start1"], row["end1"]))
        gene_ids.update(genes_in_anchor(row["chrom2"], row["start2"], row["end2"]))
    print(f"  {len(gene_ids)} unique anchor genes")

    annot, stats_df = pfam_stats_for_geneset(
        gene_ids, genes_all, emapper, background_term_counts, background_total, pfam_to_desc, cond
    )
    loop_gene_sheets[f"{cond}_genes"] = annot
    if stats_df is not None:
        loop_pfam_rows.append(stats_df)

if loop_pfam_rows:
    loop_pfam_all = pd.concat(loop_pfam_rows, ignore_index=True).sort_values(["dataset", "padj_bh"])
    loop_out_path = os.path.join(OUT_DIR, "loop_anchor_PFAM_by_dataset.xlsx")
    with pd.ExcelWriter(loop_out_path, engine="openpyxl") as writer:
        loop_pfam_all.to_excel(writer, sheet_name="PFAM_all", index=False)
        for name, gdf in loop_gene_sheets.items():
            gdf.to_excel(writer, sheet_name=name[:31], index=False)
    print(f"\nSaved {loop_out_path}")

    make_dotplot_and_heatmap(loop_pfam_all, LOOP_CONDITIONS, "loop_anchor_PFAM", "Loop anchors (Dovetail)")
else:
    print("\nNo PFAM stats computed for any loop dataset.")

print("\nDone.")
