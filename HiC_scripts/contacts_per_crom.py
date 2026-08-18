import pandas as pd

Q_THRESHOLD = 0.01
OUT_PATH = "contacts_per_chrom.xlsx"

CONTACT_SOURCES = {
    "ArimaFLF": "fithic/ArimaFLF/ArimaFLF_10kb.spline_pass2.res10000.significances.txt.gz",
    "ArimaPF": "fithic/ArimaPF/ArimaPF_10kb.spline_pass2.res10000.significances.txt.gz",
    "DovetailFLF": "fithic/DovetailFLF/DovetailFLF_10kb.spline_pass2.res10000.significances.txt.gz",
    "DovetailPF": "fithic/DovetailPF/DovetailPF_10kb.spline_pass2.res10000.significances.txt.gz",
}

summary_rows = []
sheets = {}

for name, path in CONTACT_SOURCES.items():
    df = pd.read_csv(path, sep="\t")
    sig = df[df["q-value"] < Q_THRESHOLD]
    counts = sig["chr1"].value_counts().sort_index().reset_index()
    counts.columns = ["chromosome", "n_contacts"]
    counts["source"] = name
    counts["total_contacts"] = len(sig)
    counts["pct_of_total"] = (counts["n_contacts"] / len(sig) * 100).round(1)
    sheets[name] = counts
    summary_rows.append(counts)
    print(f"{name}: {len(sig)} significant contacts (q<{Q_THRESHOLD})")

summary_df = pd.concat(summary_rows, ignore_index=True)
cols = ["source", "chromosome", "n_contacts", "pct_of_total", "total_contacts"]
summary_df = summary_df[cols]

with pd.ExcelWriter(OUT_PATH, engine="openpyxl") as writer:
    summary_df.to_excel(writer, sheet_name="Summary", index=False)
    for name, df in sheets.items():
        df[cols].to_excel(writer, sheet_name=name, index=False)

print(f"\nSaved {OUT_PATH}")
