#!/usr/bin/env python3
import pandas as pd

infile = "/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/mash/mash_dist.tsv"
outfile = "/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/mash/mash_dist_matrix.tsv"

df = pd.read_csv(
    infile, sep="\t", header=None,
    names=["query", "reference", "distance", "pvalue", "shared_hashes"]
)

def clean_name(x):
    x = str(x).split("/")[-1]
    for suf in [".fasta", ".fa", ".fna"]:
        if x.endswith(suf):
            x = x[:-len(suf)]
    return x

df["query"] = df["query"].map(clean_name)
df["reference"] = df["reference"].map(clean_name)

samples = sorted(set(df["query"]).union(df["reference"]))
mat = pd.DataFrame(0.0, index=samples, columns=samples)

for _, row in df.iterrows():
    mat.loc[row["query"], row["reference"]] = row["distance"]
    mat.loc[row["reference"], row["query"]] = row["distance"]

mat.to_csv(outfile, sep="\t")
print(f"Wrote {outfile}")