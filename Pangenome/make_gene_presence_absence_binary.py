#!/usr/bin/env python3
import pandas as pd

INPUT = "/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/phylogeny/roary/gene_presence_absence.csv"
OUTPUT = "/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/phylogeny/roary/gene_presence_absence_binary.tsv"

df = pd.read_csv(INPUT, low_memory=False)

strain_cols = list(df.columns[14:])

out = pd.DataFrame({"Gene": df["Gene"]})

for strain in strain_cols:
    out[strain] = ((df[strain].notna()) & (df[strain].astype(str) != "")).astype(int)

out.to_csv(OUTPUT, sep="\t", index=False)
print(f"Wrote {OUTPUT}")