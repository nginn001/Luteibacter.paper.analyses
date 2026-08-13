#!/usr/bin/env python3

import re
import pandas as pd

INPUT = "/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/phylogeny/roary_nosplit/gene_presence_absence.csv"
OUTPUT = "/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/phylogeny/roary_nosplit/gene_copy_number_matrix.tsv"
OUTPUT_BINARY = "/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/phylogeny/roary_nosplit/gene_presence_absence_binary.tsv"

df = pd.read_csv(INPUT, low_memory=False)
strain_cols = list(df.columns[14:])

def count_copies(cell):
    if pd.isna(cell):
        return 0
    cell = str(cell).strip()
    if cell == "":
        return 0
    # Split on semicolons OR any whitespace
    parts = [x for x in re.split(r"[;\s]+", cell) if x]
    return len(parts)

copy_df = pd.DataFrame({"Gene": df["Gene"]})

for strain in strain_cols:
    copy_df[strain] = df[strain].apply(count_copies)

binary_df = copy_df.copy()
binary_df.iloc[:, 1:] = (binary_df.iloc[:, 1:] > 0).astype(int)

copy_df.to_csv(OUTPUT, sep="\t", index=False)
binary_df.to_csv(OUTPUT_BINARY, sep="\t", index=False)

print(f"Saved copy number matrix to: {OUTPUT}")
print(f"Saved binary matrix to:      {OUTPUT_BINARY}")

all_values = copy_df.iloc[:, 1:].to_numpy().ravel()
value_counts = pd.Series(all_values).value_counts().sort_index()

print("\nValue counts across all strain-gene cells:")
print(value_counts)

multicopy_mask = (copy_df.iloc[:, 1:] > 1).any(axis=1)
n_multicopy_genes = multicopy_mask.sum()

print(f"\nGenes with copy number >1 in at least one strain: {n_multicopy_genes}")

if n_multicopy_genes > 0:
    print("\nExample multicopy genes:")
    example = copy_df.loc[multicopy_mask].head(10)
    print(example.to_string(index=False))