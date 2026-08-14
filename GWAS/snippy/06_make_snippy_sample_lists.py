#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

ROOT = Path("/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS")
META = ROOT / "structure/metadata/genome_metadata_with_files.tsv"
ASSEMBLY_DIR = ROOT / "assemblies"
OUTDIR = ROOT / "gwas/inputs"

OUTDIR.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(META, sep="\t")

# Build assembly path from strain name
meta["assembly"] = meta["strain"].apply(lambda s: str(ASSEMBLY_DIR / f"{s}.fasta"))

# Keep only strains with existing assemblies
meta = meta[meta["assembly"].apply(lambda x: Path(x).exists())].copy()

# all genomes
all_df = meta[["strain", "assembly"]].copy()

# sp1 only
sp1_df = meta.loc[meta["putative_species"] == "sp1", ["strain", "assembly"]].copy()

all_df.to_csv(OUTDIR / "snippy_all_samples.tsv", sep="\t", index=False)
sp1_df.to_csv(OUTDIR / "snippy_sp1_samples.tsv", sep="\t", index=False)

print("Wrote:")
print(OUTDIR / "snippy_all_samples.tsv")
print(OUTDIR / "snippy_sp1_samples.tsv")
print(f"All genomes: {len(all_df)}")
print(f"sp1 genomes: {len(sp1_df)}")