#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import subprocess

ROOT = Path("/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS")

META = ROOT / "structure/metadata/genome_metadata_with_files.tsv"
CLONES = ROOT / "gwas/inputs/clones_to_exclude.txt"
MASH = ROOT / "structure/mash/mash_dist_matrix.tsv"

RTAB = Path("/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS/structure/phylogeny/roary/gene_presence_absence.Rtab")

VCF_ALL = ROOT / "gwas/vcf/all_genomes_knz12_1b/core.vcf.gz"
VCF_SP1 = ROOT / "gwas/vcf/sp1_knz12_1b/core.vcf.gz"
VCF_SVR = ROOT / "gwas/vcf/SVR_genomes_svr3_8d/core.vcf.gz"

OUTROOT = ROOT / "gwas/runs"
OUTROOT.mkdir(parents=True, exist_ok=True)

def read_clone_list(path: Path) -> set[str]:
    if not path.exists():
        return set()
    with open(path) as f:
        return {x.strip() for x in f if x.strip()}

def get_vcf_samples(vcf_file: Path) -> list[str]:
    result = subprocess.run(
        ["bcftools", "query", "-l", str(vcf_file)],
        capture_output=True,
        text=True,
        check=True
    )
    return [x.strip() for x in result.stdout.splitlines() if x.strip()]

def subset_rtab(rtab_file: Path, strains: list[str], outfile: Path) -> None:
    rtab = pd.read_csv(rtab_file, sep="\t", low_memory=False)
    first_col = rtab.columns[0]
    keep_cols = [first_col] + [s for s in strains if s in rtab.columns]
    rtab = rtab[keep_cols]
    rtab.to_csv(outfile, sep="\t", index=False)

def subset_mash(mash_file: Path, strains: list[str], outfile: Path) -> None:
    mash = pd.read_csv(mash_file, sep="\t", index_col=0)
    keep = [s for s in strains if s in mash.index and s in mash.columns]
    mash = mash.loc[keep, keep]
    mash.to_csv(outfile, sep="\t")

def subset_vcf(vcf_file: Path, strains_file: Path, outfile: Path) -> None:
    subprocess.run(
        [
            "bcftools", "view",
            "-S", str(strains_file),
            "-Oz",
            "-o", str(outfile),
            str(vcf_file)
        ],
        check=True
    )
    subprocess.run(["bcftools", "index", "-f", str(outfile)], check=True)

def main():
    meta = pd.read_csv(META, sep="\t")
    meta.columns = meta.columns.str.strip()

    required_cols = {
        "strain", "site", "putative_species",
        "log10_relative_growth_nacl", "log10_relative_growth_peg"
    }
    missing = required_cols - set(meta.columns)
    if missing:
        raise ValueError(f"Missing required metadata columns: {missing}")

    clones = read_clone_list(CLONES)

    subsets = {
        "01_all": {
            "mask": meta["strain"].notna(),
            "vcf": VCF_ALL,
        },
        "02_all_noclone": {
            "mask": meta["strain"].notna() & (~meta["strain"].isin(clones)),
            "vcf": VCF_ALL,
        },
        "03_sp1": {
            "mask": meta["putative_species"] == "sp1",
            "vcf": VCF_SP1,
        },
        "04_sp1_noclone": {
            "mask": (meta["putative_species"] == "sp1") & (~meta["strain"].isin(clones)),
            "vcf": VCF_SP1,
        },
        "05_SVR": {
            "mask": meta["site"] == "SVR",
            "vcf": VCF_SVR,
        },
        "06_SVR_noclone": {
            "mask": (meta["site"] == "SVR") & (~meta["strain"].isin(clones)),
            "vcf": VCF_SVR,
        },
    }

    phenotypes = {
        "nacl": "log10_relative_growth_nacl",
        "peg": "log10_relative_growth_peg",
    }

    for subset_name, subset_info in subsets.items():
        subset_meta = meta.loc[subset_info["mask"]].copy()

        # Restrict to strains actually present in the VCF header
        vcf_samples = set(get_vcf_samples(subset_info["vcf"]))
        subset_meta = subset_meta.loc[subset_meta["strain"].isin(vcf_samples)].copy()

        for pheno_short, pheno_col in phenotypes.items():
            run_dir = OUTROOT / f"{subset_name}_{pheno_short}"
            run_dir.mkdir(parents=True, exist_ok=True)

            run_meta = subset_meta.loc[
                subset_meta[pheno_col].notna(), ["strain", pheno_col]
            ].copy()
            run_meta = run_meta.rename(columns={pheno_col: "phenotype"})

            strains = run_meta["strain"].tolist()

            strains_file = run_dir / "strains.txt"
            with open(strains_file, "w") as f:
                for s in strains:
                    f.write(f"{s}\n")

            run_meta.to_csv(run_dir / "phenotypes.tsv", sep="\t", index=False)
            subset_rtab(RTAB, strains, run_dir / "genes.Rtab")
            subset_mash(MASH, strains, run_dir / "mash.tsv")
            subset_vcf(subset_info["vcf"], strains_file, run_dir / "variants.vcf.gz")

            print(f"Prepared {run_dir} with {len(strains)} strains")

if __name__ == "__main__":
    main()