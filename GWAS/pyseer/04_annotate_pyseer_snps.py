#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re

ROOT = Path("/bigdata/roperlab/nginn001/KU_Luteibacter/microGWAS")
GWAS_RUNS = ROOT / "gwas" / "runs"

KNZ_REF_GFF = Path("/bigdata/roperlab/nginn001/KU_Luteibacter/gffs.bakta.files/KNZ12-1B.gff3")
SVR_REF_GFF = Path("/bigdata/roperlab/nginn001/KU_Luteibacter/gffs.bakta.files/SVR3-8D.gff3")

RUNS_KNZ_REF = {
    "01_all_nacl", "01_all_peg",
    "02_all_noclone_nacl", "02_all_noclone_peg",
    "03_sp1_nacl", "03_sp1_peg",
    "04_sp1_noclone_nacl", "04_sp1_noclone_peg",
}

RUNS_SVR_REF = {
    "05_SVR_nacl", "05_SVR_peg",
    "06_SVR_noclone_nacl", "06_SVR_noclone_peg",
}


def parse_variant_name(v: str):
    """
    Parse pyseer SNP variant names like:
    contig_1_1144_C_G
    """
    m = re.match(r"(.+)_(\d+)_([A-Za-z]+)_([A-Za-z]+)$", str(v))
    if not m:
        return None, None, None, None
    return m.group(1), int(m.group(2)), m.group(3), m.group(4)


def parse_gff_attributes(attr_str: str):
    attrs = {}
    if pd.isna(attr_str):
        return attrs
    for item in str(attr_str).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k] = v
    return attrs


def load_gff_features(gff_path: Path):
    rows = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts

            # prioritize CDS and gene only
            if feature_type not in {"CDS", "gene"}:
                continue

            attrs = parse_gff_attributes(attributes)

            rows.append({
                "contig": seqid,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "feature_type": feature_type,
                "locus_tag": attrs.get("locus_tag"),
                "gene_name": attrs.get("gene"),
                "product": attrs.get("product"),
            })

    gff_df = pd.DataFrame(rows)

    if not gff_df.empty:
        # Prefer CDS over gene if both overlap
        feature_priority = {"CDS": 0, "gene": 1}
        gff_df["feature_priority"] = gff_df["feature_type"].map(feature_priority).fillna(99)
        gff_df = gff_df.sort_values(
            ["contig", "start", "end", "feature_priority"]
        ).reset_index(drop=True)

    return gff_df


def choose_reference_gff(run_name: str):
    if run_name in RUNS_KNZ_REF:
        return KNZ_REF_GFF, "KNZ12-1B"
    elif run_name in RUNS_SVR_REF:
        return SVR_REF_GFF, "SVR3-8D"
    else:
        raise ValueError(f"Run name not recognized for reference selection: {run_name}")


def annotate_snps_with_gff(snps_df: pd.DataFrame, gff_df: pd.DataFrame):
    annotated_rows = []

    gff_by_contig = {
        c: d.reset_index(drop=True)
        for c, d in gff_df.groupby("contig")
    }

    for _, row in snps_df.iterrows():
        contig = row["contig"]
        pos = row["pos"]

        if contig not in gff_by_contig:
            annotated_rows.append({
                **row.to_dict(),
                "in_gene": False,
                "feature_type": None,
                "locus_tag": None,
                "gene_name": None,
                "product": None,
                "strand": None,
            })
            continue

        sub = gff_by_contig[contig]
        hits = sub[(sub["start"] <= pos) & (sub["end"] >= pos)]

        if hits.empty:
            annotated_rows.append({
                **row.to_dict(),
                "in_gene": False,
                "feature_type": None,
                "locus_tag": None,
                "gene_name": None,
                "product": None,
                "strand": None,
            })
        else:
            best = hits.iloc[0]
            annotated_rows.append({
                **row.to_dict(),
                "in_gene": True,
                "feature_type": best["feature_type"],
                "locus_tag": best["locus_tag"],
                "gene_name": best["gene_name"],
                "product": best["product"],
                "strand": best["strand"],
            })

    return pd.DataFrame(annotated_rows)


def process_run(run_dir: Path):
    snps_file = run_dir / "snps_pyseer.tsv"
    if not snps_file.exists():
        print(f"Skipping {run_dir.name}: no snps_pyseer.tsv")
        return

    gff_path, ref_strain = choose_reference_gff(run_dir.name)
    if not gff_path.exists():
        raise FileNotFoundError(f"Missing GFF for {run_dir.name}: {gff_path}")

    print(f"Processing {run_dir.name}")
    print(f"Reference strain: {ref_strain}")
    print(f"GFF: {gff_path}")

    snps = pd.read_csv(snps_file, sep="\t", low_memory=False)

    if "variant" not in snps.columns:
        raise ValueError(f"{snps_file} is missing a 'variant' column")

    parsed = snps["variant"].apply(parse_variant_name)
    snps["contig"] = [x[0] for x in parsed]
    snps["pos"] = [x[1] for x in parsed]
    snps["ref"] = [x[2] for x in parsed]
    snps["alt"] = [x[3] for x in parsed]

    gff_df = load_gff_features(gff_path)
    annot = annotate_snps_with_gff(snps, gff_df)

    annot["reference_strain"] = ref_strain
    annot["reference_gff"] = str(gff_path)

    out_file = run_dir / "snps_pyseer_annotated.tsv"
    annot.to_csv(out_file, sep="\t", index=False)

    print(f"Wrote {out_file}")


def main():
    run_dirs = sorted([p for p in GWAS_RUNS.iterdir() if p.is_dir()])

    for run_dir in run_dirs:
        if run_dir.name.startswith("."):
            continue
        process_run(run_dir)


if __name__ == "__main__":
    main()