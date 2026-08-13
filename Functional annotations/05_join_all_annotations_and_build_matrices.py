#!/usr/bin/env python3

import argparse
import csv
import sys
from collections import defaultdict


def read_tsv_dict(infile, key_col):
    data = {}
    with open(infile, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if key_col not in reader.fieldnames:
            raise ValueError(f"{key_col} not found in {infile}")
        for row in reader:
            key = row[key_col].strip()
            if key:
                data[key] = row
    return data


def read_roary_rtab(rtab_file):
    """
    Read Roary gene_presence_absence.Rtab
    first column = cluster/gene name
    remaining columns = genomes
    """
    rows = []
    with open(rtab_file, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames
        if not fieldnames:
            raise ValueError(f"Empty file: {rtab_file}")

        cluster_col = fieldnames[0]
        genomes = fieldnames[1:]

        for row in reader:
            out = {"cluster_name": row[cluster_col]}
            for g in genomes:
                val = row[g]
                try:
                    out[g] = int(val)
                except Exception:
                    try:
                        out[g] = int(float(val))
                    except Exception:
                        out[g] = 0
            rows.append(out)

    return genomes, rows


def index_besthit_by_cluster_name(besthit_rows):
    ann = {}
    for row in besthit_rows.values():
        cluster_name = row.get("cluster_name", "").strip()
        if cluster_name:
            ann[cluster_name] = row
    return ann


def join_annotations(roary_rows, besthit_by_cluster_name, brite_by_ko, pfam_by_safe_cluster):
    """
    Build full cluster-level joined rows.
    Join logic:
      roary cluster_name -> besthit row
      besthit ko -> brite row
      besthit safe_cluster_id -> pfam row
    """
    joined = []

    for row in roary_rows:
        cluster_name = row["cluster_name"]
        out = {"cluster_name": cluster_name}

        besthit = besthit_by_cluster_name.get(cluster_name, {})

        # core besthit columns
        out["safe_cluster_id"] = besthit.get("safe_cluster_id", "")
        out["representative_id"] = besthit.get("representative_id", "")
        out["ko"] = besthit.get("ko", "")
        out["trusted_hit"] = besthit.get("trusted_hit", "")
        out["score"] = besthit.get("score", "")
        out["evalue"] = besthit.get("evalue", "")
        out["threshold"] = besthit.get("threshold", "")
        out["kofamscan_definition"] = besthit.get("kofamscan_definition", "")
        out["ko_definition"] = besthit.get("ko_definition", "")
        out["module_ids"] = besthit.get("module_ids", "")
        out["module_names"] = besthit.get("module_names", "")
        out["pathway_ids"] = besthit.get("pathway_ids", "")
        out["pathway_names"] = besthit.get("pathway_names", "")

        # BRITE via KO
        ko = out["ko"]
        brite = brite_by_ko.get(ko, {})
        out["brite_A_id"] = brite.get("brite_A_id", "")
        out["brite_A_name"] = brite.get("brite_A_name", "")
        out["brite_B_id"] = brite.get("brite_B_id", "")
        out["brite_B_name"] = brite.get("brite_B_name", "")
        out["brite_C_id"] = brite.get("brite_C_id", "")
        out["brite_C_name"] = brite.get("brite_C_name", "")
        out["brite_ko_entry"] = brite.get("brite_ko_entry", "")

        # Pfam via safe_cluster_id
        safe_cluster_id = out["safe_cluster_id"]
        pfam = pfam_by_safe_cluster.get(safe_cluster_id, {})
        out["pfam_target_name"] = pfam.get("pfam_target_name", "")
        out["pfam_accession"] = pfam.get("pfam_accession", "")
        out["pfam_full_evalue"] = pfam.get("pfam_full_evalue", "")
        out["pfam_full_score"] = pfam.get("pfam_full_score", "")
        out["pfam_domain_i_evalue"] = pfam.get("pfam_domain_i_evalue", "")
        out["pfam_domain_score"] = pfam.get("pfam_domain_score", "")
        out["pfam_description"] = pfam.get("pfam_description", "")

        # Roary genome columns
        for k, v in row.items():
            if k != "cluster_name":
                out[k] = v

        joined.append(out)

    return joined


def write_joined_table(joined_rows, outfile, genomes):
    fieldnames = [
        "cluster_name",
        "safe_cluster_id",
        "representative_id",
        "ko",
        "trusted_hit",
        "score",
        "evalue",
        "threshold",
        "kofamscan_definition",
        "ko_definition",
        "module_ids",
        "module_names",
        "pathway_ids",
        "pathway_names",
        "brite_A_id",
        "brite_A_name",
        "brite_B_id",
        "brite_B_name",
        "brite_C_id",
        "brite_C_name",
        "brite_ko_entry",
        "pfam_target_name",
        "pfam_accession",
        "pfam_full_evalue",
        "pfam_full_score",
        "pfam_domain_i_evalue",
        "pfam_domain_score",
        "pfam_description",
    ] + genomes

    with open(outfile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(joined_rows)


def collapse_feature_matrix(joined_rows, genomes, feature_col, split_semicolon=False):
    """
    Feature is present in a genome if any cluster carrying that feature is present.
    """
    feature_matrix = defaultdict(lambda: {g: 0 for g in genomes})

    for row in joined_rows:
        raw = row.get(feature_col, "").strip()
        if not raw:
            continue

        if split_semicolon:
            features = [x.strip() for x in raw.split(";") if x.strip()]
        else:
            features = [raw]

        for feature in features:
            for g in genomes:
                if int(row[g]) > 0:
                    feature_matrix[feature][g] = 1

    return feature_matrix


def write_feature_matrix(feature_matrix, genomes, outfile, first_col_name):
    fieldnames = [first_col_name] + genomes
    with open(outfile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for feature in sorted(feature_matrix.keys()):
            row = {first_col_name: feature}
            for g in genomes:
                row[g] = feature_matrix[feature][g]
            writer.writerow(row)


def write_long_feature_map(joined_rows, outfile):
    """
    Optional useful long table: one row per cluster with selected annotation columns.
    """
    fieldnames = [
        "cluster_name",
        "safe_cluster_id",
        "representative_id",
        "ko",
        "trusted_hit",
        "score",
        "evalue",
        "threshold",
        "kofamscan_definition",
        "ko_definition",
        "module_ids",
        "module_names",
        "pathway_ids",
        "pathway_names",
        "brite_A_id",
        "brite_A_name",
        "brite_B_id",
        "brite_B_name",
        "brite_C_id",
        "brite_C_name",
        "brite_ko_entry",
        "pfam_target_name",
        "pfam_accession",
        "pfam_full_evalue",
        "pfam_full_score",
        "pfam_domain_i_evalue",
        "pfam_domain_score",
        "pfam_description",
    ]
    with open(outfile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in joined_rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def main():
    parser = argparse.ArgumentParser(
        description="Join Roary, KOfam/KEGG, BRITE, and Pfam annotations and build functional presence/absence matrices."
    )
    parser.add_argument("--roary_rtab", required=True, help="Roary gene_presence_absence.Rtab")
    parser.add_argument("--besthit_annotations", required=True, help="Best-hit KOfam/KEGG table")
    parser.add_argument("--brite_parsed", required=True, help="Parsed BRITE KO table")
    parser.add_argument("--pfam_besthit", required=True, help="Best-hit Pfam table")

    parser.add_argument("--out_cluster_join", required=True, help="Full joined cluster annotation table")
    parser.add_argument("--out_annotation_long", required=True, help="Long annotation-only table")

    parser.add_argument("--out_ko_matrix", required=True, help="KO presence/absence matrix")
    parser.add_argument("--out_module_matrix", required=True, help="Module presence/absence matrix")
    parser.add_argument("--out_pathway_matrix", required=True, help="Pathway presence/absence matrix")
    parser.add_argument("--out_brite_B_matrix", required=True, help="BRITE B-level presence/absence matrix")
    parser.add_argument("--out_brite_C_matrix", required=True, help="BRITE C-level presence/absence matrix")
    parser.add_argument("--out_pfam_matrix", required=True, help="Pfam presence/absence matrix")

    args = parser.parse_args()

    genomes, roary_rows = read_roary_rtab(args.roary_rtab)
    besthit_rows = read_tsv_dict(args.besthit_annotations, "safe_cluster_id")
    besthit_by_cluster_name = index_besthit_by_cluster_name(besthit_rows)

    brite_by_ko = read_tsv_dict(args.brite_parsed, "KO")
    pfam_by_safe_cluster = read_tsv_dict(args.pfam_besthit, "safe_cluster_id")

    joined_rows = join_annotations(
        roary_rows=roary_rows,
        besthit_by_cluster_name=besthit_by_cluster_name,
        brite_by_ko=brite_by_ko,
        pfam_by_safe_cluster=pfam_by_safe_cluster
    )

    write_joined_table(joined_rows, args.out_cluster_join, genomes)
    write_long_feature_map(joined_rows, args.out_annotation_long)

    ko_matrix = collapse_feature_matrix(joined_rows, genomes, "ko", split_semicolon=False)
    module_matrix = collapse_feature_matrix(joined_rows, genomes, "module_ids", split_semicolon=True)
    pathway_matrix = collapse_feature_matrix(joined_rows, genomes, "pathway_ids", split_semicolon=True)
    brite_B_matrix = collapse_feature_matrix(joined_rows, genomes, "brite_B_name", split_semicolon=False)
    brite_C_matrix = collapse_feature_matrix(joined_rows, genomes, "brite_C_name", split_semicolon=False)
    pfam_matrix = collapse_feature_matrix(joined_rows, genomes, "pfam_accession", split_semicolon=False)

    write_feature_matrix(ko_matrix, genomes, args.out_ko_matrix, "KO")
    write_feature_matrix(module_matrix, genomes, args.out_module_matrix, "module")
    write_feature_matrix(pathway_matrix, genomes, args.out_pathway_matrix, "pathway")
    write_feature_matrix(brite_B_matrix, genomes, args.out_brite_B_matrix, "brite_B")
    write_feature_matrix(brite_C_matrix, genomes, args.out_brite_C_matrix, "brite_C")
    write_feature_matrix(pfam_matrix, genomes, args.out_pfam_matrix, "pfam")

    sys.stderr.write(f"Roary clusters loaded: {len(roary_rows)}\n")
    sys.stderr.write(f"Best-hit annotation rows loaded: {len(besthit_rows)}\n")
    sys.stderr.write(f"BRITE KO rows loaded: {len(brite_by_ko)}\n")
    sys.stderr.write(f"Pfam rows loaded: {len(pfam_by_safe_cluster)}\n")
    sys.stderr.write(f"Joined cluster rows written: {len(joined_rows)}\n")
    sys.stderr.write(f"KO features written: {len(ko_matrix)}\n")
    sys.stderr.write(f"Module features written: {len(module_matrix)}\n")
    sys.stderr.write(f"Pathway features written: {len(pathway_matrix)}\n")
    sys.stderr.write(f"BRITE B features written: {len(brite_B_matrix)}\n")
    sys.stderr.write(f"BRITE C features written: {len(brite_C_matrix)}\n")
    sys.stderr.write(f"Pfam features written: {len(pfam_matrix)}\n")


if __name__ == "__main__":
    main()