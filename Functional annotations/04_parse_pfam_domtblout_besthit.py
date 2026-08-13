#!/usr/bin/env python3

import argparse
import csv
import sys


def parse_domtblout(domtbl_file):
    """
    Parse HMMER domtblout format.
    Keep one best Pfam hit per query sequence based on full-sequence bit score.
    """
    best = {}

    with open(domtbl_file, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.strip().split()
            if len(fields) < 23:
                continue

            target_name = fields[0]   # Pfam accession/name depending on db
            target_accession = fields[1]
            query_name = fields[3]    # should be cluster_000001
            query_accession = fields[4]

            full_evalue = fields[6]
            full_score = fields[7]
            full_bias = fields[8]

            domain_i_evalue = fields[12]
            domain_score = fields[13]

            description = " ".join(fields[22:]) if len(fields) > 22 else ""

            try:
                score_num = float(full_score)
            except ValueError:
                score_num = float("-inf")

            try:
                ie_num = float(domain_i_evalue)
            except ValueError:
                ie_num = float("inf")

            row = {
                "safe_cluster_id": query_name,
                "pfam_target_name": target_name,
                "pfam_accession": target_accession,
                "pfam_full_evalue": full_evalue,
                "pfam_full_score": full_score,
                "pfam_domain_i_evalue": domain_i_evalue,
                "pfam_domain_score": domain_score,
                "pfam_description": description,
                "_score_num": score_num,
                "_ie_num": ie_num
            }

            if query_name not in best:
                best[query_name] = row
            else:
                old = best[query_name]
                new_key = (row["_score_num"], -row["_ie_num"])
                old_key = (old["_score_num"], -old["_ie_num"])
                if new_key > old_key:
                    best[query_name] = row

    return best


def write_best(best, outfile):
    fieldnames = [
        "safe_cluster_id",
        "pfam_target_name",
        "pfam_accession",
        "pfam_full_evalue",
        "pfam_full_score",
        "pfam_domain_i_evalue",
        "pfam_domain_score",
        "pfam_description"
    ]
    with open(outfile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for key in sorted(best):
            row = {k: v for k, v in best[key].items() if not k.startswith("_")}
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description="Parse hmmscan domtblout and keep best Pfam hit per query.")
    parser.add_argument("--domtblout", required=True)
    parser.add_argument("--out_tsv", required=True)
    args = parser.parse_args()

    best = parse_domtblout(args.domtblout)
    write_best(best, args.out_tsv)
    sys.stderr.write(f"Best Pfam rows written: {len(best)}\n")


if __name__ == "__main__":
    main()