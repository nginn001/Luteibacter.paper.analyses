#!/usr/bin/env python3

import argparse
import csv
import sys
from collections import defaultdict


def read_cluster_map(map_file):
    cluster_map = {}
    with open(map_file, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"safe_cluster_id", "cluster_name", "representative_id"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing required columns in {map_file}: {missing}")
        for row in reader:
            cluster_map[row["safe_cluster_id"]] = row
    return cluster_map


def parse_kofamscan_results(kofam_file):
    """
    Parse KOfamScan detail output.
    Handles lines with or without leading '*'.
    """
    hits = []

    with open(kofam_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if not fields:
                continue

            if fields[0] == "*":
                trusted_hit = "yes"
                if len(fields) < 7:
                    continue
                safe_cluster_id = fields[1]
                ko = fields[2]
                threshold = fields[3]
                score = fields[4]
                evalue = fields[5]
                definition = " ".join(fields[6:])
            else:
                trusted_hit = "no"
                if len(fields) < 6:
                    continue
                safe_cluster_id = fields[0]
                ko = fields[1]
                threshold = fields[2]
                score = fields[3]
                evalue = fields[4]
                definition = " ".join(fields[5:])

            try:
                score_num = float(score)
            except ValueError:
                score_num = float("-inf")

            try:
                evalue_num = float(evalue)
            except ValueError:
                evalue_num = float("inf")

            hits.append({
                "safe_cluster_id": safe_cluster_id,
                "ko": ko,
                "threshold": threshold,
                "score": score,
                "score_num": score_num,
                "evalue": evalue,
                "evalue_num": evalue_num,
                "kofamscan_definition": definition,
                "trusted_hit": trusted_hit
            })

    return hits


def select_best_hit_per_cluster(hits):
    """
    Keep the top KO hit per representative cluster,
    whether trusted or not.
    Ranking:
      1. highest score
      2. trusted hit preferred if tied
      3. lowest evalue if still tied
    """
    best = {}

    for hit in hits:
        cid = hit["safe_cluster_id"]
        if cid not in best:
            best[cid] = hit
            continue

        old = best[cid]
        new_key = (
            hit["score_num"],
            1 if hit["trusted_hit"] == "yes" else 0,
            -hit["evalue_num"]
        )
        old_key = (
            old["score_num"],
            1 if old["trusted_hit"] == "yes" else 0,
            -old["evalue_num"]
        )

        if new_key > old_key:
            best[cid] = hit

    return best


def read_two_col_tsv(infile):
    out = {}
    with open(infile, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            key = row[0].strip()
            val = row[1].strip()
            if key:
                out[key] = val
    return out


def strip_prefix(x):
    if ":" in x:
        return x.split(":", 1)[1]
    return x


def read_link_table(infile):
    pairs = []
    with open(infile, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            left = strip_prefix(row[0].strip())
            right = strip_prefix(row[1].strip())
            if left and right:
                pairs.append((left, right))
    return pairs


def invert_links_to_multimap(pairs):
    mm = defaultdict(list)
    for a, b in pairs:
        if b not in mm[a]:
            mm[a].append(b)
    return mm


def join_names(ids, id_to_name):
    vals = []
    for x in ids:
        name = id_to_name.get(x, "")
        if name:
            vals.append(f"{x}: {name}")
        else:
            vals.append(x)
    return "; ".join(vals)


def write_besthit_table(rows, outfile):
    fieldnames = [
        "safe_cluster_id",
        "cluster_name",
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
        "pathway_names"
    ]
    with open(outfile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Create best-hit KOfam/KEGG annotation table from raw kofamscan_results.tsv"
    )
    parser.add_argument("--kofam_results", required=True)
    parser.add_argument("--cluster_map", required=True)
    parser.add_argument("--list_ko", required=True)
    parser.add_argument("--link_module_ko", required=True)
    parser.add_argument("--list_module", required=True)
    parser.add_argument("--link_pathway_ko", required=True)
    parser.add_argument("--list_pathway", required=True)
    parser.add_argument("--out_besthit", required=True)
    args = parser.parse_args()

    cluster_map = read_cluster_map(args.cluster_map)
    hits = parse_kofamscan_results(args.kofam_results)
    best = select_best_hit_per_cluster(hits)

    ko_to_def = read_two_col_tsv(args.list_ko)
    module_pairs = read_link_table(args.link_module_ko)
    pathway_pairs = read_link_table(args.link_pathway_ko)
    module_to_name = read_two_col_tsv(args.list_module)
    pathway_to_name = read_two_col_tsv(args.list_pathway)

    ko_to_modules = invert_links_to_multimap(module_pairs)
    ko_to_pathways = invert_links_to_multimap(pathway_pairs)

    rows = []
    for safe_cluster_id, hit in sorted(best.items()):
        ko = hit["ko"]
        modules = ko_to_modules.get(ko, [])
        pathways = ko_to_pathways.get(ko, [])

        row = {
            "safe_cluster_id": safe_cluster_id,
            "cluster_name": cluster_map.get(safe_cluster_id, {}).get("cluster_name", ""),
            "representative_id": cluster_map.get(safe_cluster_id, {}).get("representative_id", ""),
            "ko": ko,
            "trusted_hit": hit["trusted_hit"],
            "score": hit["score"],
            "evalue": hit["evalue"],
            "threshold": hit["threshold"],
            "kofamscan_definition": hit["kofamscan_definition"],
            "ko_definition": ko_to_def.get(ko, ""),
            "module_ids": "; ".join(modules),
            "module_names": join_names(modules, module_to_name),
            "pathway_ids": "; ".join(pathways),
            "pathway_names": join_names(pathways, pathway_to_name),
        }
        rows.append(row)

    write_besthit_table(rows, args.out_besthit)

    sys.stderr.write(f"Total KOfam hits parsed: {len(hits)}\n")
    sys.stderr.write(f"Best-hit rows written: {len(rows)}\n")


if __name__ == "__main__":
    main()