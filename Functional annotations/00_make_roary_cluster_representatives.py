#!/usr/bin/env python3

import argparse
import glob
import os
import sys
from collections import defaultdict

from Bio import SeqIO


def load_faa_sequences(faa_dir):
    """
    Load all protein sequences from Bakta .faa files.

    Returns
    -------
    seqs : dict
        key = locus tag / sequence ID (first token of FASTA header)
        value = SeqRecord
    source_file : dict
        key = locus tag
        value = source .faa filename
    duplicates : dict
        key = locus tag
        value = list of files where duplicate IDs were found
    """
    seqs = {}
    source_file = {}
    duplicates = defaultdict(list)

    faa_files = sorted(glob.glob(os.path.join(faa_dir, "*.faa")))
    if not faa_files:
        raise FileNotFoundError(f"No .faa files found in: {faa_dir}")

    for faa in faa_files:
        basename = os.path.basename(faa)
        for rec in SeqIO.parse(faa, "fasta"):
            seq_id = rec.id.split()[0]

            if seq_id in seqs:
                duplicates[seq_id].append(basename)
                continue

            seqs[seq_id] = rec
            source_file[seq_id] = basename

    return seqs, source_file, duplicates


def parse_clustered_proteins(cluster_file):
    """
    Parse Roary clustered_proteins file.

    Expected format:
    cluster_name: locus1<TAB>locus2<TAB>locus3...

    Returns
    -------
    clusters : list of tuples
        [(cluster_name, [member1, member2, ...]), ...]
    """
    clusters = []

    with open(cluster_file, "r") as fh:
        for line_num, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line.strip():
                continue

            if ":" not in line:
                raise ValueError(
                    f"Line {line_num} in {cluster_file} has no ':' separator:\n{line}"
                )

            cluster_name, members_str = line.split(":", 1)
            cluster_name = cluster_name.strip()

            members = [x.strip() for x in members_str.strip().split("\t") if x.strip()]
            if not members:
                raise ValueError(
                    f"Line {line_num} in {cluster_file} has no member IDs:\n{line}"
                )

            clusters.append((cluster_name, members))

    return clusters


def choose_representative(members_present, seqs, method="first"):
    """
    Choose representative protein from a list of member IDs present in seqs.
    """
    if method == "first":
        return members_present[0]
    elif method == "longest":
        return max(members_present, key=lambda x: len(seqs[x].seq))
    else:
        raise ValueError(f"Unknown representative method: {method}")


def main():
    parser = argparse.ArgumentParser(
        description="Create one representative protein FASTA per Roary cluster using Bakta .faa files."
    )
    parser.add_argument(
        "--clustered_proteins",
        required=True,
        help="Path to Roary clustered_proteins file"
    )
    parser.add_argument(
        "--faa_dir",
        required=True,
        help="Directory containing Bakta .faa protein files"
    )
    parser.add_argument(
        "--out_faa",
        required=True,
        help="Output representative protein FASTA"
    )
    parser.add_argument(
        "--out_map",
        required=True,
        help="Output TSV mapping cluster -> representative"
    )
    parser.add_argument(
        "--out_missing",
        required=True,
        help="Output TSV listing cluster members not found in .faa files"
    )
    parser.add_argument(
        "--pick",
        choices=["first", "longest"],
        default="first",
        help="How to choose the representative sequence per cluster (default: first)"
    )

    args = parser.parse_args()

    seqs, source_file, duplicates = load_faa_sequences(args.faa_dir)
    clusters = parse_clustered_proteins(args.clustered_proteins)

    n_clusters = 0
    n_written = 0
    n_all_missing = 0
    total_missing_members = 0

    with open(args.out_faa, "w") as faa_out, \
         open(args.out_map, "w") as map_out, \
         open(args.out_missing, "w") as miss_out:

        map_out.write(
            "safe_cluster_id\tcluster_name\trepresentative_id\trepresentative_length_aa\t"
            "n_cluster_members\tn_members_found\tn_members_missing\trepresentative_source_faa\n"
        )
        miss_out.write("cluster_name\tmissing_member_id\n")

        for cluster_name, members in clusters:
            n_clusters += 1

            members_present = [m for m in members if m in seqs]
            members_missing = [m for m in members if m not in seqs]

            total_missing_members += len(members_missing)

            for m in members_missing:
                miss_out.write(f"{cluster_name}\t{m}\n")

            if not members_present:
                n_all_missing += 1
                continue

            rep_id = choose_representative(members_present, seqs, method=args.pick)
            rep_rec = seqs[rep_id]

            safe_cluster_id = f"cluster_{n_clusters:06d}"

            rep_header = (
                f"{safe_cluster_id} "
                f'original_cluster="{cluster_name}" '
                f"rep={rep_id} "
                f"n_members={len(members)} "
                f"n_found={len(members_present)} "
                f"n_missing={len(members_missing)} "
                f"source={source_file[rep_id]}"
            )

            faa_out.write(f">{rep_header}\n{str(rep_rec.seq)}\n")

            map_out.write(
                f"{safe_cluster_id}\t{cluster_name}\t{rep_id}\t{len(rep_rec.seq)}\t"
                f"{len(members)}\t{len(members_present)}\t{len(members_missing)}\t"
                f"{source_file[rep_id]}\n"
            )

            n_written += 1

    sys.stderr.write(f"Loaded protein sequences: {len(seqs)}\n")
    sys.stderr.write(f"Clusters parsed: {n_clusters}\n")
    sys.stderr.write(f"Representative sequences written: {n_written}\n")
    sys.stderr.write(f"Clusters with no members found in .faa files: {n_all_missing}\n")
    sys.stderr.write(f"Total missing cluster-member IDs: {total_missing_members}\n")

    if duplicates:
        sys.stderr.write(
            f"WARNING: {len(duplicates)} duplicate sequence IDs were found across .faa files.\n"
            "The first occurrence was kept. Duplicate IDs may indicate a problem if unexpected.\n"
        )


if __name__ == "__main__":
    main()