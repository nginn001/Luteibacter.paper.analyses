#!/usr/bin/env python3

import argparse
import csv
import re
import sys


def parse_brite_file(brite_file):
    rows = []

    current_A_id = ""
    current_A_name = ""
    current_B_id = ""
    current_B_name = ""
    current_C_id = ""
    current_C_name = ""

    with open(brite_file, "r") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line:
                continue

            if line.startswith("+") or line.startswith("!"):
                continue

            if line.startswith("A"):
                # Example: A09100 Metabolism
                m = re.match(r"^A(\S+)\s+(.*)$", line)
                if m:
                    current_A_id = m.group(1)
                    current_A_name = m.group(2).strip()
                    current_B_id = ""
                    current_B_name = ""
                    current_C_id = ""
                    current_C_name = ""
                continue

            if line.startswith("B"):
                # Example: B  09101 Carbohydrate metabolism
                m = re.match(r"^B\s+(\S+)\s+(.*)$", line)
                if m:
                    current_B_id = m.group(1)
                    current_B_name = m.group(2).strip()
                    current_C_id = ""
                    current_C_name = ""
                continue

            if line.startswith("C"):
                # Example: C    00010 Glycolysis / Gluconeogenesis [PATH:ko00010]
                m = re.match(r"^C\s+(\S+)\s+(.*?)(?:\s+\[PATH:(ko\d+)\])?$", line)
                if m:
                    current_C_id = m.group(3) if m.group(3) else m.group(1)
                    current_C_name = m.group(2).strip()
                continue

            if line.startswith("D"):
                # Example:
                # D      K00844  HK; hexokinase [EC:2.7.1.1]
                m = re.match(r"^D\s+(K\d+)\s+(.*)$", line)
                if m:
                    ko = m.group(1)
                    ko_desc = m.group(2).strip()
                    rows.append({
                        "KO": ko,
                        "brite_A_id": current_A_id,
                        "brite_A_name": current_A_name,
                        "brite_B_id": current_B_id,
                        "brite_B_name": current_B_name,
                        "brite_C_id": current_C_id,
                        "brite_C_name": current_C_name,
                        "brite_ko_entry": ko_desc
                    })
                continue

    return rows


def write_rows(rows, outfile):
    fieldnames = [
        "KO",
        "brite_A_id",
        "brite_A_name",
        "brite_B_id",
        "brite_B_name",
        "brite_C_id",
        "brite_C_name",
        "brite_ko_entry"
    ]
    with open(outfile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Parse KEGG BRITE ko00001 hierarchy text into a KO mapping table.")
    parser.add_argument("--brite_txt", required=True, help="brite_ko00001 text file")
    parser.add_argument("--out_tsv", required=True, help="Output KO->BRITE TSV")
    args = parser.parse_args()

    rows = parse_brite_file(args.brite_txt)
    write_rows(rows, args.out_tsv)
    sys.stderr.write(f"Rows written: {len(rows)}\n")


if __name__ == "__main__":
    main()