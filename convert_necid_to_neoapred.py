"""Convert NECID peptide-prediction results to NeoaPred input CSV."""

import argparse
import csv
import re
from collections import OrderedDict
from pathlib import Path


INPUT_COLUMNS = (
    "necid_row_idx",
    "peptide_orign",
    "alleleA",
    "peptide_sequence",
)
ALLELE_PATTERN = re.compile(r"HLA-[ABC]\*(\d{2}):(\d{2})")


def convert(input_path, output_path):
    """Write paired WT and mutant peptides in NeoaPred input format."""
    pairs = OrderedDict()

    with open(input_path, newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source)
        missing_columns = set(INPUT_COLUMNS) - set(reader.fieldnames or ())
        if missing_columns:
            raise ValueError(
                "Input is missing required columns: {}".format(
                    ", ".join(sorted(missing_columns))
                )
            )

        for line_number, row in enumerate(reader, start=2):
            row_id = row["necid_row_idx"]
            allele = row["alleleA"]
            origin = row["peptide_orign"]
            sequence = row["peptide_sequence"]

            if not all((row_id, allele, origin, sequence)):
                raise ValueError(
                    "Missing required value on input line {}".format(line_number)
                )

            allele_match = ALLELE_PATTERN.fullmatch(allele)
            if allele_match is None:
                raise ValueError(
                    "Invalid allele {!r} on input line {}".format(
                        allele, line_number
                    )
                )
            if origin not in ("wt_peptide", "mut_peptide"):
                raise ValueError(
                    "Invalid peptide_orign {!r} on input line {}".format(
                        origin, line_number
                    )
                )

            key = (row_id, allele)
            peptide_pair = pairs.setdefault(key, {})
            if origin in peptide_pair:
                raise ValueError(
                    "Duplicate {} for ID {} and allele {} on input line {}".format(
                        origin, row_id, allele, line_number
                    )
                )
            peptide_pair[origin] = sequence

    output_rows = []
    for (row_id, allele), peptide_pair in pairs.items():
        missing_peptides = {"wt_peptide", "mut_peptide"} - set(peptide_pair)
        if missing_peptides:
            raise ValueError(
                "ID {} and allele {} are missing {}".format(
                    row_id, allele, " and ".join(sorted(missing_peptides))
                )
            )

        allele_match = ALLELE_PATTERN.fullmatch(allele)
        output_rows.append(
            (
                row_id,
                "A{}{}".format(*allele_match.groups()),
                peptide_pair["wt_peptide"],
                peptide_pair["mut_peptide"],
            )
        )

    with open(output_path, "w", newline="", encoding="utf-8") as destination:
        writer = csv.writer(destination)
        writer.writerow(("ID", "Allele", "WT", "Mut"))
        writer.writerows(output_rows)


def main():
    parser = argparse.ArgumentParser(
        description="Convert NECID CSV results to NeoaPred input format."
    )
    parser.add_argument("input_csv", type=Path, help="NECID results CSV file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="output CSV path (default: <input>_neoapred.csv)",
    )
    args = parser.parse_args()

    output_path = args.output or args.input_csv.with_name(
        "{}_neoapred.csv".format(args.input_csv.stem)
    )
    try:
        convert(args.input_csv, output_path)
    except ValueError as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()