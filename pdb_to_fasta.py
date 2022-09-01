"""PDB to FASTA.

This script takes a PDB file as input and returns a FASTA file with
the file name without extension as the header and the amino acid sequence
in single-letter format as the sequence.
"""

import argparse
import re

import pdb_analysis_lib as pal


def argument_parser():
    """Parse arguments for the distance_to_features scipt."""
    parser = argparse.ArgumentParser()
    required_arguments = parser.add_argument_group("Required Arguments")
    required_arguments.add_argument("-p",
                                    "--pdb_file",
                                    required=True,
                                    type=str,
                                    help="Input PDB file")
    required_arguments.add_argument("-c",
                                    "--chain",
                                    required=True,
                                    type=str,
                                    help="Structural chain to convert.")
    args = parser.parse_args()
    return args


def main():
    """Convert PDB sequence to FASTA."""
    args = argument_parser()
    # Read file
    file_object = open(args.pdb_file)
    file_data = file_object.readlines()
    file_object.close()
    # Filter for ATOM lines
    file_data = list(filter(lambda x: x.startswith("ATOM"), file_data))
    # Split rows by data type
    file_data = list(map(pal.pdb_row_to_list, file_data))
    # Filter by chain
    file_data = list(filter(lambda x: x[5] == args.chain, file_data))
    file_data = sorted(file_data, key=lambda x: int(x[6]))
    # Convert AAs listed in PDB file to single letter format
    # Initialize with first letter and residue number
    fasta_sequence = pal.AA_DICT[file_data[0][4]]
    residue = int(file_data[0][6])
    for row in file_data:
        # If the next row has a new residue number, add it to the sequence
        if int(row[6]) == residue + 1:
            fasta_sequence = fasta_sequence + pal.AA_DICT[row[4]]
            residue += 1
    # Output FASTA header and sequence
    header = re.sub(r'^.*/', '', args.pdb_file)
    header = header.rstrip(".pdb")
    print(f">{header}:{args.chain}")
    print(fasta_sequence)


if __name__ == "__main__":
    main()
