"""Distance to features.

This script takes a PDB file and a features file with parameters indicating
the chain and residue of the AA and returns a table of distances to
features outlined in uniprot.
"""

import argparse
import glob

import pdb_analysis_lib as pal


def argument_parser():
    """Parse arguments for the distance_to_features script."""
    parser = argparse.ArgumentParser()
    required_arguments = parser.add_argument_group("Required Arguments")
    required_arguments.add_argument("ember_folder",
                                    type=str,
                                    help="Input PDB file")
    required_arguments.add_argument("wt_folder",
                                    type=str,
                                    help="File of features in JSON format.")
    required_arguments.add_argument("output_folder",
                                    type=str,
                                    help="Chain the residue is part of.")
    args = parser.parse_args()
    return args


def main():
    args = argument_parser()
    ember_folder = args.ember_folder.rstrip('/')
    wt_folder = args.wt_folder.rstrip('/')
    output_folder = args.output_folder.rstrip('/')
    os.makedirs(output_folder, exist_ok=True)
    ember_files = glob.glob(f"{ember_folder}/*.pdb")
    for i in ember_files:
        file_name = i.split('/')[-1]
        pdb_id, chain = file.split('_')[0:2]
        pal.correct_ember_file(i, f"{wt_folder}/{pdb}_{chain}_WT.pdb", chain, f"{output_folder}/{file}")
