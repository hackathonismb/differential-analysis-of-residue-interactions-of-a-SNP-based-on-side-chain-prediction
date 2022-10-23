"""Distance to features.

This script takes a PDB file and a features file with parameters indicating
the chain and residue of the AA and returns a table of distances to
features outlined in uniprot.
"""

import argparse
import glob
import os
import sys

import pdb_analysis_lib as pal


def argument_parser():
    """Parse arguments for the distance_to_features script."""
    parser = argparse.ArgumentParser()
    parser.add_argument("ember_folder",
                                    type=str,
                                    help="Folder of original EMBER3D files.")
    parser.add_argument("wt_folder",
                                    type=str,
                                    help="Folder of corresponding WT pdb files.")
    parser.add_argument("output_folder",
                                    type=str,
                                    help="Output folder.")
    args = parser.parse_args()
    return args


def main():
    args = argument_parser()
    ember_folder = args.ember_folder.rstrip('/')
    wt_folder = args.wt_folder.rstrip('/')

    output_folder = args.output_folder.rstrip('/')
    os.makedirs(output_folder, exist_ok=True)
    ember_files = glob.glob(f"{ember_folder}/*.pdb")
    for ember_file_path in ember_files:
        file_name = ember_file_path.split('/')[-1]
        pdb_id, chain = file_name.split('_')[0:2]
        try:
            pal.correct_ember_file(ember_file_path, f"{wt_folder}/{pdb_id}_WT.pdb", chain, f"{output_folder}/{file_name}")
        except:
            print(f"Issue correcting residue numbers for: {ember_file_path}", file=sys.stderr)
            raise


if __name__ == "__main__":
    main()
