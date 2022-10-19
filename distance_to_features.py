"""Distance to features.

This script takes a PDB file and a features file with parameters indicating
the chain and residue of the AA and returns a table of distances to
features outlined in uniprot.
"""

import argparse
import json
import math

import pdb_analysis_lib as pal


def argument_parser():
    """Parse arguments for the distance_to_features script."""
    parser = argparse.ArgumentParser()
    required_arguments = parser.add_argument_group("Required Arguments")
    required_arguments.add_argument("pdb_file",
                                    type=str,
                                    help="Input PDB file")
    required_arguments.add_argument("-f",
                                    "--features",
                                    required=True,
                                    type=str,
                                    help="File of features in JSON format.")
    required_arguments.add_argument("-c",
                                    "--chain",
                                    required=True,
                                    type=str,
                                    help="Chain the residue is part of.")
    required_arguments.add_argument("-r",
                                    "--residue",
                                    required=True,
                                    type=str,
                                    help="Residue number to check distances from.")
    args = parser.parse_args()
    return args


def main():
    """Main function."""
    args = argument_parser()
    # Read pdb file
    pdb_file_object = open(args.pdb_file)
    pdb_lines = pdb_file_object.readlines()
    pdb_file_object.close()
    # Read features file
    features_file_object = open(args.features)
    features_dict = json.load(features_file_object)
    features_file_object.close()
    # Parse mutation information
    mut_chain, mut_residue = pal.parse_mutation_title(pdb_lines[1])
    # Filter for PDB atom data
    pdb_data = pal.read_pdb_atms(pdb_lines, ("ATOM"))
    mut_data = pal.parse_data_by_residues(pdb_data, [args.residue])
    mut_data = pal.parse_data_by_chains(mut_data, [args.chain])
    mut_coords = pal.coords_from_pdb_data(mut_data)
    # For each uniprot_id
    seen_residue = {}  # Track previously calculated distances
    result = []
    for uniprot_id in features_dict:
        for feature in features_dict[uniprot_id]:
            if feature not in ("Region"):
                for residue in features_dict[uniprot_id][feature]:
                    # If the distance was already calculated, append to result
                    if residue in seen_residue:
                        result.append([seen_residue[residue],
                                       uniprot_id,
                                       feature,
                                       str(residue)])
                    # Else, calculate the distance, append to result, and save
                    else:
                        residue_entries = pal.parse_data_by_residues(pdb_data,
                                                                     [str(residue)])
                        feature_coords = pal.coords_from_pdb_data(residue_entries)
                        dist = pal.min_distance_coords_to_coords(mut_coords,
                                                                 feature_coords)
                        result.append([dist, uniprot_id, feature, str(residue)])
                        seen_residue[residue] = dist
    # Sort result by distance
    result = sorted(result, key=lambda x: x[0])
    # Remove all infinity distances
    result = list(filter(lambda x: x[0] != math.inf, result))
    # Format the distances to a single decimal place
    result = list(map(lambda x: ["{0:.1f}".format(x[0])] + x[1:], result))
    # Output result
    print("Minimum Distance, Uniprot ID, Feature, Residue")
    for i in result:
        print(','.join(i))


if __name__ == "__main__":
    main()
