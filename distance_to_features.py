"""Distance to features.

This script parses the pdb file result of a mutation analysis from
iCn3D for a single amino acid mutation and returns a csv output to
stdout of the distance of all hetatms.
"""

import argparse
import json

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
    required_arguments.add_argument("-f",
                                    "--features",
                                    required=True,
                                    type=str,
                                    help="File of features in JSON format.")
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
    pdb_data = pal.read_pdb_atms_hetatms(pdb_lines)
    mut_data = pal.parse_data_by_residues(pdb_data, [mut_residue])
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
    result = list(map(lambda x: ["{0:.1f}".format(x[0])] + x[1:], result))
    # Output result
    print("Minimum Distance, Uniprot ID, Feature, Residue")
    for i in result:
        print(','.join(i))


if __name__ == "__main__":
    main()
