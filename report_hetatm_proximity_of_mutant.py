"""Report hetatm proximity to mutant.

This script parses the pdb file result of a mutation analysis from
iCn3D for a single amino acid mutation and returns a csv output to
stdout of the distance of all hetatms.
"""

import argparse

import pdb_analysis_lib as pal


def argument_parser():
    """Parse arguments for the hetatm_near_mutation.py scipt."""
    parser = argparse.ArgumentParser()
    required_arguments = parser.add_argument_group("Required Arguments")
    required_arguments.add_argument("-i",
                                    "--input_file",
                                    required=True,
                                    type=str,
                                    help="Input PDB file")
    parser.add_argument("-t",
                        "--threshold_distance",
                        default=8,
                        type=float,
                        help=("Distance from the mutant residue in angstroms"
                              "to search for hetatms"))
    args = parser.parse_args()
    return args


def main():
    """Print to stdout distance to mutation and atom information of hetatms."""
    args = argument_parser()
    file_object = open(args.input_file)
    file_data = file_object.readlines()
    file_object.close()
    title_line = file_data[1]
    mut_chain, mut_residue = title_line.split()[3][:-1].split('_')
    atm_results = pal.read_pdb_atms_hetatms(file_data)
    mutant_residues = list(filter(lambda x: x[5] == mut_chain and x[6] == mut_residue,
                             atm_results))
    mutant_coordiantes = list(map(lambda x: x[8:11], mutant_residues))
    hetatm_residues = list(filter(lambda x: x[0] == "HETATM", atm_results))
    # Add distance information to hetatm
    result = []
    for hetatm in hetatm_residues:
        min_dist = pal.min_distance(hetatm[8:11], mutant_coordiantes)
        result.append(["{0:.1f}".format(min_dist)] + hetatm)
    print(",".join(["Minimum Distance From Mutation"] + pal.PDB_COLUMN_NAMES))
    for i in result:
        print(",".join(i))


if __name__ == "__main__":
    main()
