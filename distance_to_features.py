"""Distance to features.

This script takes a PDB file and a features file with parameters indicating
the chain and residue of the AA and returns a table of distances to
features outlined in uniprot.
"""

import argparse

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
    result = pal.distance_to_features(args.pdb_file, args.features, args.chain, args.residue)
    # Format the distances to a single decimal place
    result = list(map(lambda x: ["{0:.1f}".format(x[0])] + x[1:], result))
    # Output result
    print("Minimum Distance, Uniprot ID, Feature, Residue")
    for i in result:
        print(','.join(i))


if __name__ == "__main__":
    main()
