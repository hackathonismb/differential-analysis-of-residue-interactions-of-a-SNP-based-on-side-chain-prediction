"""Hetatm near mutation.

Checks to see if a hetatm is near any atom from a given residue. Hetatms
represent atoms in non-standard residues"

Input:
1. File Path
2. Mutation Residue Chain
3. Mutation Residue Number

Output:
1. 1 if any mutation atoms are within the distance, else 0
"""
import argparse
import math


def argument_parser():
    """Parse arguments for the hetatm_near_mutation.py scipt."""
    parser = argparse.ArgumentParser()
    required_arguments = parser.add_argument_group("Required Arguments")
    required_arguments.add_argument("-i",
                                    "--input_file",
                                    required=True,
                                    type=str,
                                    help="Input PDB file")
    required_arguments.add_argument("-r",
                                    "--residue_number",
                                    type=int,
                                    required=True,
                                    help="Residue number of the mutation.")
    required_arguments.add_argument("-c",
                                    "--residue_chain",
                                    type=str,
                                    required=True,
                                    help="Residue chain of the mutation.")
    parser.add_argument("-t",
                        "--threshold_distance",
                        default=8,
                        type=float,
                        help=("Distance in angstroms to search for hetatms"
                              "from a given residue"))
    args = parser.parse_args()
    return args


def distance(coord_a, coord_b):
    """Return distance between 2 sets of coordinates."""
    # Check that cordinants are same dimension
    if len(coord_a) != len(coord_b):
        raise ValueError("Coordinates are not the same length.")
    # Sum the square of differences
    sum_sqr_diff = 0
    for i in range(len(coord_a)):
        sum_sqr_diff += (coord_b[i] - coord_a[i])**2
    return math.sqrt(sum_sqr_diff)


def mutant_hetatm_coordinate_parse(file_path,
                                   query_chain,
                                   query_residue_number):
    """Parse a PDB file for mutant and hetatm coordinates.

    Given a file_path for a PDB file, the chain of the mutant residue, and
    its residue number, returns a tuple of mutant and hetatm coordinates.

    PDB files information is specified by a range of columns in the file a
    category of information can be. The information for chains, residue values,
    and coordinates are spliced from the line read for the column range they are
    given in, rather than spliting by a given delimiter.
    """
    mutant_coords, hetatm_coords = [], []
    # Open a file and loop over each line
    file = open(file_path, 'r')
    for line in file:
        # For ATOM rows, validate if it matches the input residue and chain
        if line.startswith("ATOM"):
            chain = line[21]
            residue_number = int(line[22:26].strip())
            if query_chain == chain and query_residue_number == residue_number:
                coords = list(map(float, line[30:54].split()))
                mutant_coords.append(coords)
        elif line.startswith("HETATM"):
            coords = list(map(float, line[30:54].split()))
            hetatm_coords.append(coords)
    file.close()
    return (mutant_coords, hetatm_coords)


def any_within_distance_threshold(m_coords, h_coords, threshold):
    """Return boolean result for if any coordiantes are within a threshold."""
    for m_coord in m_coords:
        for h_coord in h_coords:
            if distance(m_coord, h_coord) <= threshold:
                return True
    return False


def main():
    """Print 1 if any hetatms are near a mutation atom otherwise 0."""
    # Parse arguments
    args = argument_parser()
    # Get coordinates for mutant residue and hetatms
    m_coords, h_coords = mutant_hetatm_coordinate_parse(args.input_file,
                                                        args.residue_chain,
                                                        args.residue_number)
    # Print 0 if any mutants are under the threshod, else print 1
    if any_within_distance_threshold(m_coords,
                                     h_coords,
                                     args.threshold_distance):
        print(1)
    else:
        print(0)


if __name__ == "__main__":
    main()

