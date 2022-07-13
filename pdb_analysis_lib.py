"""PDB_Analysis_Library."""

import math


PDB_INDEX_DELIMS = [0, 6, 11, 16, 17, 20, 22, 26, 27, 38, 46, 54, 60, 66, 76, 78]
PDB_COLUMN_NAMES = ["Molecule Type",
                    "Atom serial number",
                    "Atom name",
                    "Alternate location indicator",
                    "Residue name",
                    "Chain identifier",
                    "Residue sequence number",
                    "Code for insertions of residues",
                    "X",
                    "Y",
                    "Z",
                    "Occupancy",
                    "Temperature factor",
                    "Segment idenifier",
                    "Element symbol"]


def distance(coord_a, coord_b):
    """Return distance between 2 sets of coordinates."""
    # Check that cordinants are same dimension
    if len(coord_a) != len(coord_b):
        raise ValueError("Coordinates are not the same length.")
    # Sum the square of differences
    sum_sqr_diff = 0
    for i in range(len(coord_a)):
        sum_sqr_diff += (float(coord_b[i]) - float(coord_a[i]))**2
    return math.sqrt(sum_sqr_diff)


def pdb_row_to_list(atm_row: str) -> str:
    """Returns a list of each item in an atm or hetatm row."""
    delim_indicies = PDB_INDEX_DELIMS.copy()
    result = [atm_row[i: j].strip() for i, j in zip(delim_indicies, delim_indicies[1:])]
    return result


def read_pdb_atms_hetatms(pdb_lines):
    result = []
    for line in pdb_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            result.append(pdb_row_to_list(line))
    return result


def min_distance(query_coord, coord_list):
    result = math.inf
    for coord in coord_list:
        dist = distance(query_coord, coord)
        if dist <= result:
            result = dist
    return result
