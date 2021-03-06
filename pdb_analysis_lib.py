"""PDB_Analysis_Library.

A local python module for functions used for scripts analyzing PDB files.


Global variables:
PDB_INDEX_DELIMS: Indexes to split a string when parsing an ATOM or HETATM
PDB_COLUMN_NAMES: Names corresponding to each substring from PDB_INDEX_DELIMS

Functions:
distance: Return the euclidean distance between 2 sets of coordinates.
min_distance: Return the minimum distance of one coordinate to others.
pdb_row_to_list: Return a list of each item in an atm or hetatm row.
read_pdb_atms_hetatms: Return a nested list of PDB ATOM and HETATM rows.
"""

import math


PDB_INDEX_DELIMS = [0,
                    6,
                    11,
                    16,
                    17,
                    20,
                    22,
                    26,
                    27,
                    38,
                    46,
                    54,
                    60,
                    66,
                    76,
                    78]
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


def distance(coord_a: list, coord_b: list) -> float:
    """Return the euclidean distance between 2 sets of coordinates.

    Given two equal lengthnumeric iterables, calculate the euclidean distance
    between them and return that value.

    :param coord_a: List of numeric values
    :param coord_b: List of numeric values with the same length as coord_a
    :return: Euclidean distance between the sets of coordinates.
    """
    # Check that cordinants are same dimension
    if len(coord_a) != len(coord_b):
        raise ValueError("Coordinates are not the same length.")
    # Sum the square of differences
    sum_sqr_diff = 0
    for i in range(len(coord_a)):
        sum_sqr_diff += (float(coord_b[i]) - float(coord_a[i]))**2
    return math.sqrt(sum_sqr_diff)


def min_distance(query_coord: list, coord_list: list) -> float:
    """Return the minimum distance of one coordinate to a list of others.

    A one-to-many comparison of a coordinate to a list of coordinates where the
    minimum distance value is returned.

    :param query_coord: Iterable of numeric values
    :type query_coord: list
    :param coord_list: Iterable of numeric value iterables
    :type coord_list: list
    :return: Minimum distance from the one-to-many comparison
    """
    result = math.inf
    for coord in coord_list:
        dist = distance(query_coord, coord)
        if dist <= result:
            result = dist
    return result


def pdb_row_to_list(row: str) -> str:
    """Return a list of each item in an atm or hetatm row."""
    delim_idxs = PDB_INDEX_DELIMS.copy()
    result = [row[i: j].strip() for i, j in zip(delim_idxs, delim_idxs[1:])]
    return result


def read_pdb_atms_hetatms(pdb_lines: list) -> list:
    """Return a nested list of PDB ATOM and HETATM rows.

    Given a list of strings representing rows in a PDB file, return
    all ATOM and HETATM rows as parsed lists.

    :param pdb_lines: List of strings representing lines in a PDB file
    :type pdb_lines: list
    :return: Nested list of parsed ATOM and HETATM rows
    """
    result = []
    for line in pdb_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            result.append(pdb_row_to_list(line))
    return result
