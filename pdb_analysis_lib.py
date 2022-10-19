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
import re


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

AA_DICT = {"ALA":"A",
           "CYS":"C",
           "ASP":"D",
           "GLU":"E",
           "PHE":"F",
           "GLY":"G",
           "HIS":"H",
           "ILE":"I",
           "LYS":"K",
           "LEU":"L",
           "MET":"M",
           "ASN":"N",
           "PRO":"P",
           "GLN":"Q",
           "ARG":"R",
           "SER":"S",
           "THR":"T",
           "VAL":"V",
           "TRP":"W",
           "TYR":"Y"}


def distance(coord_a: list, coord_b: list) -> float:
    """Return the euclidean distance between 2 sets of coordinates.

    Given two equal length numeric iterables, calculate the euclidean distance
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


def min_distance_coord_to_coords(query_coord: list, coord_list: list) -> float:
    """Return the minimum distance of one coordinate to a group of others.

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


def min_distance_coords_to_coords(query_coords: list, coord_list: list) -> float:
    """Return the minimum distance of one group of coordinates to another."""
    result = math.inf
    for query_coord in query_coords:
        dist = min_distance_coord_to_coords(query_coord, coord_list)
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

def read_pdb_atms(pdb_lines: list, mol_types: tuple = ("ATOM", "HETATM")) -> list:
    """Return a nested list of PDB ATOM and HETATM rows.

    Given a list of strings representing rows in a PDB file, return
    all rows of specified molecule types.

    :param pdb_lines: List of strings representing lines in a PDB file
    :type pdb_lines: list
    :return: Nested list of parsed molecule type rows
    """
    result = []
    for line in pdb_lines:
        if any([line.startswith(mtype) for mtype in mol_types]):
            result.append(pdb_row_to_list(line))
    return result

# Convert these into a single function 
def parse_data_by_residues(pdb_data, residues):
    """Return subset of pdb data containing specified residues."""
    result = list(filter(lambda x: x[6] in residues, pdb_data))
    return result

def parse_data_by_chains(pdb_data, chains):
    """Return subset of pdb data containing specified chains"""
    result = list(filter(lambda x: x[5] in chains, pdb_data))
    return result
###
def coords_from_pdb_data(pdb_data):
    """Return just coordinates from pdb data."""
    result = list(map(lambda x: x[8:11], pdb_data))
    return result

def parse_mutation_title(title_line):
    """Parse the mutation title for the chain and residue."""
    chain, residue = title_line.split()[3][:-1].split('_')
    return (chain, residue)

def pdb_to_fasta(pdb_file, chain):
    """Convert PDB sequence to FASTA."""
    # Read file
    file_object = open(pdb_file)
    file_data = file_object.readlines()
    file_object.close()
    # Filter for ATOM lines
    file_data = list(filter(lambda x: x.startswith("ATOM"), file_data))
    # Split rows by data type
    file_data = list(map(pdb_row_to_list, file_data))
    # Filter by chain
    file_data = list(filter(lambda x: x[5] == chain, file_data))
    # Check that the chain exists
    if len(file_data) == 0:
        raise ValueError("No AA to convert.")
    #file_data = sorted(file_data, key=lambda x: int(x[6]))
    # Convert AAs listed in PDB file to single letter format
    # Initialize with first letter and residue number
    fasta_sequence = AA_DICT[file_data[0][4]]
    completed_residues = set()
    residue = file_data[0][6]
    for row in file_data:
        if row[6] in completed_residues: # If there is a repeat residue number
            raise ValueError("Repeated residues found.")
        if row[6] != residue: # If there is a new residue
            fasta_sequence = fasta_sequence + AA_DICT[row[4]]
            completed_residues.add(residue)
            residue = row[6]
    # Output FASTA header and sequence
    header = re.sub(r'^.*/', '', pdb_file)
    header = header.rstrip(".pdb")
    return (header, fasta_sequence)
