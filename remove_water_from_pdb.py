"""Remove water atoms from a pdb file."""

import sys
import os


def main():
    """Remove HOH atoms from pdn file and return results."""
    pdb_file = sys.argv[1]
    if pdb_file.endswith('.pdb'):
        output = pdb_file[:-4]
    else:
        output = pdb_file
    os.system(f"grep -v HOH {pdb_file} > {output}_clean.pdb")


if __name__ == "__main__":
    main()
