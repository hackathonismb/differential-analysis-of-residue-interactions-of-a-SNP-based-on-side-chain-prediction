"""Given a folder with WT and mutations, return fasta files."""

import argparse
import glob
import logging

import pdb_analysis_lib as pal


def argument_parser():
    parser = argparse.ArgumentParser()
    required_arguments = parser.add_argument_group("Required Arguments")
    parser.add_argument('-v',
                        "--verbose",
                        action="store_true")
    required_arguments.add_argument('-i',
                                    "--input_folder",
                                    required=True,
                                    type=str,
                                    help="Input folder of PDB files.")
    required_arguments.add_argument('-o',
                                    "--output_folder",
                                    required=True,
                                    type=str,
                                    help="Output folder for FASTA files.")
    args = parser.parse_args()
    return args


def write_fasta_file(file_path, header, sequence, max_sequence_line_length=79):
    with open(file_path, 'w') as file_object:
        file_object.write(f">{header}\n")
        start_idx = 0
        end_idx = max_sequence_line_length
        while end_idx < len(sequence):
            file_object.write(sequence[start_idx:end_idx] + '\n')
            start_idx += max_sequence_line_length
            end_idx += max_sequence_line_length
        file_object.write(sequence[start_idx:] + '\n')


def main():
    args = argument_parser()
    input_folder = args.input_folder.rstrip('/')
    output_folder = args.output_folder.rstrip('/')
    files = glob.glob(f"{input_folder}/*.pdb")
    mutation_files = list(filter(lambda x: "WT" not in x.split('/')[-1],
                                 files))
    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.WARNING
    logging.basicConfig(level=log_level)
    id_chain_pairs = set()
    logging.info(f"Input of {len(files)} files")
    logging.info(f"Converting {len(mutation_files)} mutation files")
    for mfile in mutation_files:
        mfile_name = mfile.split('/')[1]
        pdb_id, chain, residue, mutation = mfile_name.rstrip('.pdb').split('_')
        id_chain_pairs.add((pdb_id, chain))
        try:
            header, sequence = pal.pdb_to_fasta(mfile, chain)
        except KeyError as keyerror:
            logging.warning(f"{mfile_name} not converted, unconventional AA: {keyerror}")
        except ValueError as valerr:
            logging.warning(f"{mfile_name} not converted: {valerr}")
        except:
            logging.critical(f"{mfile_name} cannot be processed.")
            raise
        else:
            # Output file
            output_file_name = f"{pdb_id}_{chain}_{residue}_{mutation}.fasta"
            output_file = f"{output_folder}/{output_file_name}"
            write_fasta_file(output_file, header, sequence)
    logging.info("Done!")
    # Convert WT files
    wt_files = list(filter(lambda x: "WT" in x.split('/')[-1], files))
    wt_ids = list(map(lambda x: x.split('/')[-1].split('_')[0], wt_files))
    wt_file_dict = dict(zip(wt_ids, wt_files))
    id_chain_pairs = list(id_chain_pairs)
    logging.info(f"Converting {len(wt_files)} wildtype files...")
    for pdb_id, chain in id_chain_pairs:
        wtfile_name = wt_file_dict[pdb_id]
        try:
            header, sequence = pal.pdb_to_fasta(wtfile_name, chain)
        except KeyError as keyerror:
            logging.warning(f"{wtfile_name} not converted, unconventional AA: {keyerror}")
        except ValueError as valerr:
            logging.warning(f"{wtfile_name} not converted: {valerr}")
        except:
            logging.critical(f"{wtfile_name} cannot be processed.")
            raise
        else:
            output_file_name = f"{pdb_id}_{chain}_WT.fasta"
            output_file = f"{output_folder}/{output_file_name}"
            write_fasta_file(output_file, header, sequence)
    logging.info("Done!")


if __name__ == "__main__":
    main()
