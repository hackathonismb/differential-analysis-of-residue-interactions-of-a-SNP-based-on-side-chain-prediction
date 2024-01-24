import os


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fasta_directory',
        type=str,
        help='directory with the fasta files to process',
    )
    parser.add_argument(
        '--output_directory',
        type=str,
        help='output directory',
    )
    parser.add_argument(
        '--alphafold_dataset',
        type=str,
        help='alphafold dataset',
    )
    parser.add_argument(
        '--alphafold_directory',
        type=str,
        help='alphafold directory',
    )
    args = parser.parse_args()
    return args


def main():
    args = argument_parser()
    fasta_directory = args.fasta_directory
    output_directory = args.output_directory
    alphafold_dataset = args.alphafold_dataset
    alphafold_directory = args.alphafold_directory
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    filenames_list = [
        filename
        for filename in os.listdir(fasta_directory)
        if '.fasta' in filename
    ]
    identifiers_list = [
        filename.rsplit('.fasta')[0]
        for filename in filenames_list
    ]
    args_list = [
        identifier
        for identifier in identifiers_list
        if not os.path.exists(f'{output_directory}/{identifier}')
    ]
    if args_list:
        for identifier_filename in filenames_list:
            identifier = identifier_filename.rsplit('.fasta')[0]
            input_path = os.path.join(output_directory, identifier)
            out_path = os.path.join(fasta_directory, identifier_filename)
            # You need to adapt the following command to your
            # alphafold installation, the following is tailored for running
            # at uppmax
            alphafold_command = [
                'bash',
                os.path.join(alphafold_directory, 'run_alphafold.sh'),
                f'-d {alphafold_dataset}',
                f'--fasta_paths={input_path}',
                f'--output_dir={out_path}',
                '--model_preset=monomer',
                '--db_preset=full_dbs'
            ]
            os.system(' '.join(alphafold_command), shell=True)


if __name__ == "__main__":
    main()
