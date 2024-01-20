from gromacs_protocol import GromacsPprotocol


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--pdb_directory',
        type=str,
        help='directory with the fasta files to process',
    )
    parser.add_argument(
        '--output_directory',
        type=str,
        help='output directory',
    )
    parser.add_argument(
        '--mdp_directory',
        type=str,
        help='files: ions.mdp, minim.mdp,'
             f' nvt.mdp,npt.mdp should be in the mdp directory',
    )
    parser.add_argument(
        '--exe',
        type=str,
        help='executable',
        default='gmx',
    )
    args = parser.parse_args()
    return args


def main():
    args = argument_parser()
    gromacs_prot = GromacsProtocol(
        pdb_directory=args.pdb_directory,
        output_directory=args.output_directory,
        mdp_directory=args.mdp_directory,
        exe=args.exe
    )
    gromacs_prot.main(
        args_list=gromacs_prot.identifiers_list
    )


if __name__ == '__main__':
    main()
