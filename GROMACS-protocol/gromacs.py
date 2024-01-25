from gromacs_protocol import GromacsProtocol
import argparse


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
        '--hard_force',
        action='store_true',
        default=False,
        help='Remove all files in output directory if it exists.'
    )
    parser.add_argument(
        '--GMX',
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
        GMX=args.GMX,
        hard_force=args.hard_force,
    )
    gromacs_prot.main(
        args_list=gromacs_prot.identifiers_list
    )


if __name__ == '__main__':
    main()

# python3 gromacs.py --pdb_directory ../test_directory_pdb --output_directory ../test_directory_output --mdp_directory /proj/naiss2023-22-736/msarrias/projects/ismb/differential-analysis-of-residue-interactions-of-a-SNP-based-on-side-chain-prediction/GROMACS-protocol