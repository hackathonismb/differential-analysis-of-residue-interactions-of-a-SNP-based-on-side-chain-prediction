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
    parser.add_argument(
        '--ntmpi',
        type=str,
        help='The number of MPI threads to use.'
             ' This sets how many separate processes will be '
             'spawned for parallel computation. Typically, '
             'this should not exceed the number of physical '
             'or logical CPU cores available.',
        default=""
    )
    parser.add_argument(
        '--ntomp',
        type=str,
        help='The number of OpenMP threads per MPI process.'
             ' This controls the parallelism within each MPI process.'
             ' It is often set to the number of CPU cores per MPI '
             'task or to suit the specific computational workload.',
        default='0'
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
        ntomp=args.ntomp,
        ntmpi=args.ntmpi
    )
    gromacs_prot.main(
        args_list=gromacs_prot.identifiers_list
    )


if __name__ == '__main__':
    main()
