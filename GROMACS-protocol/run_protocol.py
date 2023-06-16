from protocol import GROMACS_protocol 

PROJ_DIR = '/proj/snic2022-22-771/users/x_maher'
DATA_DIR = f'{PROJ_DIR}/data/WT_Benchmark_PDB'
WORKING_DIR = f'{PROJ_DIR}/GROMACS'
DUMP_DIR = f'{WORKING_DIR}/WT_Benchmark_PDB_GROMACS'

if __name__ == '__main__':
    batch = GROMACS_protocol(DATA_DIR,PROJ_DIR,DUMP_DIR,WORKING_DIR, exe='gmx')
    batch.thread_simulations(batch.identifiers_list)