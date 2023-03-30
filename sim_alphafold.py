import os

PROJ_DIR = '/proj/snic2022-22-771/users/x_maher'
DATA_DIR = f'{PROJ_DIR}/data/test_fasta'
WORKING_DIR = f'{PROJ_DIR}/ALPHAFOLD'
DUMP_DIR = f'{WORKING_DIR}/Benchmark_ALPHAFOLD'
filenames_list = [filename for filename in os.listdir(DATA_DIR) if '.fasta' in filename]
identifiers_list =  [filename.rsplit('.fasta')[0] for filename in filenames_list]
args_list = [identifier for identifier in identifiers_list if not os.path.exists(f'{DUMP_DIR}/{identifier}')]  

if args_list:
    for identifier_filename in filenames_list:
        identifier = identifier_filename.rsplit('.fasta')[0]
        sim_direct = f'{DUMP_DIR}/{identifier}'
        file_path = f'{DATA_DIR}/{identifier_filename}'
        alphafold_call = ['run_singularity_af',
                          f'--fasta_paths={file_path}',
                          f'--output_dir={sim_direct}',
                          '--model_preset=monomer', 
                          '--db_preset=full_dbs',
                          '--bfd_database_path=/data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt',
                          '--pdb70_database_path=/data/pdb70/pdb70',
                          '--uniclust30_database_path=/data/uniclust30/uniclust30_2018_08/uniclust30_2018_08',
                          '--max_template_date=2020-05-14',
                          '--use_gpu_relax=True']
        os.system(' '.join(alphafold_call)) 
    