import os, sys, time
import tempfile 
from subprocess import PIPE, DEVNULL, STDOUT, Popen, run, TimeoutExpired, call
from multiprocessing import Process


class GROMACS_protocol:
    def __init__(self,
                 DATA_DIR,
                 PROJ_DIR,
                 DUMP_DIR,
                 WORKING_DIR):
        self.DATA_DIR = DATA_DIR
        self.PROJ_DIR = PROJ_DIR
        self.DUMP_DIR = DUMP_DIR
        self.WORKING_DIR = WORKING_DIR
        self.filenames_list = [filename for filename in os.listdir(self.DATA_DIR) if 'pdb' in filename]
        self.identifiers_list =   [filename.rsplit('.pdb')[0] for filename in self.filenames_list]
        self.files_req = [f'{self.WORKING_DIR}/ions.mdp',
                          f'{self.WORKING_DIR}/minim.mdp',
                          f'{self.WORKING_DIR}/nvt.mdp',
                          f'{self.WORKING_DIR}/npt.mdp']
        if not all([os.path.exists(file_path) for file_path in self.files_req]):
            raise ValueError(f'files: ions.mdp, minim.mdp, nvt.mdp,npt.mdp should be in working directory: {self.WORKING_DIR}')
        if len(self.identifiers_list) != len(set(self.identifiers_list)):
            raise ValueError('all identifiers should be distinct')
            
            
    @staticmethod
    def suprocess_call(command):
        ec = call(command, stdout=DEVNULL, stderr=STDOUT)
        if ec == 1:
            comm_ = ' '.join(command)
            sys.exit(f'error while calling: {comm_}')
            
            
    @staticmethod
    def interactive_process(command, inpt):
        proc = Popen(command, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        try:
            outs, errs = proc.communicate(input= inpt, timeout=120)
            proc.kill()
        except TimeoutExpired:
            proc.kill()
            comm_ = ' '.join(command)
            sys.exit(f'error while calling: {comm_}')
            
            
    @staticmethod
    def batch(iterable, n=1):
        l = len(iterable)
        for ndx in range(0, l, n):
            yield iterable[ndx:min(ndx + n, l)]
            
            
    def run_protocol(self, identifier):
        sim_direct = f'{self.DUMP_DIR}/{identifier}'
        os.mkdir(sim_direct)
        os.chdir(sim_direct)
        time.sleep(5)
        node = '/software/sse/easybuild/prefix/software/nodejs/16.15.1-GCCcore-11.3.0-nsc1/bin/node'
        # STEP 1: The following step will remove all hetero atoms, use one chain name, and add missing atoms
        if not os.path.exists(f'{identifier}_nohet.pdb'):
            rmhet_command = [node, f'{self.WORKING_DIR}/rmhet.js', f'{self.DATA_DIR}/{identifier}.pdb']
            with open(f'{identifier}_nohet.pdb', 'w') as out:
                run(rmhet_command, stdout=out)
        if not os.path.exists(f'{identifier}_clean.pdb'):
            addmissingatom_command = [node,
                                      f'{self.WORKING_DIR}/addmissingatoms.js',
                                      f'{identifier}_nohet.pdb']
            with open(f'{identifier}_clean.pdb', 'w') as out:
                run(addmissingatom_command, stdout=out)
        # STEP 2: Here we select AMBER99SB-ILDN force-field (#6 in the list) 
        # for describing atom-atom interaction energies to be used in the GROMACS simulations.
        if not os.path.exists(f'{identifier}_processed.gro'):
            pdb2gmx_command = ['gmx_mpi_d', 'pdb2gmx',
                               '-f', f'{identifier}_clean.pdb',  
                               '-o', f'{identifier}_processed.gro',
                               '-water', 'spce']
            self.interactive_process(pdb2gmx_command, b'6')
        # STEP 3: Here we create a simulation box 0.5nm from the protein edge in all 6 directions.     
        if not os.path.exists(f'{identifier}_newbox.gro'):
            editconf_command = ['gmx_mpi_d', 'editconf',
                                '-f', f'{identifier}_processed.gro',  
                               '-o', f'{identifier}_newbox.gro',
                               '-c', '-d', '0.5', '-bt', 'cubic']
            self.suprocess_call(editconf_command)
        # STEP 4: Here we add water molecules in the simulation box and new topology file is written.
        if not os.path.exists(f'{identifier}_solv.gro'):
            solvate_command = ['gmx_mpi_d', 'solvate',
                               '-cp', f'{identifier}_newbox.gro',  
                               '-cs', 'spc216.gro',
                               '-o', f'{identifier}_solv.gro',
                               '-p', f'topol.top']
            self.suprocess_call(solvate_command)
        # STEP 5: Generate parameter file ions.tpr for all atoms.
        if not os.path.exists(f'ions.tpr'):
            grompp_command = ['gmx_mpi_d', 'grompp',
                              '-f', f'{self.WORKING_DIR}/ions.mdp',  
                              '-c', f'{identifier}_solv.gro',
                              '-p', 'topol.top', 
                              '-o', 'ions.tpr',
                              '-maxwarn', '2']
            self.suprocess_call(grompp_command)
        # STEP 6: This will add ions to your simulation box and replace 
        # the solvent molecules if there is a clash.
        if not os.path.exists(f'{identifier}_solv_ions.gro'):
            genion_command = ['gmx_mpi_d', 'genion',
                              '-s', f'ions.tpr',  
                              '-o', f'{identifier}_solv_ions.gro',
                              '-p', f'topol.top',
                              '-pname', 'NA',
                              '-nname', 'CL', '-neutral']
            self.interactive_process(genion_command, b'13')
        # STEP 7: Energy minimization
        if not os.path.exists('em.tpr'):
            em_command = ['gmx_mpi_d', 'grompp',
                          '-f', f'{self.WORKING_DIR}/minim.mdp',  
                          '-c', f'{identifier}_solv_ions.gro',
                          '-p', 'topol.top', 
                          '-o', 'em.tpr']
            self.suprocess_call(em_command)
        # STEP 8: This runs the energy minimization.
        run_em_command = ['gmx_mpi_d', 'mdrun', '-v', '-deffnm', 'em']
        self.suprocess_call(run_em_command)
        # STEP 9: This creates input parameters for NVT molecular dynamics.
        if not os.path.exists(f'nvt.tpr'):
            nvt_md_command = ['gmx_mpi_d', 'grompp',
                              '-f', f'{self.WORKING_DIR}/nvt.mdp',
                              '-c', 'em.gro',
                              '-r', 'em.gro',
                              '-p', 'topol.top',
                              '-o', 'nvt.tpr']
            self.suprocess_call(nvt_md_command)
        # STEP 10: This runs the NVT MD.
        tic = time.time()
        run_nvt_md_command = ['gmx_mpi_d', 'mdrun', '-deffnm', 'nvt']
        self.suprocess_call(run_nvt_md_command)
        tac = time.time()
        print(f'NVT MD {identifier} - took {round((tac-tic)/60,2)} minutes')
        tic = time.time()
        # STEP 11: This creates input parameters for NPT molecular dynamics.
        if not os.path.exists(f'npt.tpr'):
            npt_md_command = ['gmx_mpi_d', 'grompp',
                              '-f', f'{self.WORKING_DIR}/npt.mdp',
                              '-c', 'nvt.gro',
                              '-r', 'nvt.gro',
                              '-t', 'nvt.cpt',
                              '-p', 'topol.top',
                              '-o', 'npt.tpr']
            self.suprocess_call(npt_md_command)
        # STEP 12: This runs the NPT MD for 100ps.
        run_npt_md_command = ['gmx_mpi_d', 'mdrun', '-deffnm', 'npt']
        self.suprocess_call(run_npt_md_command)
        tac = time.time()
        print(f'NPT MD {identifier} - took {round((tac-tic)/60,2)} minutes')
        print(' ')
        
        
    def thread_simulations(self, args_list, threads=6):
        args_list = [identifier for identifier in args_list 
                     if not os.path.exists(f'{self.DUMP_DIR}/{identifier}')]  
        if args_list:
            batches_list = [i for i in self.batch(args_list, 10)]
            for batch in batches_list:
                children = [Process(target=self.run_protocol, args=(identifier,)) for identifier in batch]
                for child in children:
                    child.start()
                for child in children:
                    child.join()
                    
                    