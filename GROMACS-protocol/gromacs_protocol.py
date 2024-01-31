# Before running:
# a) Set up node by install nvm. The node version should be around 16 or later
# b) Install npm icn3d by following the "Installation" instruction at 
# https://github.com/ncbi/icn3d/tree/master/icn3dnode
# c) Download the script "rmhet.js" and "addmissingatoms.js" to your directory from
# https://github.com/ncbi/icn3d/tree/master/icn3dnode.

import os
import sys
import time
import shutil
from subprocess import PIPE, DEVNULL, STDOUT, Popen, run, TimeoutExpired, call
from multiprocessing import Process


class GromacsProtocol:
    def __init__(
            self,
            pdb_directory: str,
            output_directory: str,
            mdp_directory: str,
            ntmpi: str,
            ntomp: str,
            GMX='gmx',
            hard_force: bool = False,
    ):
        self.pdb_directory = pdb_directory
        self.output_directory = output_directory
        self.mdp_directory = mdp_directory
        self.GMX = GMX
        self.ntmpi = ntmpi
        self.ntomp = ntomp
        self.hard_force = hard_force
        self.identifiers_list = [
            filename.rsplit('_NoHOH.pdb')[0]
            if "NoHOH" in filename
            else filename.rsplit('.pdb')[0]
            for filename in os.listdir(self.pdb_directory)
        ]
        if self.hard_force:
            if os.path.exists(self.output_directory):
                shutil.rmtree(self.output_directory)
        if not os.path.exists(self.output_directory):
            os.mkdir(self.output_directory)
        self.files_req = [
            os.path.join(self.mdp_directory, 'ions.mdp'),
            os.path.join(self.mdp_directory, 'minim.mdp'),
            os.path.join(self.mdp_directory, 'nvt.mdp'),
            os.path.join(self.mdp_directory, 'npt.mdp')
        ]
        if not all([os.path.exists(file_path) for file_path in self.files_req]):
            raise ValueError(
                f'files: ions.mdp, minim.mdp,'
                f' nvt.mdp,npt.mdp should be in'
                f' mdp directory: {self.mdp_directory}'
            )
        if len(self.identifiers_list) != len(set(self.identifiers_list)):
            raise ValueError(
                'all identifiers should be distinct'
            )
        # remove incomplete runs
        for identifier in os.listdir(self.output_directory):
            dir_path = os.path.join(
                self.output_directory,
                identifier
            )
            contains_edr = any(
                file.endswith('.edr') for file in os.listdir(dir_path)
            )
            if not contains_edr:
                shutil.rmtree(dir_path)

    @staticmethod
    def subprocess_call(
            command_list: list,
    ):
        try:
            # Run the command, redirect stdout to DEVNULL, capture stderr
            result = run(
                command_list,
                stdout=DEVNULL,
                stderr=PIPE
            )
            # Check the return code of the command
            if result.returncode != 0:
                command_string = ' '.join(command_list)
                # Decode stderr to convert from bytes to string
                error_message = result.stderr.decode('utf-8') if result.stderr else 'Unknown error'
                raise ValueError(f'Error while calling: {command_string}\nError Message: {error_message}')
        except Exception as e:
            print(e)

    def protocol(
            self,
            identifier: str,
    ):
        simulation_directory = os.path.join(self.output_directory, identifier)
        if not os.path.exists(simulation_directory):
            os.mkdir(simulation_directory)
        os.chdir(simulation_directory)
        # STEP 0:remove any water molecules to create a clean PDB:
        pdb_file_direct = os.path.join(self.pdb_directory, f'{identifier}_NoHOH.pdb')
        if not os.path.exists(pdb_file_direct):
            pdb_file_direct = os.path.join(self.pdb_directory, f'{identifier}.pdb')
            molecules = ["HOH", "PO4"]
            pattern = '|'.join(molecules)
            rmv_hoh_command = [
                "grep",
                "-v",
                "-E",  # Use extended regular expressions
                pattern,
                pdb_file_direct
            ]
            with open(f"{identifier}_NoHOH.pdb", 'w') as out:
                run(rmv_hoh_command, stdout=out)
        # STEP 1: The following step will remove all hetero atoms,
        # use one chain name, and add missing atoms
        rmhet_command = [
            "node",
            os.path.join(self.mdp_directory, 'rmhet.js'),
            os.path.join(self.pdb_directory, f'{identifier}_NoHOH.pdb')
        ]
        with open(f'{identifier}_nohet.pdb', 'w') as out:
            run(rmhet_command, stdout=out)
        addmissingatom_command = [
            "node",
            os.path.join(self.mdp_directory, 'addmissingatoms.js'),
            f'{identifier}_nohet.pdb'
        ]
        with open(f'{identifier}_clean.pdb', 'w') as out:
            run(addmissingatom_command, stdout=out)

        # STEP 2: Here we select AMBER99SB-ILDN force-field
        # (#6 in the list)
        # for describing atom-atom interaction energies
        # to be used in the GROMACS simulations.
        if not os.path.exists(f'{identifier}_processed.gro'):
            pdb2gmx_command = [
                self.GMX,
                'pdb2gmx',
                '-f',
                os.path.join(self.pdb_directory, f'{identifier}_clean.pdb'),
                '-o',
                f'{identifier}_processed.gro',
                '-water',
                'spce',
                '-ff',
                'amber99sb-ildn'
            ]
            self.subprocess_call(pdb2gmx_command)
        # STEP 3: Here we create a simulation box 0.5nm
        # from the protein edge in all 6 directions.
        if not os.path.exists(f'{identifier}_newbox.gro'):
            editconf_command = [
                self.GMX,
                'editconf',
                '-f',
                f'{identifier}_processed.gro',
                '-o',
                f'{identifier}_newbox.gro',
                '-c',
                '-d',
                '0.5',
                '-bt',
                'cubic'
            ]
            self.subprocess_call(editconf_command)
        # STEP 4: Here we add water molecules in the
        # simulation box and new topology file is written.
        if not os.path.exists(f'{identifier}_solv.gro'):
            solvate_command = [
                self.GMX,
                'solvate',
                '-cp',
                f'{identifier}_newbox.gro',
                '-cs',
                'spc216.gro',
                '-o',
                f'{identifier}_solv.gro',
                '-p',
                f'topol.top'
            ]
            self.subprocess_call(solvate_command)
        # STEP 5: Generate parameter file ions.tpr for all atoms.
        if not os.path.exists(f'ions.tpr'):
            grompp_command = [
                self.GMX,
                'grompp',
                '-f',
                os.path.join(self.mdp_directory, "ions.mdp"),
                '-c',
                f'{identifier}_solv.gro',
                '-p',
                'topol.top',
                '-o',
                'ions.tpr',
                '-maxwarn',
                '2'
            ]
            self.subprocess_call(grompp_command)
        # STEP 6: This will add ions to your simulation box and replace
        # the solvent molecules if there is a clash.
        if not os.path.exists(f'{identifier}_solv_ions.gro'):
            genion_command = [
                'echo',
                'SOL',
                '|',
                self.GMX,
                'genion',
                '-s',
                f'ions.tpr',
                '-o',
                f'{identifier}_solv_ions.gro',
                '-p',
                f'topol.top',
                '-pname',
                'NA',
                '-nname',
                'CL',
                '-neutral',
                '> /dev/null 2>&1'  # supress stdout
            ]
            os.system(' '.join(genion_command))
            # STEP 7: Energy minimization
        if not os.path.exists('em.tpr'):
            em_command = [
                self.GMX,
                'grompp',
                '-f',
                f'{self.mdp_directory}/minim.mdp',
                '-c',
                f'{identifier}_solv_ions.gro',
                '-p',
                'topol.top',
                '-o',
                'em.tpr'
            ]
            self.subprocess_call(em_command)
        # STEP 8: This runs the energy minimization.
        run_em_command = [
            self.GMX,
            'mdrun',
            '-ntomp',
            self.ntomp,
            '-deffnm',
            'em'
        ]
        if self.ntmpi:
            run_em_command.extend(['-ntmpi', str(self.ntmpi)])
        self.subprocess_call(run_em_command)
        # STEP 9: This creates input parameters
        # for NVT molecular dynamics.
        if not os.path.exists(f'nvt.tpr'):
            nvt_md_command = [
                self.GMX,
                'grompp',
                '-f',
                f'{self.mdp_directory}/nvt.mdp',
                '-c',
                'em.gro',
                '-r',
                'em.gro',
                '-p',
                'topol.top',
                '-o',
                'nvt.tpr'
            ]
            self.subprocess_call(nvt_md_command)
        # STEP 10: This runs the NVT MD.
        tic = time.time()
        run_nvt_md_command = [
            self.GMX,
            'mdrun',
            '-ntomp',
            self.ntomp,
            '-deffnm',
            'nvt'
        ]
        if self.ntmpi:
            run_nvt_md_command.extend(['-ntmpi', str(self.ntmpi)])
        self.subprocess_call(run_nvt_md_command)
        tac = time.time()
        sys.stdout.write(
            f'NVT MD {identifier} '
            f'- took {round((tac - tic) / 60, 2)} minutes \n'
        )
        tic = time.time()
        # STEP 11: This creates input parameters for NPT molecular dynamics.
        if not os.path.exists(f'npt.tpr'):
            npt_md_command = [
                self.GMX,
                'grompp',
                '-f',
                f'{self.mdp_directory}/npt.mdp',
                '-c',
                'nvt.gro',
                '-r',
                'nvt.gro',
                '-t',
                'nvt.cpt',
                '-p',
                'topol.top',
                '-o',
                'npt.tpr'
            ]
            self.subprocess_call(npt_md_command)
        # STEP 12: This runs the NPT MD for 100ps.
        run_npt_md_command = [
            self.GMX,
            'mdrun',
            '-ntomp',
            self.ntomp,
            '-deffnm',
            'npt'
        ]
        if self.ntmpi:
            run_npt_md_command.extend(['-ntmpi', str(self.ntmpi)])
        self.subprocess_call(run_npt_md_command)
        tac = time.time()
        sys.stdout.write(
            f'NPT MD {identifier} '
            f'- took {round((tac - tic) / 60, 2)} minutes \n'
        )
        sys.stdout.write(' ')

    def main(
            self,
            args_list
    ):
        args_list = [
            identifier
            for identifier in args_list
            if not os.path.exists(
                os.path.join(self.output_directory, identifier)
            )
        ]
        if args_list:
            for arg in args_list:
                try:
                    self.protocol(
                        identifier=arg
                    )
                except ValueError as e:
                    sys.stdout.write(
                        f'**** ERROR PROCESSING: {arg} ****\n'
                        f'{e}\n'
                    )
                    continue
