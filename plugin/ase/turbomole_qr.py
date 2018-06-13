from __future__ import print_function
"""This module defines a QR-specific ASE interface to Turbomole

command line define:
 set "basis=cefine" -> set "method=<cefine command>" 
eg: qr.refine [...] quatum.method='cefine -func pbe0 -fon -bas minix -ri -noopt -d3 ' quantum.basis='cefine'


http://www.turbomole.com/
"""
import os
import sys

import numpy as np
from ase.units import mol, kcal
from ase.units import Hartree, Bohr
from ase.io import read, write
from ase.calculators.general import Calculator
from subprocess import Popen, PIPE, STDOUT
import copy

key_parameters = {'maxit': 200}
#                   'basis': None,
#                   'coordinates': None,
#                   'charge': None,
#                   'method': None,
#                   'dftd': None,
#                   'pointcharges': None,}

class Turbomole(Calculator):
    def __init__(self, label='turbomole',
                 calculate_energy='ridft',
                 calculate_forces='rdgrad',
                 charge='0',
                 method='b-p',
                 basis='def2-SV(P)',
                 define_str =   '\n\na coord\n*\nno\nb all def-SV(P)\n*\neht\n\n'+str(0)+'\n\nscf\niter\n300\n\ncc\nmemory\n4000\n*\ndft\non\nfunc\nb-p\n*\nri\non\nm\n1000\n*\n* ',
                 post_HF=False,
                 pointcharges = None):

        self.label = label
        self.converged = False
         # set calculators for energy and forces
        self.calculate_energy = calculate_energy
        self.calculate_forces = calculate_forces

        self.key_parameters = copy.deepcopy(key_parameters)
        self.key_parameters['basis'] = basis
        self.key_parameters['method'] = method
        self.key_parameters['charge'] = charge
        self.command = calculate_energy
        # turbomole has no stress
        self.stress = np.empty(6)

        # storage for energy and forces
        self.e_total = None
        self.energy_free = None
        self.forces = None
        self.updated = False

        # atoms must be set
        self.atoms = None

        # POST-HF method
        self.post_HF = post_HF
        self.pointcharges = pointcharges

    def initialize(self, atoms):
        self.numbers = atoms.get_atomic_numbers().copy()
        self.species = []
        for a, Z in enumerate(self.numbers):
            self.species.append(Z)
        self.converged = False

    def execute(self, command):

        try:
            if self.pointcharges is not None:
              f = open(self.pointcharges, "r")
              point_charges = f.readlines()
              f.close()
              f = open("control", "r")
              contents = f.readlines()
              f.close()
              contents = contents[:-1] + ['$point_charges\n'] + point_charges +[contents[-1]]
              f = open("control", "w")
              f.writelines( contents )
              f.close()
            # the sub process gets started here
            proc = Popen([command], shell=True, stderr=PIPE)
            error = proc.communicate()[1]
            # check the error output
            if 'abnormally' in error:
                raise OSError(error)
            #print('TM command: ', command, 'successfully executed')
        except OSError as e:
            print('Execution failed:', e, file=sys.stderr)
            sys.exit(1)

    def get_potential_energy(self, atoms):
        # update atoms
        self.updated = self.e_total is None
        self.set_atoms(atoms)
        # if update of energy is necessary
        if self.update_energy:
            # calculate energy
            self.execute(self.calculate_energy + ' > ASE.TM.energy.out')
            # check for convergence of dscf cycle
            if os.path.isfile('dscf_problem'):
                print('Turbomole scf energy calculation did not converge')
                raise RuntimeError(
                    'Please run Turbomole define and come thereafter back')
            # read energy
            self.read_energy()
        self.update_energy = False
        return self.e_total


    def get_forces(self, atoms):
        # update atoms
        self.updated = self.forces is None
        self.set_atoms(atoms)
        # complete energy calculations
        if self.update_energy:
            self.get_potential_energy(atoms)
        # if update of forces is necessary
        if self.update_forces:
            # calculate forces
            self.execute(self.calculate_forces + ' > ASE.TM.forces.out')
            # read forces
            self.read_forces()
        self.update_forces = False
        return self.forces.copy()



    def get_stress(self, atoms):
        return self.stress

    def set_atoms(self, atoms):
        # Delete old  coord control, ... files, if exist
        for f in ['coord',
          'basis',
          'energy',
          'gradients',
          'alpha',
          'beta',
          'mos',
          'forceapprox',
          'statistics',
          'dscf_problem',
          'control']:
                if os.path.exists(f):
                        os.remove(f)
#        if self.atoms == atoms:
#            if (self.updated and os.path.isfile('coord')):
#                self.updated = False
#                a = read('coord').get_positions()
#                if np.allclose(a, atoms.get_positions(), rtol=0, atol=1e-13):
#                    return
#            else:
#                return
        # performs an update of the atoms
        write('coord', atoms)

        if not self.key_parameters["basis"] == 'cefine':
            string='\n\na coord\n*\nno\nb all '+self.key_parameters['basis']+'\n*\neht\n\n'+str(self.key_parameters['charge'])+'\n\nscf\niter\n'+str(self.key_parameters['maxit'])+'\n\ncc\nmemory\n1000\n*\ndft\non\nfunc\n'+self.key_parameters['method']+'\n*\nri\non\nm\n3000\n*\n*' 
            self.define_str = string
            with open('def.inp', 'w') as f:
                f.write(self.define_str)
            command = 'define < def.inp > define.out'
        else:
            self.define_str = self.key_parameters['method']+'-chrg '+str(self.key_parameters['charge'])
            command = self.define_str+' > cefine.out'
            
        # run define
        proc = Popen([command], shell=True, stderr=PIPE)
        error = proc.communicate()[1]
        if 'abnormally' in error:
            raise OSError(error)
        Calculator.set_atoms(self, atoms)
        # energy and forces must be re-calculated
        self.update_energy = True  
        self.update_forces = True


    def read_energy(self):
        """Read Energy from Turbomole energy file."""
        text = open('energy', 'r').read().lower()
        lines = iter(text.split('\n'))

        # Energy:
        for line in lines:
            if line.startswith('$end'):
                break
            elif line.startswith('$'):
                pass
            else:
                energy_tmp = float(line.split()[1])
                if self.post_HF:
                    energy_tmp += float(line.split()[4])
        # update energy units
        # print("E(TM): "+str(energy_tmp)) #debug
        self.e_au=energy_tmp
        self.e_total = energy_tmp * Hartree/(kcal / mol)
        self.energy_free = self.e_total

    def read_forces(self):
        """Read Forces from Turbomole gradient file."""
        file = open('gradient', 'r')
        lines = file.readlines()
        file.close()

        forces = np.array([[0, 0, 0]])

        nline = len(lines)
        iline = -1

        for i in range(nline):
            if 'cycle' in lines[i]:
                iline = i

        if iline < 0:
            raise RuntimeError('Please check TURBOMOLE gradients')

        # next line
        iline += len(self.atoms) + 1
        # $end line
        nline -= 1
        # read gradients
        for i in range(iline, nline):
            line = lines[i].replace('D', 'E')
            if "***" in line:  # is this a real thing?
                 raise RuntimeError('Please check TURBOMOLE gradients')  
            tmp = np.array([[float(f) for f in line.split()[0:3]]])
            forces = np.concatenate((forces, tmp))
        # Note the '-' sign for turbomole, to get forces
        self.forces = (-np.delete(forces, np.s_[0:1], axis=0)) * (Hartree / Bohr)/(kcal / mol)

    def calculation_required(self, atoms, properties):
        if self.atoms != atoms:
            return True
        for prop in properties:
            if prop == 'energy' and self.e_total is None:
                return True
            elif prop == 'forces' and self.forces is None:
                return True
        return False

    def set_modules(self):
        with open('control', 'r') as out:
            for line in out:
                if '$rij' in line:
                    self.calculate_energy='ridft'
                    self.calculate_forces='rdgrad'
                    break
                else:
                    self.calculate_energy='dscf'
                    self.calculate_forces='grad'


    def update(self, atoms_new):
        if not self.atoms_are_equal(atoms_new):
            self.atoms = atoms_new.copy()
            self.qr_run()

    def get_command(self):
        """Return command string if program installed, otherwise None.  """
        command = None
        if self.command is not None:
          command = self.command
        elif ('TM_COMMAND' in os.environ):
          command = os.environ['TM_COMMAND']
        return command

    # Q|R requirements
    def set_charge(self, charge):
      self.key_parameters['charge'] = charge

    def set_basis(self, basis):
      self.key_parameters['basis'] = basis

    def set_method(self, method):
      self.key_parameters['method'] = method

    def set_label(self, label):
      self.label = label

    def set(self, **kwargs):
        for key, value in kwargs.items():
           if key in key_parameters:
                self.key_parameters[str(key)]=value

    def run_qr(self,
               atoms,
               define_str,
               coordinates,
               charge,
               pointcharges,
               command=None):
        import subprocess
        self.atoms = atoms
        self.key_parameters['charge'] = charge
        
        self.coordinates = coordinates
        working_dir = os.getcwd()
        if self.pointcharges is not None:
          self.pointcharges = os.path.abspath(self.pointcharges)
        if not  os.path.isdir(self.label):
          os.mkdir(self.label)
        turbomole_dir = self.label
        os.chdir(turbomole_dir)
        self.set_atoms(self.atoms)
        try:
            if self.pointcharges is not None:
              f = open(self.pointcharges, "r")
              point_charges = f.readlines()
              f.close()
              f = open("control", "r")
              contents = f.readlines()
              f.close()
              contents = contents[:-1] + ['$point_charges\n'] + point_charges +[contents[-1]]
              f = open("control", "w")
              f.writelines( contents )
              f.close()

            self.set_modules() # ridft or dscf  
            command = self.calculate_energy + ' > ASE.TM.energy.out'
            print(command) #debug
            proc = subprocess.Popen([command], shell=True, stderr=PIPE)
            error = proc.communicate()[1]
            proc.wait()
            exitcode = proc.returncode
        
            
            if exitcode != 0:
                print(str(exitcode))
                raise RuntimeError('Turbomole exited with error code')
            if 'abnormally' in error:
                raise OSError(error)
            # check for convergence of dscf cycle
            if os.path.isfile('dscf_problem'):
                print('Turbomole scf energy calculation did not converge')
                raise RuntimeError(
                'Please run Turbomole define and come thereafter back')
            # read energy
            self.read_energy()
            print(self.label+':  '+str(self.e_au)) #debug

            # calculate forces
            command = self.calculate_forces + ' > ASE.TM.forces.out'
            proc = subprocess.Popen([command], shell=True, stderr=PIPE)
            error = proc.communicate()[1]
            proc.wait()
            exitcode = proc.returncode
            if exitcode != 0:
                print(str(exitcode))
                raise RuntimeError('Turbomole exited with error code')
            if 'abnormally' in error:
                raise OSError(error)
            # read forces
            self.read_forces()
        except OSError as e:
            print('Execution failed:', e, file=sys.stderr)
            sys.exit(1)
        os.chdir(working_dir)
