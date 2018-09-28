"""
  based on ASE script for Mopac and then for orca.


"""
import os
import string
import numpy as np
from ase.io import write
from ase.units import kcal, mol
from ase.units import Hartree, Bohr
from ase.calculators.general import Calculator

import copy

key_parameters = {
                  'coordinates': None,
                  'charge': None,
                  'method': None,
                  'pointcharges': None,
                  }

class GFNxTB(Calculator):
    name = 'gfn-xtb'

    def __init__(self,
                 label='ase_plugin',
                 coordinates='xtb_tmp.xyz',
                 charge='0',
                 version=2,
                 method='-gfn2',
                 atoms=None,
                 command=None,
                 **kwargs):

        self.coordinates = coordinates
        self.key_parameters = copy.deepcopy(key_parameters)
        self.key_parameters['coordinates'] = coordinates
        self.key_parameters['charge'] = charge
        self.key_parameters['method'] = method
        self.version=version
        # save label
        self.label = label
        # set atoms
        self.atoms = atoms
        # initialize the results
        self.energy_zero = None
        self.energy_free = None
        self.forces = None
        self.stress = None
        self.command = command


    def run_command(self,command):
        """
        execute <command> in a subprocess and check error code
        """
        from subprocess import Popen, PIPE, STDOUT
        if command == '':
            raise RuntimeError('no command for run_command :(')
        print 'Running: ', command #debug
        proc = Popen([command], shell=True, stderr=PIPE)
        proc.wait()
        exitcode = proc.returncode
        if exitcode != 0:
            print exitcode
            raise RuntimeError(command+' exited with error code')
        return 0

    def write_input(self,atoms):
        
        atoms = copy.deepcopy(self.atoms)
        write(self.coordinates, atoms)
        fname=self.coordinates #+'/xtb_tmp.xyz'
        finput = open(fname, "w")
        finput.write(" %s \n \n" % len(atoms))
        for index in range(len(atoms)):
            finput.write(str(atoms.get_chemical_symbols()[index]) +  " "
                            + str(atoms[index].position[0])+  " "
                            + str(atoms[index].position[1])+  " "
                            + str(atoms[index].position[2])+ "\n" )
        finput.close()

    def get_command(self):
        """Return command string if program installed, otherwise None.  """
        command = None
        if self.command is not None:
          command = self.command
        elif ('XTBHOME' in os.environ):
          command = os.environ['XTBHOME']+'/xtb '
        return command

    def run_qr(self,
               atoms,   
               coordinates,
               charge,
               pointcharges,
               command=None,
               define_str=None
        ):
        import subprocess
        """
        Handels GFN-xTB calculations
        """
        # set the input file name
        self.atoms = atoms
        method=self.key_parameters['method']
        self.coordinates = coordinates
        self.key_parameters['charge'] = charge
        foutput = self.label + '.out'
        
        self.coordinates = coordinates
        working_dir = os.getcwd()
        if not  os.path.isdir(self.label):
          os.mkdir(self.label)
        calc_dir = os.path.join(working_dir,self.label)
        os.chdir(calc_dir)
        self.coordinates = 'xtb_tmp.xyz'
        print 'current work dir ',os.getcwd()
        # print('coords:',coordinates)
        self.write_input(self.atoms)


        command = self.get_command()
        command=command +str(self.coordinates)+' -chrg '+str(self.key_parameters['charge'])+' -grad '+str(method)+' > xtb.out'
        if command is None:
            raise RuntimeError('$XTBHOME not set')

        self.run_command(command)
            
        self.read_energy()
        self.read_forces()
        self.energy_zero= self.energy_free
        # print 'debug xtb energy:',self.energy_zero
        os.chdir(working_dir)

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
        # update energy units
        # print("E(TM): "+str(energy_tmp)) #debug
        self.e_au=energy_tmp
        self.e_total = energy_tmp * Hartree/(kcal / mol)
        self.energy_free = self.e_total


    def read_forces(self):
        """xTB uses turbomole format gradients. We read forces and energy from it"""
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
            tmp = np.array([[float(f) for f in line.split()[0:3]]])
            forces = np.concatenate((forces, tmp))
        # Note the '-' sign for turbomole, to get forces
        self.forces = (-np.delete(forces, np.s_[0:1], axis=0)) * (Hartree / Bohr)/(kcal / mol)

    def set_pointcharges(self):
        if self.pointcharges is not None:
              f = open(self.pointcharges, "r")
              point_charges = f.readlines()
              f.close()
            

    def set(self, **kwargs):
        for key, value in kwargs.items():
            self.key_parameters[str(key)] = value

    # Q|R requirements
    def set_charge(self, charge):
      self.key_parameters['charge'] = charge

    def set_method(self, method):
      self.key_parameters['method'] = method

    def set_label(self, label):
      self.label = label

