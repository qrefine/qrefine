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
                 nproc='1',
                 atoms=None,
                 command=None,
                 pointcharges=None,
                 **kwargs):

        self.coordinates = coordinates
        self.pointcharges = pointcharges
        self.key_parameters = copy.deepcopy(key_parameters)
        self.key_parameters['coordinates'] = coordinates
        self.key_parameters['charge'] = charge
        self.key_parameters['method'] = method
        self.key_parameters['nproc'] = nproc
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
        self.calc_dir = None


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
            # print exitcode,'label:', self.calc_dir
            error='%s exited with error code %i in %s' % (
                           command,exitcode,self.calc_dir)
            raise RuntimeError(error)
        return 0

    def write_input(self,atoms):
        
        atoms = copy.deepcopy(self.atoms)
        write(self.coordinates, atoms)
        fname=self.coordinates #+'/xtb_tmp.xyz'
        finput = open(fname, "w")        
        symbols = atoms.get_chemical_symbols()
        coordinates = atoms.get_positions()
        finput.write(" %s \n \n" % len(atoms))
        for i in range(len(atoms)):
            finput.write('%-10s' % symbols[i])
            for j in range(3):
                finput.write('%20.10f' % coordinates[i, j])
            finput.write('\n')
        finput.close()

    def get_command(self):
        """Return command string if program installed, otherwise None.  """
        command = None
        if ('XTBHOME' in os.environ):
          command = os.environ['XTBHOME']+'/bin/xtb '
        if command is None:
            raise RuntimeError('$XTBHOME not set')
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

        # directory
        working_dir = os.getcwd()
        if not  os.path.isdir(self.label):
          os.mkdir(self.label)
        self.calc_dir = os.path.join(working_dir,self.label)
        os.chdir(self.calc_dir)
        self.coordinates = 'xtb_tmp.xyz'

        # debug statements
        # print 'current work dir ',os.getcwd()
        # print('coords:',coordinates)
        self.write_input(self.atoms)
        
        #point charges
        self.pointcharges=pointcharges
        if self.pointcharges is not None:
          self.pointcharges = os.path.abspath(self.pointcharges)
          self.set_pointcharges()


        binary = self.get_command()
        if (self.key_parameters['nproc'] > 1):
            nproc=self.key_parameters['nproc']
        else:
            nproc=1


        command='%s %s --chrg %s --grad %s --parallel %s > xtb.out' % (
                binary,
                str(self.coordinates),
                str(self.key_parameters["charge"]),
                str(method),
                str(nproc))

        #clean up
        for f in ['energy','gradient','xtbrestart']:
            if os.path.exists(f):
                os.remove(f)

        self.run_command(command)

        # probably not worth it.
        # if (not self.check_scf_conv):
        #     print 'restarting failed xtb job'
        #     self.run_command(command)

        self.read_energy()
        self.read_forces()
        self.energy_zero= self.energy_free
        os.chdir(working_dir)

    def check_scf_conv(self):
        text = open('energy', 'r').read().lower()
        lines = iter(text.split('\n'))
        scf_conv=False
        for line in lines:
            if 'scf converged':
                scf_conv=True
        return scf_conv

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
            raise RuntimeError('Please check xTB gradients')

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
        f = open(self.pointcharges, "r")
        pchrg = f.readlines()
        f.close()
        f = open("pcharge","w")
        f.writelines(pchrg[0])
        f.writelines(pchrg[2:])
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

    def set_nproc(self, nproc):
      self.key_parameters['nproc'] = str(nproc)