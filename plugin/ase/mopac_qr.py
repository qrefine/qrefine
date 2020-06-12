"""
Version 2012/08/20, Torsten Kerber
Contributors:
  Torsten Kerber, Ecole normale superieure de Lyon:
  Paul Fleurat-Lessard, Ecole normale superieure de Lyon
  based on a script by Rosa Bulo, Ecole normale superieure de Lyon
This work is supported by Award No. UK-C0017, made by King Abdullah
University of Science and Technology (KAUST), Saudi Arabia
See accompanying license files for details.
"""
from __future__ import print_function
import os
import string
import numpy as np
import platform

from ase.units import kcal, mol
from ase.calculators.general import Calculator

str_keys = ['functional', 'job_type']
int_keys = ['restart', 'spin', 'charge']
bool_keys = ['OPT']
float_keys = ['RELSCF']


class Mopac(Calculator):
    name = 'MOPAC'
    def __init__(self,
                nproc=1,
                 label='ase',
                 **kwargs):
        # define parameter fields
        self.str_params = {}
        self.int_params = {}
        self.bool_params = {}
        self.float_params = {}

        # initials parameter fields
        for key in str_keys:
            self.str_params[key] = None
        for key in int_keys:
            self.int_params[key] = None
        for key in bool_keys:
            self.bool_params[key] = None
        for key in float_keys:
            self.float_params[key] = None

        # set initial values
        functional = 'PM7'
        env = os.environ.get('MOPAC_FUNCTIONAL', None)
        if env: functional = env
        self.set(restart=0,
                 spin=0,
                 OPT=False,
                 functional=functional,
                 job_type=' 1SCF  GRADIENTS AUX(0,PRECISION=9) ',
                 RELSCF= None)
        # set user values
        self.set(**kwargs)

        # save label
        self.label = label

        #set atoms
        self.atoms = None
        # initialize the results
        self.version = None
        self.energy_zero = None
        self.energy_free = None
        self.forces = None
        self.stress = None
        self.calc_dir = None
        # initialize the results
        self.occupations = None

        # command
        self.command = self.get_command()

    def set(self, **kwargs):
        """
        Sets the parameters on the according keywords
        Raises RuntimeError when wrong keyword is provided
        """
        for key in kwargs:
            if key in self.bool_params:
                self.bool_params[key] = kwargs[key]
            elif key in self.int_params:
                self.int_params[key] = kwargs[key]
            elif key in self.str_params:
                self.str_params[key] = kwargs[key]
            elif key in self.float_params:
                self.float_params[key] = kwargs[key]
            else:
                raise RuntimeError('MOPAC calculator: unknown keyword: ' + key)

    def get_version(self):
        return self.version

    def initialize(self, atoms):
        pass

    def write_input(self, fname, atoms):
        """
        Writes the files that have to be written each timestep
        """
        # start the input
        mopac_input = ''

        #write functional and job_type
        for key in 'functional', 'job_type':
            if self.str_params[key] != None:
                mopac_input += self.str_params[key] + ' '

        if self.float_params['RELSCF'] != None:
            mopac_input += 'RELSCF=' + str(self.float_params['RELSCF']) + ' '

        #write charge/
        # charge = sum(atoms.get_initial_charges())
        #if charge != 0:
        #   mopac_input += 'CHARGE=%i ' % (charge)
        charge=self.int_params['charge']
        mopac_input += 'CHARGE= ' + str(charge)+'  '

        if (self.int_params['nproc'] > 1):
            nproc=self.int_params['nproc']
        else:
            nproc=1

        # threads should be specified by user
        mopac_input += ' THREADS=%i' %(nproc)

        # add solvent
        mopac_input += '  EPS=78.4'

        #write spin
        spin = self.int_params['spin']
        if spin == 1.:
            mopac_input += 'DOUBLET '
        elif spin == 2.:
            mopac_input += 'TRIPLET '

        #input down
        mopac_input += '\n'
        mopac_input += 'Title: ASE job\n\n'

        f = 1
        # write coordinates
        for iat in range(len(atoms)):
            atom = atoms[iat]
            xyz = atom.position
            mopac_input += ' %2s' % atom.symbol
            # write x, y, z
            for idir in range(3):
                mopac_input += '    %16.5f %i' % (xyz[idir], f)
            mopac_input += '\n'

        if atoms.pbc.any():
            for v in atoms.get_cell():
                mopac_input += 'Tv %8.3f %8.3f %8.3f\n' % (v[0], v[1], v[2])

        # write input
        myfile = open(fname, 'w')
        myfile.write(mopac_input)
        myfile.close()
        self.mopac_input = mopac_input

    def get_command(self):
        """Return command string if program installed, otherwise None.  """
        command = None
        if ('MOPAC_COMMAND' in os.environ):
          command = os.environ['MOPAC_COMMAND']
        return command

    def set_command(self, command):
      self.command = command


    def run_command(self,command):
        """
        execute <command> in a subprocess and check error code
        """
        from subprocess import Popen, PIPE, STDOUT
        if command == '':
            raise RuntimeError('no command for run_command :(')
        # print 'Running: ', command #debug
        proc = Popen([command], shell=True, stderr=PIPE)
        proc.wait()
        exitcode = proc.returncode
        if exitcode != 0:
            # print exitcode,'label:', self.calc_dir
            error='%s exited with error code %i in %s' % (
                           command,exitcode,self.calc_dir)
            stdout,stderr = proc.communicate()
            print('shell output: ',stdout,stderr)
            raise RuntimeError(error)
        return 0

    def run(self):
        import subprocess, shlex
        from threading import Timer

        def run_timeout(cmd, timeout_sec):
          proc = subprocess.Popen(shlex.split(cmd),
          # proc = subprocess.Popen(cmd, 
                                  stdout=subprocess.PIPE,
                                  shell=True,
                                  stderr=subprocess.PIPE)
          kill_proc = lambda p: p.kill()
          timer = Timer(timeout_sec, kill_proc, [proc])
          try:
            timer.start()
            stdout,stderr = proc.communicate()
            print(stdout,stderr)
          finally:
            timer.cancel()

        """
        Writes input in label.mop
        Runs MOPAC
        Reads Version, Energy and Forces
        """
        # set the input file name
        finput = self.label + '.mop'
        foutput = self.label + '.out'
        self.write_input(finput, self.atoms)

         # directory
        self.calc_dir = os.getcwd()

        command = self.command
        if command is None:
          raise RuntimeError('MOPAC_COMMAND is not specified')

        WhatOS=platform.system()
        if "Linux" in WhatOS:
            if ('MOPAC_DIR' in os.environ):
                mdir = os.environ['MOPAC_DIR']
            else:
                raise RuntimeError('MOPAC_DIR is not specified')
            command_exc= "LD_PRELOAD=%s/libiomp5.so %s  %s" % (mdir,command,finput)
        if "Darwin" in WhatOS:
            command_exc= "  ".join([command , finput])

        # run_timeout(command_exc ,72000)# 20hours
        self.run_command(command_exc)
#        exitcode = os.system('%s %s' % (command, finput)+ '  > /dev/null 2>&1    ')

#        if exitcode != 0:
#            raise RuntimeError('MOPAC exited with error code')

        self.version = self.read_version(foutput)
        energy = self.read_energy(foutput)
        self.energy_zero = energy
        self.energy_free = energy
        self.forces = self.read_forces(foutput)

    def read_version(self, fname):
        """
        Reads the MOPAC version string from the second line
        """
        version = 'unknown'
        lines = open(fname).readlines()
        for line in lines:
            if "  Version" in line:
                version = line.split()[-2]
                break
        return version

    def read_energy(self, fname):
        """
        Reads the ENERGY from the output file (HEAT of FORMATION in kcal / mol)
        Raises RuntimeError if no energy was found
        """
        outfile = open(fname)
        lines = outfile.readlines()
        outfile.close()

        energy = None
        for line in lines:
            if line.find('HEAT OF FORMATION') != -1:
                words = line.split()
                energy = float(words[5])
            if line.find('H.o.F. per unit cell') != -1:
                words = line.split()
                energy = float(words[5])
            if line.find('UNABLE TO ACHIEVE SELF-CONSISTENCE') != -1:
                energy = None
        if energy is None:
            raise RuntimeError('MOPAC: could not find total energy')
### do not change unit for mopac
        energy *= (kcal / mol)
        return energy

    def read_forces(self, fname):
        """
        Reads the FORCES from the output file
        search string: (HEAT of FORMATION in kcal / mol / AA)
        """
        outfile = open(fname)
        lines = outfile.readlines()
        outfile.close()
        nats = len(self.atoms)
        forces = np.zeros((nats, 3), float)
        infinite_force="*****"
        if 'mozyme' in self.str_params['job_type'].lower():
            for i, line in enumerate(lines):
                if line.find('FINAL  POINT  AND  DERIVATIVES') != -1:
                    for j in range(nats):
                        gline = lines[i + j + 5]
                        pre_force=gline[8:35]
                        if(infinite_force in pre_force):
                            forces[j] = [999999999.9999,999999999.9999,999999999.9999]
                        else:
                            forces[j] =  [float( pre_force[0:9].strip()),float( pre_force[9:18].strip()),float( pre_force[18:27].strip())]
        else:
          for i, line in enumerate(lines):
            if line.find('GRADIENT\n') != -1:
                for j in range(nats * 3):
                    gline = lines[i + j + 1]
                    pre_force=gline[49:62]
                    if(infinite_force in pre_force):
                        forces[int(j/3), int(j%3)] =999999999.9999
                    else:
                        forces[int(j/3), int(j%3)] = float(pre_force)
                break
#do not change unit for mopac
        forces *= - (kcal / mol)
        return forces

    def atoms_are_equal(self, atoms_new):
        ''' (adopted from jacapo.py)
        comparison of atoms to self.atoms using tolerances to account
        for float/double differences and float math.
        '''

        TOL = 1.0e-6  # angstroms

        # check for change in cell parameters
        test = len(atoms_new) == len(self.atoms)
        if test is not True:
            return False

        # check for change in cell parameters
        test = (abs(self.atoms.get_cell() - atoms_new.get_cell()) <= TOL).all()
        if test is not True:
            return False

        old = self.atoms.arrays
        new = atoms_new.arrays

        # check for change in atom position
        test = (abs(new['positions'] - old['positions']) <= TOL).all()
        if test is not True:
            return False

        # passed all tests
        return True

    def update(self, atoms_new, **kwargs):
        self.set(**kwargs)
        if not self.atoms_are_equal(atoms_new):
            self.atoms = atoms_new.copy()
            self.run()

    def run_qr(self, atoms_new, **kwargs):
        for key in kwargs:
          if key in self.bool_params:
              self.bool_params[key] = kwargs[key]
          elif key in self.int_params:
              self.int_params[key] = kwargs[key]
          elif key in self.str_params:
              self.str_params[key] = kwargs[key]
          elif key in self.float_params:
              self.float_params[key] = kwargs[key]
        self.atoms = atoms_new.copy()
        self.run()

    # Q|R requirements
    def set_charge(self, charge):
      self.int_params['charge'] = charge

    def set_method(self, method):
      self.str_params['functional'] = method

    def set_label(self, label):
      self.label = label

    def set_nproc(self, nproc):
      self.int_params['nproc'] = int(nproc)
