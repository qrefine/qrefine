"""
collection of add-ons for QM calculations.
contains:
 - gCP
"""
from __future__ import print_function
import os
import sys
from ase.io import read, write

import numpy as np
from ase.units import mol, kcal
from ase.units import Hartree, Bohr

def run_command(command):
    """
    execute <command> in a subprocess and check error code
    """
    from subprocess import Popen, PIPE, STDOUT
    if command == '':
       raise RuntimeError('no command for run_command :(')
    print('Running: ', command) #debug
    proc = Popen([command], shell=True, stderr=PIPE)
    proc.wait()
    exitcode = proc.returncode
    if exitcode != 0:
        print(exitcode)
        raise RuntimeError(command+' exited with error code')
    return 0

def run_gcp(atoms,level):
    write('gcp_tmp.xyz', atoms)
    outfile='gcp.out'
    avail = ['HF/MINIS', 'DFT/MINIS', 'HF/MINIX', 'DFT/MINIX',
            'HF/SV', 'DFT/SV', 'HF/def2-SV(P)', 'DFT/def2-SV(P)', 'HF/def2-SVP',
            'DFT/def2-SVP', 'HF/DZP', 'DFT/DZP', 'HF/def-TZVP', 'DFT/def-TZVP',
            'HF/def2-TZVP', 'DFT/def2-TZVP', 'HF/631Gd', 'DFT/631Gd',
            'HF/def2-TZVP', 'DFT/def2-TZVP', 'HF/cc-pVDZ', 'DFT/cc-pVDZ',
            'HF/aug-cc-pVDZ', 'DFT/aug-cc-pVDZ', 'DFT/SV(P/h,c)', 'DFT/LANL',
            'DFT/pobTZVP', 'TPSS/def2-SVP', 'PW6B95/def2-SVP', 'hf3c', 'pbeh3c','hf/631g']
    avail = [f.lower() for f in avail]
    if level.lower() not in avail:
             print("Warning: selected gCP level not standard! Beware of what you are doing!") # during development times
           # raise RuntimeError("""%s. gCP level not avaiable:  %r""" % (level.upper(), avail))
    exe='gcp gcp_tmp.xyz -grad -v -l '+level+' > '+outfile
    run_command(command=exe)
    return read_gcp(outfile,atoms)

def read_gcp(outfile,atoms):
    """
    reads gcp output and returns energy+gradient in kcal/mol[/A]
    """
    out = open(outfile, 'r')
    nats = len(atoms)
    gradient = np.zeros((nats,3), float)
    energy = 0.0
    with open(outfile, 'r') as out:
        lines = out.readlines()
        for i, line in enumerate(lines):
            if '  Egcp:' in line:
                sline = lines[i].split()
                energy = float(sline[1])
            if 'gradient: Ggcp' in line:
                for  j in range(nats):
                    nline= lines[1+i+j]
                    tmp = np.array([[float(f) for f in nline.split()[:3]]])
                    gradient[j,:3]=tmp[:3]

                    
    print('E(GCP): ', energy)
    energy*=(Hartree)/(kcal / mol)
    gradient*=(Hartree/Bohr)/(kcal / mol)
    # print gradient
    if energy == 0:
        raise RuntimeError("gcp failed")
    return energy,gradient

def run_dftd3(atoms,level):
    write('dftd3_tmp.tmol', atoms,format='turbomole') # the odd ASE XYZ format does not work
    outfile='dftd3.out'
    exe='dftd3 dftd3_tmp.tmol -grad -func '+level+' > '+outfile
    run_command(command=exe)
    return read_dftd3(outfile,atoms)

def read_dftd3(outfile,atoms):
    """
    reads dftd3 output and returns energy+gradient in kcal/mol[/A]
    """
    out = open(outfile, 'r')
    nats = len(atoms)
    gradient = np.zeros((nats,3), float)
    energy = 0.0
    with open(outfile, 'r') as out:
        lines = out.readlines()
        for i, line in enumerate(lines):
            if 'Edisp /kcal' in line:
                sline = lines[i].split()
                energy = float(sline[3])

    with open('dftd3_gradient','r') as out:
        lines = out.readlines()
        for  j in range(nats):
            nline= lines[j]
            tmp = np.array([[float(f) for f in nline.split()[:3]]])
            gradient[j,:3]=tmp[:3]

                    
    print('E(D3): ', energy)
    energy*=(Hartree)/(kcal / mol)
    gradient*=(Hartree/Bohr)/(kcal / mol)
    # print gradient
    if energy == 0:
        raise RuntimeError("dftd3 failed")
    return energy,gradient

def qm_toolbox(atoms,charge,pointcharges,label,addon,addon_method):

    # mandatory selection of qm_addon_method
    if addon_method is None:
       raise RuntimeError('missing: qm_addon_method !')

    # change to work directory (location of current coordinates)
    cwd = os.getcwd()
    wdir = os.path.join(os.getcwd(),label)
    if not os.path.exists(wdir):
        os.mkdir(wdir)
    os.chdir(wdir)
    
    # select helper program 
    # return E/G in kcal/mol/Angstrom
    if 'gcp-d3' in addon.lower():
        print()
        egcp,ggcp=run_gcp(atoms,addon_method.split("+")[0])
        ed3,gd3=run_dftd3(atoms,addon_method.split("+")[1])
        energy=egcp+ed3
        gradient=ggcp+gd3
    elif 'gcp' in addon.lower():
        energy,gradient=run_gcp(atoms,addon_method)
    elif 'dftd3' in addon.lower():
        energy,gradient=run_dftd3(atoms,addon_method)
    
    else:
       raise RuntimeError('invalid qm_addon')

    os.chdir(cwd)
    return energy, gradient 
