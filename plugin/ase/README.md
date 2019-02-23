# qr-plugin-ase
A plugin for the Atomic Simulation Environment(ASE) (https://wiki.fysik.dtu.dk/ase/)

Quantum refinement code Q|R (https://github.com/qrefine/qrefine/blob/master/restraints.py) calls QM calculator through ASE.

Now Q|R can use Mopac, Terachem, Turbomole, XTB and pyscf as QM engines. 

We are also using the ANI and TorchANI calculators for Artificial Neural Network (ANN) engines.

These files should be located under plugin/ase/ directory of your qr-core code.

The wrappers of those QM and ANN calculators were adapted from ASE.
