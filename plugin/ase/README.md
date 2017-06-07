# qr-plugin-ase
A plugin for the Atomic Simulation Environment(ASE) (https://wiki.fysik.dtu.dk/ase/)

Quantum refinement code Q|R (https://github.com/qrefine/qr-core) calls QM calculator through ASE.

Now Q|R can use Mopac, Terachem, Turbomole and pyscf as QM engine. 

These files should be located under plugin/ase/ directory of your qr-core code.

The wrappers of those four QM calculators were slightly modified.
