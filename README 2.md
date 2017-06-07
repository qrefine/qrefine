# Quantum Refinement Project                                                                                                                                                                                           

Quantum refinement uses quantum chemistry to better refine bio-macromolecules.

# Quick start

### Download qr and add the path of qr-core to the PYTHONPATH
```
git clone https://github.com/qrefine/qr-core.git
```
Add the absolute path of qr-core to the PYTHONPATH
```
export PYTHONPATH=$PYTHONPATH:absolute_path_of_qr-core
```

### Download and install nightly build Phenix
http://phenix-online.org/download/nightly_builds.cgi

after Phenix installation, set the environment value of PHENIX_TRUST_OTHER_ENV
```
export PHENIX_TRUST_OTHER_ENV=1
```

### Install ase, jpype and yoink
To not worry about the PYTHONPATH set up, install ase and jpype using phenix.python.
```
phenix.python -m  pip install -U pip
phenix.python -m  pip install -r requirements.txt 
```
Otherwise add the installtion paths of ase and jpype to PYTHONPATH
```
export PYTHONPATH=$PYTHONPATH:absolute_path_of_ase
export PYTHONPATH=$PYTHONPATH:absolute_path_of_jpype
```
yoink is available here https://github.com/qrefine/qr-plugin-yoink. After get it from GitHub, move it to the folder plugin under qr-core
```
git clone https://github.com/qrefine/qr-plugin-yoink.git && mv ./qr-plugin-yoink/* absolute_path_of_qr-core/plugin/
```

QM engine in ASE*

Now Q|R can use Mopac, Terachem, Turbomole and pyscf as QM engine via ASE.

Those four QM engines are called by four files mopac_qr.py, terachem_qr.py, turbomole_qr.py and pyscf_qr.py.

They are not in ASE, available https://github.com/qrefine/qr-plugin-ase. To get those files:
```
git clone https://github.com/qrefine/qr-plugin-ase.git
```
Then get the directory of ase.calculators 

       >>> import os.path,ase.calculators

       >>> os.path.dirname(ase.calculators.__file__)
and copy those files in this directory

Tips for jpype*

jpype could not set up proper environment for JAVA always. Better to set up JAVA environment manually. jdk version > 1.8.
```
export JAVA_HOME=absolute_path_of_java_home
export JAVA_LIB_PATH=absolute_path_of_java_lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$JAVA_LIB_PATH/server

```


### Run qr tests

set up the environmental value QR_REPO_PARENT as the parent folder of qr_core
```
export QR_REPO_PARENT=parent_folder_of_qr-core
cd qr-core/tests
cctbx.python run_tests.py
```
Please ensure all tests pass, if not please [raise an issue](https://github.com/qrefine/qr-core/issues)

### Run qr
```
cctbx.python qr.py input.pdb input.mtz restraints=qm qm_calculator=terachem/pyscf/mopac/turbomole
```

### Contact us 

The best way to get a hold of us is via the  [issue tracker](https://github.com/qrefine/qr-core/issues)

