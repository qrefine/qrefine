# Quantum Refinement Module

[![Build Status](https://travis-ci.org/qrefine/qrefine.svg?branch=master)](https://travis-ci.org/qrefine/qrefine)

### Quickstart

cd into the modules subdirectory of Phenix:

```
 source ../phenix_env.sh
 
 sudo git clone --recursive https://github.com/qrefine/qrefine.git
 
 sudo libtbx.configure qrefine
 ```
 
 or cd into the modules subdirectory of CCTBX (plus a few dependencies)
 
```
 source ../build/setpaths.sh
 
 sudo git clone --recursive https://github.com/qrefine/qrefine.git
 
 sudo libtbx.configure qrefine
 ```
 
 
 ### Run Tests 

``` 
 qrefine.test
 
```
If any of the tests fail, please raise and issue here:

### Run Example 

If tests run sucessfully, then try and run an example: 

```
 qrefine.run 1uso 
```

### Help 

If you run into any trouble please ask for help:
```
 qrefine.help
```

### Contact us 

The best way to get a hold of us is via the  [issue tracker](https://github.com/qrefine/qr-core/issues)
