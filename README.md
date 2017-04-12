# Quantum Refinement Module

[![Build Status](https://travis-ci.org/qrefine/qrefine.svg?branch=master)](https://travis-ci.org/qrefine/qrefine)

### Quickstart

Please install PHENIX first, see https://www.phenix-online.org/
 
Once you have PHENIX installed, cd into the modules subdirectory of Phenix:

``` 
 sudo git clone --recursive https://github.com/qrefine/qrefine.git
 chmod +x patch.sh
 ./patch.sh
 ```
 you may need to use sudo depending on the permissions of your PHENIX installation.
 
 ### Run Tests 

``` 
 qr.test
 
```
If any of the tests fail, please raise and issue here:

### Run Example 

If tests run successfully, then try and run an example: 

```
 qr.example refine 
 
 '''
 or 
 
 '''
  
 qr.example cluster 

```

### Help 

If you run into any trouble please ask for help:
```
 qr.help
```

### Contact us 

The best way to get a hold of us is via the  [issue tracker](https://github.com/qrefine/qr-core/issues)
