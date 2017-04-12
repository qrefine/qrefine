# Quantum Refinement Module

[![Build Status](https://travis-ci.org/qrefine/qrefine.svg?branch=master)](https://travis-ci.org/qrefine/qrefine)

### Quickstart

Please first install PHENIX, see https://www.phenix-online.org/
 
Once you have PHENIX installed, go to the directory where you installed phenix.

```
 source phenix_env.sh
 cd modules 
 git clone https://github.com/qrefine/qrefine.git
 cd qrefine
 ./patch.sh
 ```
 Note: you may need to use sudo depending on the permissions of your PHENIX installation.
 
 ### Run Tests 

``` 
 qr.run_tests
 
```
If any of the tests fail, please raise and issue here:

### Run Example 

If tests run successfully, then try and run an example: 

```
 qr.example refine 
``` 

 or 
 
```
 qr.example cluster 

```

### Help 

If you run into any trouble please ask for help:
```
 qr.help
```

### Contact us 

The best way to get a hold of us is via the  [issue tracker](https://github.com/qrefine/qr-core/issues)
