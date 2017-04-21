Quantum|Refinement: final-stage refinement with restraints from Quantum Chemistry.

Authors(s)
 -Min Zheng, Pavel Afonine, Mark Waller, Nigel Moriarty

Purpose
qr is a command line tool for refining bio-macromolecules using restraints from Quantum Chemistry. 

Usage
How qr works:
qr is a new open-source module that carries out refinement of bio-macromolecules. 
To maintain a small and agile code-base, qr is built on top of cctbx and Terachem.
The cctbx library provides most of the routines needed for x-ray refinement.
The key feature of the qr code is that it interfaces to Terachem to obtain
chemical restraints using ab initio methods.
 
In principle, qr only needs a data file (e.g. mtz) and a model(e.g. pdb).

qr.start input.pdb input.mtz 

Sensible default options are selected.

Literature:
https://journals.iucr.org/d/issues/2017/01/00/lp5021/lp5021.pdf

List of all available keywords

Options and keywords in qrefine:
       - qm_calculator
       - macro_cycles
       - micro_cycles
       - max_bond_rmsd
       - refine_sites
       - refine_adp
       - cluster_qm
       - charge_embedding
       - cluster

  refine

	base_path = None
       base path for files (default is current working directory)

    temp_dir = None
      temporary directory (it must exist)

    clean_up = None 
       At the end of the entire run the TEMP directories will be removed if clean_up is True.
       Files listed in keep_files will not be deleted.

    jobid =
       select which job to pause

    jobid =
       select which job to kill

    jobs =
       prints a list of jobs currently in memory (started,paused, and stopped)

    queue_commands =  None
	  build a command for your queueing system.
	  appended to the pbs submit file.
	  For exmaple:
	   queue_commands='#PBS -N qr'
	   queue_commands='#PBS -j oe'
	   queue_commands='#PBS -l walltime=03:00:00'
	   queue_commands='#PBS -l nodes=1:ppn=16'
         NOTE: If you set run_command=qsub (or otherwise submit to a batch queue),
         then you should set background=False, so that the batch queue can keep track of your runs.
         There is no need to use background=True in this case because all the runs go as controlled by your batch system.


	run_command = "sh "
      When you specify nproc=nn, you can run the subprocesses as jobs in background with sh (default)
        If you have a multi-processor machine, use sh.
        If nproc is greater than 1 and you use run_command='sh '(or similar, sh is default) then normally you will use background=True so that all the jobs run simultaneously.


    citations = 
        prints the citations of the relevant papers

  help
    keywords = 
        
    commands = 

  example
    refine

    cluster

  run_tests
     unit = 
        all unit tests will run

     regression = 
        all regression tests will run   

     pdb =
        You can run tests on the entire PDB, but it takes a very long time

"""qr.py is the entry point to Q|R that takes user inputs, and then constructs all of the objects
   needed to carry out the quantum refinement. “””

“””driver.py is used to drive the macro/micro-cycles for either refinement, or optimization.
   The optimization option is used only for comparison/validation, and is not intended to
   be useful as a standalone optimizer. This class also takes care of the convergence criteria.
   It delegates the minimization procedure itself to the L-BFGS implementation in CCTBX.
   The driver.py requires the target and gradient (energy and force from QM) in order to minimize.”””

“””calculator.py  handles the weight factors, scaling, and is used to convert input parameters
   to which can be used for either quantum or traditional refinement.
   An adaptive restraints weight factor calculator is implemented, whereby the weight factor is
   doubled if a sufficiently large bond-RMSD is observed. Conversely, if a sufficiently small
   bond-RMSD is observed, then the weight factor is halved.””

“””restraints.py contains two classes for either quantum refinement, or for standard refinement.
   The calculation of the restraints are delegated to either ASE
   for quantum-based, or CCTBX for standard refinement.”””

“””results.py  stores and handles all of the data needed for logging the results of the refinement”””

We then need to create a set of objects to carry out the computation:

- fmodel (crystallographic information)

- calculator (composite object)
  - restraints_manager (computes energy and gradients using either qm codes or cctbx (standard))
  - geometery_restraints_manager (analyses geometry e.g. bond RMSDs)
  - weights (scale factors needed to scale up or down data versus restraints contributions)

Then we process them by the refinement/optimization engine, driver.py:

```
 driver.refine(params     = params,
               fmodel     = fmodel,
               calculator = calculator_manager,
               results    = results)

```

- results_manager (store all reportable infomation, and write it out as a log, and also write our final pdb structure.)

```
 results.finalize()
 ```
