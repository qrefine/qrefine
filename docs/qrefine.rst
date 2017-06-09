Quantum|Refinement: final-stage refinement with restraints from Quantum Chemistry.

Authors
 -Min Zheng, Pavel Afonine, Mark Waller, Nigel Moriarty

Purpose
qr is a command line tool for refining bio-macromolecules using restraints from Quantum Chemistry.

Usage
qr is a new open-source module that carries out refinement of bio-macromolecules.
To maintain a small and agile code-base, qr is built on top of cctbx and Terachem.
The cctbx library provides most of the routines needed for x-ray refinement.
The key feature of the qr code is that it interfaces to Terachem to obtain
chemical restraints using ab initio methods.

In principle, qr only needs a data file (e.g. mtz) and a model (e.g. pdb).

qr.refine input.pdb input.mtz

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
      When you specify nproc=n, you can run the subprocesses as jobs in background with sh (default)
        If you have a multi-processor machine, use sh.
        If nproc is greater than 1 and you use run_command='sh '(or similar, sh is default) then normally you will use background=True so that all the jobs run simultaneously.


    citations =
        prints the citations of the relevant papers

  fragment


  restraint


  finalise


  test
     unit =
        all unit tests will run

     regression =
        all regression tests will run

     pdb =
        You can run tests on the entire PDB, but it takes a very long time
