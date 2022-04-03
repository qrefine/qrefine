from __future__ import division
from __future__ import print_function
from tokenize import Floatnumber

# LIBTBX_SET_DISPATCHER_NAME qr.granalyse
import numpy as np
import os
import sys
from qrefine.utils.mathbox import get_grad_mad, get_grad_angle
from libtbx.utils import Usage
import iotbx.pdb


def get_help():
    raise Usage(
        """
    qr.granalyse analyses gradients obtained from 'qr.refine mode=gtest' runs.
    Will write the gradient or difference gradient into the PDB.
    The reference gradient is found automatically if not explicitly specified.

    Examples:
    i)  qr.granalyse model.pdb  
    ii) qr.granalyse model.pdb --ref 3-15.npy

    Options:
      --ref <npy files> (set reference gradient)
      --occ write to occupancy field (instead of beta) 
      --grad (write gradient instead of difference gradient into pdb)
      --help  (print this help)
    """
    )
    sys.exit(0)
    return


def id_file(filename):
    g_mode = int(filename[0])
    clustersize = int(filename[2:].split(".")[0])
    return g_mode, clustersize


def get_deviations(log, ref_grad, grad):
    import ase.units as ase_units

    """
  calculates several deviation metrics for two gradient vectors
  """
    # convert from unit less to kcal/mol
    unit_convert = ase_units.mol / ase_units.kcal  # ~ 23.06
    ref_grad = ref_grad / unit_convert
    grad = grad / unit_convert
    ref_max = max(abs(i) for i in ref_grad)
    ref_min = min(abs(i) for i in ref_grad)
    ref_gnorm = np.linalg.norm(ref_grad)
    # grad=np.array(grad)

    print(" d(angle)  %f" % (get_grad_angle(grad, ref_grad)), file=log)
    print(" d(gnorm)  %f" % (abs(np.linalg.norm(grad) - ref_gnorm)), file=log)
    print(" d(max_g)  %f" % (abs(max(abs(i) for i in grad) - ref_max)), file=log)
    print(" d(min_g)  %f" % (abs(min(abs(i) for i in grad) - ref_min)), file=log)
    print(" RMSD      %f" % (rmsd(grad, ref_grad)), file=log)
    print(" MAD       %f" % (get_grad_mad(grad, ref_grad)), file=log)
    print(" ", file=log)


def get_grad_delta(ref_grad, grad):
    return np.abs(ref_grad - grad)


# error sccore aka weighted delta gradient (final paper version)
def get_grad_wdelta(ref, g, name, do_debug):
    dim3 = int(ref.shape[0])
    dim = int(dim3 / 3)
    d = np.zeros(dim)
    norm = np.zeros(dim)
    gshape = np.reshape(g, (dim, 3))
    rshape = np.reshape(ref, (dim, 3))
    data = []
    # atom-wise norm of reference vector
    for i in range(dim):
        # index offset for cartesian components
        ii = 0 + 3 * i
        norm[i] = np.linalg.norm(ref[ii : ii + 3])

    # median of reference vector norms
    median = np.median(norm)

    # atom-wise error score (wdelta)
    for i in range(dim):
        ii = 0 + 3 * i
        i_g = g[ii : ii + 3]
        i_ref = ref[ii : ii + 3]
        dgrad = i_g - i_ref
        i_norm_delta = np.linalg.norm(dgrad)
        i_max = max(norm[i], median)
        # final weighted difference gradient for atom i
        val = 100 * (i_norm_delta / i_max)
        # limit max value to 100 for rare cases.
        d[i] = 100 if val > 100 else val
        data.append(
            [
                i,
                d[i],
                norm[i],
                i_norm_delta,
                median,
                i_max,
                g[ii],
                g[ii + 1],
                g[ii + 2],
                ref[ii],
                ref[ii + 1],
                ref[ii + 2],
            ]
        )
    if do_debug:
        np.savetxt(
            name + "-delta.txt",
            data,
            fmt="atom=%i  w.delta=%6.2f  |ref|=%5.2f  |g|=%5.2f  median(|ref|)=%5.2f max=%5.2f gradient=%f %f %f  reference=%f %f %f",
        )
        np.savetxt(name + ".txt", gshape, fmt="%12.8f")
        np.savetxt("reference_gradient.txt", rshape, fmt="%12.8f")
    return d


# apply custom weight based on gtest ID to allow easy sorting
def sorting_weight(name):
    g_mode, size = id_file(name)
    weight = g_mode * 1e4 + size
    return weight


def set_ph_field(ph, val, field):
    # sets val to respective PDB field.
    nat_hierachy = len(ph.hierarchy.atoms())
    nat_numpy = val.shape[0]
    assert nat_numpy - nat_hierachy == 0, "wrong number of atoms!"

    count = 0
    if "occ" in field.lower():
        for a in ph.hierarchy.atoms():
            a.set_occ(val[count])
            count += 1
    elif "beta" in field.lower():
        for a in ph.hierarchy.atoms():
            a.set_b(val[count])
            count += 1
    else:
        assert mode > 0, "wrong value in set_ph_field function"


def atomic_mean_deviation(x):
    nat = int(x.shape[0] / 3)
    new = np.zeros(nat)
    for i in range(nat):
        ii = 0 + 3 * i
        new[i] = (x[ii] + x[ii + 1] + x[ii + 2]) / 3
    return new


def rmsd(V, W):
    dim = V.shape[0]
    d = 0.0
    for i in range(dim):
        d += (V[i] - W[i]) ** 2.0
    return np.sqrt(d / dim)


#######################################
def run(args, log):
    field = "beta"
    ref_name = None
    do_grad = False
    do_delta = False
    do_wdelta = True
    do_expan = False
    do_debug = False
    max_val = 30

    args = sys.argv[1:]
    if ("--ref") in args:
        ref_name = args[args.index("--ref") + 1]
    if ("-ref") in args:
        ref_name = args[args.index("-ref") + 1]
    if ("--occ") in args:
        field = "occ"
    if ("--grad") in args:
        do_grad = True
    if ("--delta") in args:
        do_delta = True
    if ("--help") in args:
        get_help()
    if ("--debug") in args:
        do_debug = True
    if len(args) < 1:
        get_help()

    pdbname = args[0]
    print("Reading: %s" % (pdbname))
    pdb_inp = iotbx.pdb.input(file_name=pdbname)
    pdb_obj = iotbx.pdb.hierarchy.input(file_name=pdbname)
    ph = pdb_inp.construct_hierarchy()
    nat = len(pdb_obj.hierarchy.atoms())

    # vectors
    ref_grad = None
    # idl=[]
    iname = []
    igrad = []

    # find all npy
    files = []
    for f in os.listdir("."):
        if f.endswith(".npy"):
            files.append(f)

    # sort & label
    sfiles = sorted(files, key=sorting_weight)
    n_files = len(sfiles)
    print("Found the following files:", file=log)
    print(sfiles, file=log)
    print("", file=log)
    for i in range(0, n_files):
        iname.append(sfiles[i])
        igrad.append(np.load(iname[i]))
    # set or read reference gradient
    if ref_name is None:
        if "0-0.npy" in iname:
            print("Found supercell calc")
            ref_idx = 0
        else:
            ref_idx = n_files - 1  # last gradient should have best parameters
        ref_name = iname[ref_idx]
        ref_grad = np.array(igrad[ref_idx])
    else:
        ref_grad = np.load(ref_name)
    print("reference gradient taken from: %s \n" % (ref_name), file=log)

    # loop over available gradients and compare
    for i in range(n_files):
        print("     ~g_mode - max_res:", iname[i][:-4], file=log)
        get_deviations(log, ref_grad, igrad[i])
        if do_delta:
            dg = np.abs(ref_grad - igrad[i])
            delta = atomic_mean_deviation(dg)
        elif do_grad:
            delta = atomic_mean_deviation(igrad[i])
        elif do_wdelta:
            delta = get_grad_wdelta(ref_grad, igrad[i], iname[i][:-4], do_debug)
        set_ph_field(ph=pdb_obj, val=delta, field=field)
        output_pdb = "%s.pdb" % (iname[i][:-4])
        pdb_obj.hierarchy.write_pdb_file(file_name=output_pdb)

    return 0


if __name__ == "__main__":
    log = sys.stdout
    run(args=sys.argv[1:], log=log)
