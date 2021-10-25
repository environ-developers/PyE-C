import sys

import numpy as np
import qepy
from edftpy.engine.engine_environ import EngineEnviron

try:
    from mpi4py import MPI
    mpi_comm = MPI.COMM_WORLD
    comm = mpi_comm.py2f()
except Exception:
    comm = None


def printt(s):
    if VERBOSE and ionode:
        print(s)
        sys.stdout.flush()


VERBOSE = False
ionode = mpi_comm.rank == 0

printt(f'comm:{comm}')

# INPUT FILES
pw_file = 'neutral.in'
env_file = 'environ.in'

use_environ = True

qepy.qepy_pwscf(pw_file, comm)

embed = qepy.qepy_common.embed_base()

if use_environ:

    # QEPY PARAMETERS NEEDED IN ENVIRON
    nat = qepy.ions_base.get_nat()  # number of atoms
    ntyp = qepy.ions_base.get_nsp()  # number of species
    nelec = qepy.klist.get_nelec()  # number of electrons
    atom_label = qepy.ions_base.get_array_atm().view('c')  # atom labels
    atom_label = atom_label[:, :ntyp]  # discard placeholders

    alat = qepy.cell_base.get_alat()  # lattice parameter
    at = qepy.cell_base.get_array_at()  # 3 x 3 lattice in atomic units
    gcutm = qepy.gvect.get_gcutm()  # G-vector cutoff

    ityp = qepy.ions_base.get_array_ityp()  # species indices
    zv = qepy.ions_base.get_array_zv()  # ionic charges
    tau = qepy.ions_base.get_array_tau()  # ion positions

    # PRINT PARAMETERS
    printt(f'nat={nat}')
    printt(f'nelec={nelec}')
    printt(f'ntyp={ntyp}')
    printt(f'atom_label={atom_label}')
    printt(f'alat={alat}')
    printt(f'at={at}')
    printt(f'gcutm={gcutm}')
    printt(f'ityp={ityp}')
    printt(f'zv={zv}')
    printt(f'tau={tau}')
    printt('')

    inputs = {
        'nat': nat,
        'ntyp': ntyp,
        'nelec': nelec,
        'atom_label': atom_label,
        'ityp': ityp,
        'alat': alat,
        'at': at,
        'gcutm': gcutm,
        'zv': zv,
        'tau': tau,
        'use_internal_pbc_corr': False,
    }

    # ENVIRON INIT
    printt("initializing Environ")
    environ = EngineEnviron()
    environ.initial(comm=mpi_comm, **inputs)

    # GET INITIAL DENSITY
    printt("retrieving density")
    nnt = environ.get_nnt()  # total number of grid points in Environ
    rho = np.zeros((nnt, 1), order='F')  # density register used in each step
    qepy.qepy_mod.qepy_get_rho(rho, True)
    printt(f'rho={np.sum(rho)}')

    # PRE-SCF UPDATES
    environ.pre_scf(rho)

# SCF CYCLE
printt("starting scf cycle")
embed.iterative = True
maxsteps = 60

for i in range(maxsteps):

    # QE SCF
    if i > 0: embed.initial = False

    printt("running qepy scf")
    embed.mix_coef = -1.0
    qepy.qepy_electrons_scf(0, 0, embed)
    embed.mix_coef = 0.7
    qepy.qepy_electrons_scf(2, 0, embed)

    if use_environ:

        # UPDATE ENERGY
        environ_energy = environ.calc_energy()
        printt(f'environ energy = {environ_energy}')
        printt("passing environment energy contribution to qepy")
        embed.extene = environ_energy
        printt(f"corrected energy = {embed.etotal}")

        # ENVIRON OUTPUT
        environ.print_energies()  #TODO figure out one-electron term

    conv_elec = qepy.control_flags.get_conv_elec()

    if use_environ:

        # QEPY -> ENVIRON
        printt("retrieving density")
        rho[:] = 0.0
        qepy.qepy_mod.qepy_get_rho(rho, True)
        printt(f'rho = {np.sum(rho)}')

        # ENVIRON SCF
        printt("running Environ scf")
        threshold = environ.get_threshold()
        update = not conv_elec or embed.dnorm < threshold
        environ.scf(rho, update)

        # ENVIRON -> QEPY
        printt("passing environment potential to qepy")
        qepy.qepy_mod.qepy_set_extpot(embed, environ.get_potential())

    if conv_elec: break  # TODO can this move up?

if use_environ:
    printt("potential shift")
    environ.print_potential_shift()

# FORCES
if use_environ:
    printt("calculating environment forces")
    environ_force = environ.get_force()
    printt(environ_force)

printt("calculating qepy forces")
qepy.qepy_forces(embed=embed)
forces = qepy.force_mod.get_array_force().T
printt(forces)

# SUM UP FORCES
forces += environ_force

# NORMALIZE FORCES TO REMOVE FORCE ON CELL
# TODO remove the average of the forces from the force

printt(forces)

# CLEAN ALLOCATIONS AND EXIT
if use_environ: environ.stop_run()

qepy.punch('all')
qepy.qepy_stop_run(0, what='no')
