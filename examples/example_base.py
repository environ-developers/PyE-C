import sys

import numpy as np
import qepy
import pyec.setup_interface as environ_setup
import pyec.control_interface as environ_control
import pyec.calc_interface as environ_calc
import pyec.output_interface as environ_output

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

# QEPY SETUP
qepy.qepy_pwscf(pw_file, comm)

embed = qepy.qepy_common.embed_base()

if use_environ:

    # QEPY PARAMETERS NEEDED IN ENVIRON
    nat = qepy.ions_base.get_nat()  # number of atoms
    ntyp = qepy.ions_base.get_nsp()  # number of species
    nelec = qepy.klist.get_nelec()  # number of electrons
    atom_label = qepy.ions_base.get_array_atm().view('c')  # atom labels
    atom_label = atom_label[:, :ntyp]  # discard placeholders*

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

    # ENVIRON INIT
    printt('setting up io')
    environ_setup.init_io(ionode, 0, comm, 6)

    printt("reading Environ input")
    environ_setup.read_input(env_file)

    printt("initializing Environ")
    environ_setup.init_environ(comm, nelec, nat, ntyp, atom_label, ityp, zv,
                               False, alat, at, gcutm)

    # PRE-SCF UPDATES
    printt("updating ions")
    environ_control.update_ions(nat, tau, alat)
    printt("updating cell")
    environ_control.update_cell(at, alat)

    # CELL-DEPENDENT PARAMATERS
    nnt = environ_control.get_nnt()  # total number of grid points in Environ
    rho = np.zeros((nnt, 1), order='F')  # density register used in each step
    dvtot = np.zeros(nnt, order='F')  # potential contribution per scf step

    # GET INITIAL DENSITY
    printt("retrieving density")
    qepy.qepy_mod.qepy_get_rho(rho, True)
    printt(f'rho = {np.sum(rho)}')

    # UPDATE ELECTRONS
    printt("updating electrons")
    environ_control.update_electrons(rho, True)

    # CALCULATE ENVIRONMENT CONTRIBUTION TO THE POTENTIAL
    printt("calculating environment contribution to the potential")
    restart = environ_control.is_restart()
    environ_calc.calc_potential(restart, dvtot, lgather=True)

# SCF CYCLE
printt("starting scf cycle")
embed.iterative = True
maxsteps = 60

for i in range(maxsteps):

    # QE SCF
    if i > 0: embed.initial = False

    printt("running qepy scf")
    embed.mix_coef = -1.0
    qepy.qepy_electrons_scf(2, 0, embed)
    embed.mix_coef = 0.7
    qepy.qepy_electrons_scf(2, 0, embed)

    if use_environ:

        # UPDATE ENERGY
        environ_energy = environ_calc.calc_energy()
        printt(f'environ energy = {environ_energy}')
        printt("passing environment energy contribution to qepy")
        embed.extene = environ_energy
        printt(f"corrected energy = {embed.etotal}")

        # ENVIRON OUTPUT
        environ_output.print_energies()  #TODO figure out one-electron term

    conv_elec = qepy.control_flags.get_conv_elec()

    if use_environ:

        # QEPY -> ENVIRON
        printt("retrieving density")
        rho[:] = 0.0
        qepy.qepy_mod.qepy_get_rho(rho, True)
        printt(f'rho = {np.sum(rho)}')

        # ENVIRON SCF
        printt("running Environ scf")
        printt("updating electrons")
        environ_control.update_electrons(rho, True)

        # check if Environ should compute its potential contribution
        # (either the threshold has been met, or this is a restarted job)
        if embed.dnorm > 0.0:
            threshold = environ_control.get_threshold()
            update = not conv_elec or embed.dnorm < threshold

        printt("calculating environment contribution to the potential")
        environ_calc.calc_potential(update, dvtot, lgather=True)
        printt(f'dvtot = {np.sum(dvtot)}')

        # ENVIRON -> QEPY
        printt("passing environment potential to qepy")
        qepy.qepy_mod.qepy_set_extpot(embed, dvtot)

    if conv_elec: break  # TODO can this move up?

if use_environ:
    printt("potential shift")
    environ_output.print_potential_shift()

# FORCES
if use_environ:
    printt("calculating environment forces")
    environ_force = np.zeros((3, nat), dtype=float, order='F')
    environ_calc.calc_force(environ_force)
    printt(environ_force.T)
    printt("passing environment forces to qepy")
    qepy.qepy_mod.qepy_set_extforces(embed, environ_force)

printt("calculating qepy forces")
qepy.qepy_forces(embed=embed)
forces = qepy.force_mod.get_array_force().T

printt(forces)

# CLEAN ALLOCATIONS AND EXIT
if use_environ: environ_setup.clean_environ()

qepy.punch('all')
qepy.qepy_stop_run(0, what='no')
