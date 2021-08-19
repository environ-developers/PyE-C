import sys

import numpy as np
import qepy
import pyec.setup_interface as environ_setup
import pyec.control_interface as environ_control
import pyec.calc_interface as environ_calc
import pyec.output_interface as environ_output

from mpi4py import MPI

VERBOSE=True
def printt(s):
    if VERBOSE:
        print(s)
        sys.stdout.flush()

comm = MPI.COMM_WORLD
comm = comm.py2f()
printt(comm)

fname = 'dielectric.in'

qepy.qepy_pwscf(fname, comm)

embed = qepy.qepy_common.embed_base()
embed.exttype = 0
embed.finish = False

## CONTAINER FOR ENVIRON THINGS
nat = qepy.ions_base.get_nat()
ntyp = qepy.ions_base.get_nsp()
nelec = qepy.klist.get_nelec()
atom_label = qepy.ions_base.get_array_atm()

alat = qepy.cell_base.get_alat()
at = qepy.cell_base.get_array_at()
me = 0
root = 0
gcutm = qepy.gvect.get_gcutm()
e2_in = qepy.constants.e2

ityp = qepy.ions_base.get_array_ityp()
zv = qepy.ions_base.get_array_zv()
tau = qepy.ions_base.get_array_tau()

vltot = qepy.scf.get_array_vltot()
nnr = MPI.COMM_WORLD.allreduce(vltot.size, op=MPI.SUM)

printt(f'nnr={nnr}')

rho = np.zeros((nnr, 1), order='F')
rhohist = np.zeros((nnr, 1), order='F')
dvtot = np.zeros((nnr), order='F')
printt("qepy get rho")
qepy.qepy_mod.qepy_get_rho(rho, True)

printt(f'nat={nat}')
printt(f'nelec={nelec}')
printt(f'ntyp={ntyp}')
printt(f'atom_label={atom_label[:, :ntyp]}')
printt(f'alat={alat}')
printt(f'at={at}')
printt(f'gcutm={gcutm}')
printt(f'ityp={ityp}')
printt(f'zv={zv}')
printt(f'tau={tau}')
printt(f'vltot={vltot.shape}')
printt(f'rho={np.sum(rho)}')

# ENVIRON INIT
printt('io')
environ_setup.init_io('PW', True, 0, comm, 6)
printt("base 1")
environ_setup.init_base_first(nelec, nat, ntyp, atom_label[:, :ntyp], False)
printt("base 2")
environ_setup.init_base_second(alat, at, comm, me, root, gcutm, e2_in)

# update functions
printt("ions")
environ_control.update_ions(nat, ntyp, ityp, zv[:ntyp], tau, alat)
printt("cell")
environ_control.update_cell(at, alat)
printt("potential")
environ_control.update_potential(vltot)
printt("electrons")
environ_control.update_electrons(rho, True)

# calculator interface
printt("calcpotential")
environ_calc.calc_potential(False, dvtot, lgather=True)

printt("scf")
nstep = 2
for i in range(nstep):
    # QE SCF
    if i == 0 :
        initial = True
    else :
        initial = False
    embed.initial = initial
    embed.mix_coef = -1.0
    qepy.qepy_electrons_scf(0, 0, embed)
    embed.mix_coef = 0.7
    qepy.qepy_electrons_scf(2, 0, embed)

    # QEPY -> ENVIRON
    qepy.qepy_calc_energies(embed)
    rho[:] = 0.0
    printt("rho gather")
    qepy.qepy_mod.qepy_get_rho(rho, True)

    # ENVIRON SCF
    printt("rho scatter")
    environ_control.update_electrons(rho, True)
    printt("v gather")
    environ_calc.calc_potential(True, dvtot, lgather=True)

    # UPDATE ENERGY
    environ_energy = environ_calc.calc_energy()
    embed.etotal += environ_energy

    # ENVIRON OUTPUT
    environ_output.print_potential_shift()
    environ_output.print_energies()

    # ENVIRON -> QEPY
    printt("v scatter")
    qepy.qepy_mod.qepy_set_extpot(embed, dvtot, True)

    printt(f"corrected energy = {embed.etotal}")


qepy.qepy_calc_energies(embed)
etotal = embed.etotal

qepy.qepy_forces(0)
forces = qepy.force_mod.get_array_force().T

force_environ = np.zeros((3, nat), dtype=float, order='F')
environ_calc.calc_force(force_environ)
printt(force_environ.T)

# TODO clean up environ stuff

qepy.punch('all')
qepy.qepy_stop_run(0, what = 'no') 

