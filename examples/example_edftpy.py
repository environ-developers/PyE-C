import sys

import numpy as np
import qepy
from edftpy.engine.engine_environ import EngineEnviron

from mpi4py import MPI

VERBOSE=True
def printt(s):
    if VERBOSE:
        print(s)
        sys.stdout.flush()

comm = MPI.COMM_WORLD
comm = comm.py2f()
printt(f'comm:{comm}')

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
nnr = vltot.size

rho = np.zeros((nnr, 1), order='F')
rhohist = np.zeros((nnr, 1), order='F')
dvtot = np.zeros((nnr), order='F')
qepy.qepy_mod.qepy_get_rho(rho, False)

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
printt(f'nnr={nnr}')
printt(f'vltot={vltot.shape}')
printt(f'rho={np.sum(rho)}')

inputs = {
        'nat': nat,
        'ntyp': ntyp,
        'nelec': nelec,
        'atom_label': atom_label,
        'ityp': ityp,
        'alat': alat,
        'at': at,
        'gcutm': gcutm,
        'e2': e2_in,
        'zv': zv,
        'tau': tau,
        'vltot': vltot,
        'nnr': nnr,
        'rho': rho,
}

# ENVIRON INIT
printt('init')
environ = EngineEnviron()
environ.initial(comm, **inputs)

nstep = 3
for i in range(nstep):
    # QE SCF
    printt('scf')
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
    qepy.qepy_mod.qepy_get_rho(rho, False)

    # ENVIRON SCF
    environ.scf(rho)

    # UPDATE ENERGY
    environ_energy = environ.get_energy()
    embed.etotal += environ_energy * 2.0

    # ENVIRON -> QEPY
    qepy.qepy_mod.qepy_set_extpot(embed, environ.get_potential())

    printt(f"corrected energy = {embed.etotal}")


qepy.qepy_calc_energies(embed)
etotal = embed.etotal

qepy.qepy_forces(0)
forces = qepy.force_mod.get_array_force().T

forces += environ.get_force() * 2.0
printt(forces)

environ.clean()

qepy.punch('all')
qepy.qepy_stop_run(0, what = 'no')
