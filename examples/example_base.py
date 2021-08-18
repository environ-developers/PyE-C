import numpy as np
import qepy
import pyec.setup_interface as environ_setup
import pyec.control_interface as environ_control
import pyec.calc_interface as environ_calc
import pyec.output_interface as environ_output

from mpi4py import MPI

comm = MPI.COMM_WORLD
comm = comm.py2f()
print(comm)

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

print(f'nat={nat}')
print(f'nelec={nelec}')
print(f'ntyp={ntyp}')
print(f'atom_label={atom_label[:, :ntyp]}')
print(f'alat={alat}')
print(f'at={at}')
print(f'gcutm={gcutm}')
print(f'ityp={ityp}')
print(f'zv={zv}')
print(f'tau={tau}')
print(f'nnr={nnr}')
print(f'vltot={vltot.shape}')
print(f'rho={np.sum(rho)}')

# ENVIRON INIT
#print('io')
environ_setup.init_io('PW', True, 0, comm, 6)
#print("base 1")
environ_setup.init_base_first(nelec, nat, ntyp, atom_label[:, :ntyp], False)
#print("base 2")
environ_setup.init_base_second(alat, at, comm, me, root, gcutm, e2_in)

# update functions (TODO: rename the interface functions)
#print("ions")
environ_control.init_ions(nat, ntyp, ityp, zv[:ntyp], tau, alat)
#print("cell")
environ_control.init_cell(at, alat)
#print("potential")
environ_control.init_potential(nnr, vltot)
#print("electrons")
environ_control.init_electrons(nnr, rho, nelec)

# calculator interface
environ_calc.calc_potential(False, nnr, dvtot)

nstep = 3
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
    qepy.qepy_mod.qepy_get_rho(rho, False)

    # ENVIRON SCF
    environ_control.init_electrons(nnr, rho, nelec)
    environ_calc.calc_potential(True, nnr, dvtot)

    # UPDATE ENERGY
    environ_energy = np.zeros((1), dtype=float)
    environ_calc.calc_energy(environ_energy)
    embed.etotal += environ_energy[0]

    # ENVIRON OUTPUT
    environ_output.print_potential_shift()
    environ_output.print_energies()

    # ENVIRON -> QEPY
    qepy.qepy_mod.qepy_set_extpot(embed, dvtot, False)

    print(f"corrected energy = {embed.etotal}")


qepy.qepy_calc_energies(embed)
etotal = embed.etotal

qepy.qepy_forces(0)
forces = qepy.force_mod.get_array_force().T

force_environ = np.zeros((3, nat), dtype=float, order='F')
environ_calc.calc_force(nat, force_environ)
print(force_environ.T)

# TODO clean up environ stuff

qepy.punch('all')
qepy.qepy_stop_run(0, what = 'no')
