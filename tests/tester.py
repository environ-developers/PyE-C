import sys
import numpy as np
from typing import Any
from configparser import ConfigParser
import pyec.environ_interface as environ

try:
    from mpi4py import MPI
    mpi_comm = MPI.COMM_WORLD
    comm = mpi_comm.py2f()
    ionode = 0
    is_ionode = mpi_comm.rank == ionode
except Exception:
    comm = None
    ionode = 0
    is_ionode = True


def pprint(item: Any = '', end: str = '\n') -> str:
    """Parallel-safe printer."""
    if is_ionode:
        print(item, end=end)
        sys.stdout.flush()


def format_labels(n: int, lables: np.ndarray) -> np.ndarray:
    """Converts atom labels to expected Fortran format."""
    flabels = np.array([[[ord(s)] for s in atom] for atom in lables])
    return flabels.view('c').T[0]


def get_system_input(filename: str) -> None:
    """docstring"""
    my_input = ConfigParser()
    my_input.read('system.ini')


nelec = 0
nat = 1
ntyp = 1
ityp = (1, )
zv = (-1.0, )
atom_labels = format_labels(ntyp, ('H  ', ))

gcutrho = 300
gcutm = gcutrho / np.pi**2

mt_corr = True

at = np.zeros((3, 3), order='F')
tau = np.zeros((3, nat), order='F')
alat = 15.0

for i in range(3):
    at[i, i] = 1.0

at *= alat

for i in range(nat):
    for j in range(3):
        tau[j, i] = 0.5

tau *= alat

environ.init_io(is_ionode, ionode, comm, 6, False)
environ.read_input()

pprint("Initializing Environ...")
environ.init_environ(comm, nelec, nat, ntyp, atom_labels, ityp, zv, mt_corr,
                     at, gcutm)

pprint("\nUpdating ions...")
environ.update_ions(nat, tau)

pprint("\nUpdating cell...")
environ.update_cell(at)

nnt = environ.get_nnt()
dvtot = np.zeros(nnt, order='F')
forces = np.zeros((3, nat), order='F')

pprint("\nCalculating potential...")
environ.calc_potential(True, dvtot, lgather=True)
pprint(f"dvtot = {np.sum(dvtot):.3e} Ry")

pprint("\nCalculating energy...")
energy = environ.calc_energy()
pprint(f"energy = {energy:.5f} Ry")

pprint("\nCalculating forces... ")
environ.calc_force(forces)
for i, force in enumerate(forces.T, 1):
    pprint("atom {}: {:.8f} {:.8f} {:.8f} Ry/bohr".format(i, *force))
