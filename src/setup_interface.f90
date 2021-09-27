!----------------------------------------------------------------------------------------
!
! Copyright (C) 2021 ENVIRON (www.quantum-environ.org)
!
!----------------------------------------------------------------------------------------
!
!     This file is part of Environ version 2.0
!
!     Environ 2.0 is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 2 of the License, or
!     (at your option) any later version.
!
!     Environ 2.0 is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more detail, either the file
!     `License' in the root directory of the present distribution, or
!     online at <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------------------
!
! Authors: Edan Bainglass   (Department of Physics, UNT)
!          Matthew Truscott (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE setup_interface
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE environ_param, ONLY: DP
    !
    USE environ_input, ONLY: read_environ_input
    USE env_base_input
    !
    USE env_global_objects, ONLY: env, setup
    USE class_destructor, ONLY: clean
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: init_io, read_input, init_environ, clean_environ
    !
    !------------------------------------------------------------------------------------

    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_io(ionode, ionode_id, comm, program_unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: ionode
        INTEGER, INTENT(IN) :: ionode_id
        INTEGER, INTENT(IN) :: comm
        INTEGER, INTENT(IN) :: program_unit
        !
        !--------------------------------------------------------------------------------
        !
        CALL io%init(ionode, ionode_id, comm, program_unit, ionode)
        !
        io%debug_unit = io%find_free_unit()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_io
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE read_input(filename)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80), INTENT(IN), OPTIONAL :: filename
        !
        CHARACTER(LEN=80) :: sub_name = 'read_input'
        !
        !--------------------------------------------------------------------------------
        !
        CALL read_environ_input(filename)
        !
        io%verbosity = verbose
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE read_input
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ(comm, nelec, nat, ntyp, atom_label, ityp, zv, &
                            use_internal_pbc_corr, alat, at, gcutm)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm
        INTEGER, INTENT(IN) :: nelec, nat, ntyp
        CHARACTER(LEN=3), INTENT(IN) :: atom_label(ntyp)
        INTEGER, INTENT(IN) :: ityp(nat)
        REAL(DP), INTENT(IN) :: zv(ntyp)
        LOGICAL, INTENT(IN) :: use_internal_pbc_corr
        REAL(DP), INTENT(IN) :: alat
        REAL(DP), INTENT(IN) :: at(3, 3)
        REAL(DP), INTENT(IN) :: gcutm
        !
        REAL(DP), ALLOCATABLE :: at_scaled(:, :)
        REAL(DP) :: gcutm_scaled
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ'
        !
        !--------------------------------------------------------------------------------
        !
        CALL setup%init() ! flags and solvers
        !
        CALL env_allocate_mp_buffers()
        !
        ALLOCATE (at_scaled(3, 3))
        at_scaled = at * alat
        gcutm_scaled = gcutm / alat**2
        !
        CALL setup%init_cell(gcutm_scaled, comm, at_scaled)
        !
        DEALLOCATE (at_scaled)
        !
        CALL setup%init_cores(gcutm_scaled, use_internal_pbc_corr)
        !
        CALL env%init(setup, nelec, nat, ntyp, atom_label, ityp, zv)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE clean_environ()
        !--------------------------------------------------------------------------------
        !
        CALL clean%all()
        !
        CALL env_deallocate_mp_buffers()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE clean_environ
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE setup_interface
!----------------------------------------------------------------------------------------
