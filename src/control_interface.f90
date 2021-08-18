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
! Authors: Edan Bainglass (Department of Physics, UNT)
! // TODO update names of these subroutines
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE control_interface
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP
    !
    USE class_environ, ONLY: env
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: init_potential, init_cell, init_ions, init_electrons, init_response
    !
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               UPDATE METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Called at every ionic step
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_potential(nnr, vltot)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: vltot(nnr)
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%init_potential(nnr, vltot)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_potential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_cell(at, alat)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: at(3, 3), alat
        !
        REAL(DP), ALLOCATABLE :: at_scaled(:, :)
        !
        !--------------------------------------------------------------------------------
        !
        ALLOCATE(at_scaled(3, 3))
        at_scaled = at * alat
        !
        CALL env%init_cell(at_scaled)
        !
        DEALLOCATE(at_scaled)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_ions(nat, ntyp, ityp, zv, tau, alat)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat, ntyp
        INTEGER, INTENT(IN) :: ityp(nat)
        REAL(DP), INTENT(IN) :: zv(ntyp), tau(3, nat), alat
        !
        REAL(DP), ALLOCATABLE :: tau_scaled(:, :)
        !
        !--------------------------------------------------------------------------------
        !
        ALLOCATE(tau_scaled(3, nat))
        tau_scaled = tau * alat
        !
        CALL env%init_ions(nat, ntyp, ityp, zv, tau_scaled)
        !
        DEALLOCATE(tau_scaled)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_ions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_electrons(nnr, rho, nelec)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: rho(nnr)
        REAL(DP), INTENT(IN), OPTIONAL :: nelec
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%init_electrons(nnr, rho, nelec)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_electrons
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_response(nnr, drho)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        REAL(DP), INTENT(IN) :: drho(nnr)
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%init_response(nnr, drho)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_response
    !------------------------------------------------------------------------------------    
    !>
    !!
    !------------------------------------------------------------------------------------
    ! SUBROUTINE add_potential(nnr, added_potential, label)
    !     !--------------------------------------------------------------------------------
    !     !
    !     IMPLICIT NONE
    !     !
    !     INTEGER, INTENT(IN) :: nnr
    !     REAL(DP), INTENT(IN) :: added_potential(nnr)
    !     CHARACTER(LEN=80), INTENT(IN), OPTIONAL :: label
    !     !
    !     CHARACTER(LEN=80) :: local_label = 'additional_potential'
    !     !
    !     CHARACTER(LEN=80) :: sub_name = 'add_potential'
    !     !
    !     !--------------------------------------------------------------------------------
    !     !
    !     IF (PRESENT(label)) local_label = label
    !     !
    !     CALL env%additional_potential%init(env%environment_cell, local_label)
    !     !
    !     env%need_additional_potential = .TRUE.
    !     !
    !     !--------------------------------------------------------------------------------
    ! END SUBROUTINE add_potential
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE control_interface
!----------------------------------------------------------------------------------------
