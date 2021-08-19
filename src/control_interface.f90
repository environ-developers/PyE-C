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
    USE env_base_scatter, ONLY: env_scatter_grid
    !
    USE environ_param, ONLY: DP
    !
    USE class_environ, ONLY: env
    !
    USE env_types_fft, ONLY: env_fft_type_descriptor
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: update_potential, update_cell, update_ions, update_electrons, update_response
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
    SUBROUTINE update_potential(vltot, lscatter)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: vltot(env%system_cell%nnr)
        LOGICAL, INTENT(IN), OPTIONAL :: lscatter
        !
        REAL(DP) :: aux(env%system_cell%nnr)
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__MPI)
        IF (PRESENT(lscatter)) THEN
            IF (lscatter) THEN
                CALL env_scatter_grid(env%system_cell%dfft, vltot, aux)
            ELSE
                aux = vltot
            END IF
        ELSE
            aux = vltot
        END IF
#else
        aux = vltot
#endif
        !
        CALL env%init_potential(env%system_cell%nnr, aux)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_potential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_cell(at, alat)
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
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_cell
    !------------------------------------------------------------------------------------
    !>
    !! // TODO consider whether ions are ever changed or just need updated positions
    !------------------------------------------------------------------------------------
    SUBROUTINE update_ions(nat, ntyp, ityp, zv, tau, alat)
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
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_ions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_electrons(rho, lscatter)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: rho(env%system_cell%dfft%nnt)
        LOGICAL, INTENT(IN), OPTIONAL :: lscatter
        !
        REAL(DP) :: aux(env%system_cell%dfft%nnr)
        REAL(DP) :: nelec
        !
        !--------------------------------------------------------------------------------
        !
        !
#if defined(__MPI)
        IF (PRESENT(lscatter)) THEN
            IF (lscatter) THEN
                CALL env_scatter_grid(env%system_cell%dfft, rho, aux)
            ELSE
                aux = rho
            END IF
        ELSE
            aux = rho
        END IF
#else
        aux = rho
#endif
        !
        nelec = REAL(env%system_electrons%number, DP)
        !
        CALL env%init_electrons(env%system_cell%dfft%nnr, aux, nelec)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_electrons
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_response(drho)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        ! // TODO may need scatter
        REAL(DP), INTENT(IN) :: drho(env%system_cell%dfft%nnr)
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%init_response(env%system_cell%dfft%nnr, drho)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_response
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
