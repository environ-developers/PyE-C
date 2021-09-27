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
MODULE control_interface
    !------------------------------------------------------------------------------------
    !
    USE env_base_scatter, ONLY: env_scatter_grid
    USE env_types_fft, ONLY: env_fft_type_descriptor
    !
    USE environ_param, ONLY: DP
    !
    USE env_global_objects, ONLY: env, setup
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: update_cell, update_ions, update_electrons, update_response, &
              add_mbx_charges
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   UPDATE METHODS
    !
    !------------------------------------------------------------------------------------
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
        ALLOCATE (at_scaled(3, 3))
        at_scaled = at * alat
        !
        CALL setup%update_cell(at_scaled)
        !
        CALL env%update_cell_dependent_quantities()
        !
        CALL setup%end_cell_update()
        !
        DEALLOCATE (at_scaled)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_ions(nat, tau, alat)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        REAL(DP), INTENT(IN) :: tau(3, nat), alat
        !
        REAL(DP), ALLOCATABLE :: tau_scaled(:, :)
        !
        !--------------------------------------------------------------------------------
        !
        ALLOCATE (tau_scaled(3, nat))
        tau_scaled = tau * alat
        !
        CALL env%update_ions(nat, tau_scaled)
        !
        DEALLOCATE (tau_scaled)
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
        REAL(DP), INTENT(IN) :: rho(setup%system_cell%dfft%nnt)
        LOGICAL, INTENT(IN), OPTIONAL :: lscatter
        !
        REAL(DP) :: aux(setup%system_cell%dfft%nnr)
        REAL(DP) :: nelec
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__MPI)
        IF (PRESENT(lscatter)) THEN
            !
            IF (lscatter) THEN
                CALL env_scatter_grid(setup%system_cell%dfft, rho, aux)
            ELSE
                aux = rho
            END IF
            !
        ELSE
            aux = rho
        END IF
        !
#else
        aux = rho
#endif
        nelec = REAL(env%system_electrons%number, DP)
        !
        CALL env%update_electrons(setup%system_cell%dfft%nnr, aux, nelec)
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
        REAL(DP), INTENT(IN) :: drho(setup%system_cell%dfft%nnr) ! # TODO may need scatter
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%update_response(setup%system_cell%dfft%nnr, drho)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_response
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE add_mbx_charges(rho, lscatter)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: rho(:)
        LOGICAL, INTENT(IN), OPTIONAL :: lscatter
        !
        REAL(DP) :: aux(setup%system_cell%dfft%nnr)
        !
        CHARACTER(LEN=80) :: local_label = 'mbx_charges'
        !
        !--------------------------------------------------------------------------------
        !
#if defined(__MPI)
        IF (PRESENT(lscatter)) THEN
            !
            IF (lscatter) THEN
                CALL env_scatter_grid(setup%system_cell%dfft, rho, aux)
            ELSE
                aux = rho
            END IF
            !
        ELSE
            aux = rho
        END IF
        !
#else
        aux = rho
        !
#endif
        CALL env%add_charges(setup%system_cell%dfft%nnr, aux, local_label)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE add_mbx_charges
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE control_interface
!----------------------------------------------------------------------------------------
