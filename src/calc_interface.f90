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
!
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE calc_interface
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: DP
    !
    USE env_base_scatter, ONLY: env_gather_grid
    !
    USE class_environ, ONLY: env
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: calc_potential, calc_energy, calc_force, calc_denergy
    !------------------------------------------------------------------------------------
    
    !
    !------------------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               COMPUTATION METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_potential(update, nnr, dvtot, local_verbose, lgather)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: update
        INTEGER, INTENT(IN) :: nnr
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        LOGICAL, INTENT(IN), OPTIONAL :: lgather
        !
        REAL(DP), INTENT(OUT) :: dvtot(nnr)
        !
        REAL(DP) :: aux(nnr)
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%potential(update, local_verbose)
        !
        aux = env%dvtot%of_r
        !
#if defined(__MPI)
        IF (PRESENT(lgather)) THEN
            IF (lgather) THEN
                CALL env_gather_grid(env%system_cell%dfft, aux, dvtot)
            ELSE
                dvtot = aux
            END IF
        ELSE
            dvtot = aux
        END IF
#else
        dvtot = aux
#endif
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_potential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_energy(total_energy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(INOUT) :: total_energy
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%energy(total_energy)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_energy
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_force(nat, force_environ)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        !
        REAL(DP), INTENT(INOUT) :: force_environ(3, nat)
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%force(nat, force_environ)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_force
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_denergy(total_energy)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(INOUT) :: total_energy
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%denergy(total_energy)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_denergy
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE calc_interface
!----------------------------------------------------------------------------------------
