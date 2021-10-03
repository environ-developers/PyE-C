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
MODULE calc_interface
    !------------------------------------------------------------------------------------
    !
    USE env_base_scatter, ONLY: env_gather_grid
    !
    USE environ_param, ONLY: DP
    !
    USE env_global_objects, ONLY: env, setup
    USE class_calculator, ONLY: calc
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: calc_potential, calc_energy, calc_force, calc_denergy
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                COMPUTATION METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_potential(update, potential, local_verbose, lgather)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: update
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        LOGICAL, INTENT(IN), OPTIONAL :: lgather
        !
        REAL(DP), INTENT(OUT) :: potential(setup%system_cell%dfft%nnt)
        !
        REAL(DP) :: aux(setup%system_cell%dfft%nnr)
        !
        !--------------------------------------------------------------------------------
        !
        CALL calc%potential(env, update, local_verbose)
        !
        aux = env%dvtot%of_r
        !
#if defined(__MPI)
        IF (PRESENT(lgather)) THEN
            !
            IF (lgather) THEN
                CALL env_gather_grid(setup%system_cell%dfft, aux, potential)
            ELSE
                potential = aux
            END IF
            !
        ELSE
            potential = aux
        END IF
#else
        !
        potential = aux
#endif
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
        REAL(DP), INTENT(OUT) :: total_energy
        !
        !--------------------------------------------------------------------------------
        !
        CALL calc%energy(env, total_energy)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_energy
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_force(force_environ)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(OUT) :: force_environ(3, env%system_ions%number)
        !
        !--------------------------------------------------------------------------------
        !
        CALL calc%force(env, env%system_ions%number, force_environ)
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
        REAL(DP), INTENT(OUT) :: total_energy
        !
        !--------------------------------------------------------------------------------
        !
        CALL calc%denergy(env, total_energy)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_denergy
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE calc_interface
!----------------------------------------------------------------------------------------
