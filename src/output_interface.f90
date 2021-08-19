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
MODULE output_interface
    !------------------------------------------------------------------------------------
    !
    USE environ_param, ONLY: BOHR_RADIUS_SI, RYDBERG_SI, RYTOEV
    !
    USE env_base_io, ONLY: prog, ionode, program_unit
    !
    USE class_environ, ONLY: env
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: print_energies, print_potential_shift, &
              print_potential_warning, print_summary, &
              print_clocks, update_output_program_unit
    !------------------------------------------------------------------------------------
    
    !
    !------------------------------------------------------------------------------------
CONTAINS
!------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   OUTPUT METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Write out the different Environ contributions to the energy.
    !! Called by electrons.f90
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_energies()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=80) :: sub_name = 'print_energies'
        !
        !--------------------------------------------------------------------------------
        !
        IF (ionode) THEN
            !
            SELECT CASE (prog)
                !
            CASE ('PW')
                !
                IF (env%lelectrostatic) WRITE (program_unit, 1000) env%eelectrostatic
                !
                IF (env%lsurface) WRITE (program_unit, 1001) env%esurface
                !
                IF (env%lvolume) WRITE (program_unit, 1002) env%evolume
                !
                IF (env%lconfine) WRITE (program_unit, 1003) env%econfine
                !
                IF (env%lelectrolyte) WRITE (program_unit, 1004) env%eelectrolyte
                !
                WRITE (program_unit, 1005) env%deenviron
                !
            CASE DEFAULT
                CALL env_errore(sub_name, 'Wrong program calling Environ', 1)
                !
            END SELECT
            !
        END IF
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT('     electrostatic embedding   =', F17.8, ' Ry')
1001    FORMAT('     cavitation energy         =', F17.8, ' Ry')
1002    FORMAT('     PV energy                 =', F17.8, ' Ry')
1003    FORMAT('     confinement energy        =', F17.8, ' Ry')
1004    FORMAT('     electrolyte free energy   =', F17.8, ' Ry')
1005    FORMAT('     correction to one-el term =', F17.8, ' Ry')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_energies
    !------------------------------------------------------------------------------------
    !>
    !! If Gaussian nuclei are used, write out the corresponding potential shift
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_potential_shift()
        !--------------------------------------------------------------------------------
        !
        IF (env%lsmearedions) &
            WRITE (program_unit, 1100) env%environment_ions%potential_shift * RYTOEV
        !
1100    FORMAT(/, 5(' '), &
                'the potential shift due to the Gaussian-smeared nuclei is ', &
                F10.4, ' ev')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_potential_shift
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_potential_warning()
        !--------------------------------------------------------------------------------
        !
        IF (env%need_pbc_correction) WRITE (program_unit, 1200)
        !
1200    FORMAT(/, &
                5(' '), 'WARNING: you are using the parabolic pbc correction;', /, &
                5(' '), '         the potential shift above must be added to ', /, &
                5(' '), '         band and Fermi energies.')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_potential_warning
    !------------------------------------------------------------------------------------
    !!
    !> Write out the main parameters of Environ calculations, summarizing
    !! the input keywords (some info also on internal vs input units).
    !! Called by summary.f90
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_summary()
        !--------------------------------------------------------------------------------
        !
        IF (ionode .AND. prog == 'PW') CALL env%summary()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_summary
    !------------------------------------------------------------------------------------
    !>
    !! Write out the time informations of the Environ dependent calculations.
    !! Called by print_clock_pw.f90
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_clocks(passed_unit)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN), OPTIONAL :: passed_unit
        !
        INTEGER :: actual_unit
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(passed_unit)) THEN
            actual_unit = passed_unit
        ELSE
            actual_unit = program_unit
        END IF
        !
        WRITE (actual_unit, *)
        WRITE (actual_unit, '(5X,"Environ routines")')
        !
        !--------------------------------------------------------------------------------
        ! Dielectric subroutines
        !
        IF (env%lelectrostatic) THEN
            !
            CALL env_print_clock('calc_eelect')
            !
            CALL env_print_clock('calc_velect')
            !
            CALL env_print_clock('calc_vgcs')
            !
            CALL env_print_clock('dielectric')
            !
            CALL env_print_clock('electrolyte')
            !
            CALL env_print_clock('calc_felect')
            !
        END IF
        !
        IF (env%lsemiconductor) CALL env_print_clock('calc_vms')
        !
        !--------------------------------------------------------------------------------
        ! TDDFT
        !
        IF (env%ltddfpt) CALL env_print_clock('calc_vsolvent_tddfpt')
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_clocks
    !------------------------------------------------------------------------------------
    !>
    !! Sets the output file target #TODO do we need this routine?
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_output_program_unit(program_unit_in)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: program_unit_in
        !
        !--------------------------------------------------------------------------------
        !
        program_unit = program_unit_in
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_output_program_unit
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE output_interface
!----------------------------------------------------------------------------------------
