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
MODULE output_interface
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE env_global_objects, ONLY: env, setup
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: print_energies, print_potential_shift, print_potential_warning, &
              print_summary, print_clocks
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
        CALL env%print_energies('PW')
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
        CALL env%print_potential_shift()
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
        CALL setup%print_potential_warning()
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
        CALL setup%print_summary()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_summary
    !------------------------------------------------------------------------------------
    !>
    !! Write out the time informations of the Environ dependent calculations.
    !! Called by print_clock_pw.f90
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_clocks()
        !--------------------------------------------------------------------------------
        !
        CALL setup%print_clocks()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_clocks
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE output_interface
!----------------------------------------------------------------------------------------
