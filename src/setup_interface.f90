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
    USE env_char_ops, ONLY: env_uppercase
    !
    USE env_base_io, ONLY: prog, ionode, ionode_id, comm, program_unit, environ_unit, &
                           lstdout, verbose_ => verbose
    !
    USE env_io, ONLY: env_find_free_unit
    !
    USE environ_param, ONLY: DP, e2
    !
    USE environ_input, ONLY: read_environ_input
    !
    USE env_base_input
    !
    USE class_environ, ONLY: env
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: init_io, init_base_first, init_base_second, environ_clean
    !
    !------------------------------------------------------------------------------------

    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               INITIALIZATION METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Set global I/O constants
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_io(ionode_, ionode_id_, comm_, program_unit_)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: ionode_
        INTEGER, INTENT(IN) :: ionode_id_
        INTEGER, INTENT(IN) :: comm_
        INTEGER, INTENT(IN) :: program_unit_
        !
        !--------------------------------------------------------------------------------
        !
        ionode = ionode_
        ionode_id = ionode_id_
        comm = comm_
        !
        program_unit = program_unit_
        environ_unit = env_find_free_unit()
        !
        prog = 'PW'
        !
        lstdout = .TRUE. .AND. ionode
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_io
    !------------------------------------------------------------------------------------
    !>
    !! Set global base
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_base_first(nelec, nat, ntyp, atom_label, use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: use_internal_pbc_corr
        INTEGER, INTENT(IN) :: nelec, nat, ntyp
        INTEGER, INTENT(IN) :: atom_label(3, ntyp)
        !
        INTEGER :: i, j
        !
        CHARACTER(LEN=3) :: atom_chr(ntyp)
        !
        CHARACTER(LEN=80) :: sub_name = 'init_base_first'
        !
        !--------------------------------------------------------------------------------
        !
        CALL read_environ_input() ! read namelists and cards from environ.in
        !
        verbose_ = verbose ! set internal verbosity from input
        !
        !--------------------------------------------------------------------------------
        ! # TODO move this out into a util function or move into python
        !
        DO i = 1, ntyp
            !
            DO j = 1, 3
                atom_chr(i) (j:j + 1) = CHAR(atom_label(j, i))
            END DO
            !
        END DO
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%init_first(nelec, nat, ntyp, atom_chr, use_internal_pbc_corr)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_base_first
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_base_second(alat, at, comm_in, gcutm, e2_in)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm_in
        REAL(DP), INTENT(IN) :: alat
        REAL(DP), INTENT(IN) :: at(3, 3)
        REAL(DP), INTENT(IN) :: gcutm
        REAL(DP), OPTIONAL, INTENT(IN) :: e2_in
        !
        REAL(DP), ALLOCATABLE :: at_scaled(:, :)
        REAL(DP) :: gcutm_scaled
        !
        CHARACTER(LEN=80) :: sub_name = 'init_base_second'
        !
        !--------------------------------------------------------------------------------
        ! Allocate buffers used by env_mp_sum
        !
        CALL env_allocate_mp_buffers()
        !
        !--------------------------------------------------------------------------------
        ! PW uses Ryderg units (2.D0 * AU)
        ! CP uses Hartree units (e2_in = 1.D0)
        !
        IF (PRESENT(e2_in)) THEN
            e2 = e2_in
        ELSE
            e2 = 2.D0
        END IF
        !
        !--------------------------------------------------------------------------------
        !
        IF (alat < 1.D-8) CALL env_errore(sub_name, 'Wrong alat', 1)
        !
        IF (alat < 1.0_DP) CALL env_warning('strange lattice parameter')
        !
        ALLOCATE (at_scaled(3, 3))
        at_scaled = at * alat
        gcutm_scaled = gcutm / alat**2
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%init_second(at_scaled, comm_in, gcutm_scaled)
        !
        DEALLOCATE (at_scaled)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_base_second
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  CLEANING METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Clean up all the Environ related allocated variables, and call
    !! clean up subroutines of specific Environ modules.
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_clean(lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        !--------------------------------------------------------------------------------
        !
        CALL environ_clean_first(lflag)
        !
        CALL environ_clean_second(lflag)
        !
        CALL env_deallocate_mp_buffers()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_clean
    !------------------------------------------------------------------------------------
    !>
    !! Clean up all the Environ related allocated variables, and call clean up
    !! subroutines of specific Environ modules.
    !!
    !! The structure of this subroutine mirrors the one of init_environ subroutines
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_clean_first(lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        !--------------------------------------------------------------------------------
        ! Deallocate environment variables
        !
        IF (ASSOCIATED(env%vzero%cell)) CALL env%vzero%destroy()
        !
        IF (ASSOCIATED(env%dvtot%cell)) CALL env%dvtot%destroy()
        !
        !--------------------------------------------------------------------------------
        ! base_environ variables
        !
        IF (env%lelectrostatic .AND. ASSOCIATED(env%vreference%cell)) &
            CALL env%vreference%destroy()
        !
        IF (env%lsoftcavity .AND. ASSOCIATED(env%vsoftcavity%cell)) &
            CALL env%vsoftcavity%destroy()
        !
        IF (env%lconfine .AND. ASSOCIATED(env%vconfine%cell)) CALL env%vconfine%destroy()
        !
        !--------------------------------------------------------------------------------
        ! Destroy derived types which were allocated in input
        !
        IF (env%lelectrostatic .OR. env%lconfine) THEN
            !
            CALL env%system_charges%destroy(lflag)
            !
            CALL env%environment_charges%destroy(lflag)
            !
        END IF
        !
        IF (env%lexternals) CALL env%externals%destroy(lflag)
        !
        IF (env%lstatic) CALL env%static%destroy(lflag)
        !
        IF (env%lelectrolyte) CALL env%electrolyte%destroy(lflag)
        !
        IF (env%lsemiconductor) CALL env%semiconductor%destroy(lflag)
        !
        IF (env%laddcharges) CALL env%additional_charges%destroy()
        !
        CALL env%system_electrons%destroy(lflag)
        !
        CALL env%system_ions%destroy(lflag)
        !
        CALL env%system_system%destroy(lflag)
        !
        CALL env%environment_electrons%destroy(lflag)
        !
        CALL env%environment_ions%destroy(lflag)
        !
        CALL env%environment_system%destroy(lflag)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_clean_first
    !------------------------------------------------------------------------------------
    !>
    !! Clean up all the Environ-related allocated variables and call clean up
    !! subroutines of specific Environ modules. These are quantities that may
    !! be needed by TDDFPT, thus may need to be cleaned later
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE environ_clean_second(lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        LOGICAL :: opnd
        !
        !--------------------------------------------------------------------------------
        !
        IF (lflag) THEN
            INQUIRE (unit=environ_unit, opened=opnd)
            !
            IF (opnd) CLOSE (unit=environ_unit)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! base_environ variables
        !
        IF (env%lelectrostatic) THEN
            !
            IF (ASSOCIATED(env%velectrostatic%cell)) CALL env%velectrostatic%destroy()
            !
            CALL electrostatic_clean(lflag)
            !
        END IF
        !
        !--------------------------------------------------------------------------------
        ! Destroy derived types which were allocated in input
        !
        IF (env%loptical) THEN
            !
            CALL env%environment_response_charges%destroy(lflag)
            !
            CALL env%environment_response_electrons%destroy(lflag)
            !
            CALL env%system_response_charges%destroy(lflag)
            !
            CALL env%system_response_electrons%destroy(lflag)
            !
            CALL env%optical%destroy(lflag)
            !
        END IF
        !
        IF (env%lsolvent) CALL env%solvent%destroy(lflag)
        !
        IF (env%lboundary) CALL env%derivatives%destroy(lflag)
        !
        IF (env%ldoublecell) THEN
            !
            CALL env%mapping%destroy(lflag)
            !
            CALL env%environment_cell%destroy(lflag)
            !
        END IF
        !
        CALL env%system_cell%destroy(lflag)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE environ_clean_second
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE electrostatic_clean(lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: lflag
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%outer%destroy(lflag)
        !
        CALL env%reference%destroy(lflag)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE electrostatic_clean
    !------------------------------------------------------------------------------------
    !
    !
    !------------------------------------------------------------------------------------
END MODULE setup_interface
!----------------------------------------------------------------------------------------
