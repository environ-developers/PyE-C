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
! Authors: Oliviero Andreussi (Department of Physics, UNT)
!          Edan Bainglass     (Department of Physics, UNT)
!
!----------------------------------------------------------------------------------------
!>
!! An interface module for all communications between Environ and the calling program
!!
!----------------------------------------------------------------------------------------
MODULE pyec_interface
    !------------------------------------------------------------------------------------
    !
    USE env_char_ops, ONLY: env_uppercase
    !
    USE env_base_io, ONLY: prog, ionode, ionode_id, comm, program_unit, environ_unit, &
                           lstdout, depth, verbose_ => verbose
    !
    USE env_io, ONLY: env_find_free_unit
    !
    USE environ_param, ONLY: DP, BOHR_RADIUS_SI, RYDBERG_SI, RYTOEV, e2
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
    PUBLIC :: init_io, init_base_first, init_base_second, &
              init_potential, init_cell, init_ions, &
              init_electrons, init_response
    !
    PUBLIC :: calc_potential, calc_energy, calc_force, calc_denergy
    !
    PUBLIC :: print_energies, print_potential_shift, &
              print_potential_warning, print_summary, &
              print_clocks, update_output_program_unit
    !
    PUBLIC :: environ_clean, environ_clean_first, environ_clean_second
    !
    PUBLIC :: is_tddfpt, is_restart
    !
    PUBLIC :: get_threshold, set_restart, get_nskip, &
              get_dvtot_of_r, get_verbose
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
    SUBROUTINE init_io(prog_, ionode_, ionode_id_, comm_, program_unit_)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: prog_
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
        prog = env_uppercase(prog_(1:2))
        !
        SELECT CASE (prog)
            !
        CASE ('PW', 'CP')
            lstdout = .TRUE.
            !
        CASE DEFAULT
            lstdout = .FALSE.
            !
        END SELECT
        !
        lstdout = lstdout .AND. ionode
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_io
    !------------------------------------------------------------------------------------
    !>
    !! Set global base
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_base_first(nelec, nat, ntyp, atom_label, &
                                       use_internal_pbc_corr)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: use_internal_pbc_corr
        INTEGER, INTENT(IN) :: nelec, nat, ntyp
        INTEGER, INTENT(IN) :: atom_label(3, ntyp)
        CHARACTER(LEN=3) :: atom_chr(ntyp)
        !
        INTEGER :: i, j
        !
        CHARACTER(LEN=80) :: sub_name = 'init_base_first'
        !
        ! //TODO move this out into a util function or move into python
        DO i = 1, ntyp
            DO j = 1, 3
                atom_chr(i)(j:j+1) = CHAR(atom_label(j, i))
            ENDDO
        ENDDO
        !
        !--------------------------------------------------------------------------------
        !
        CALL read_environ_input() ! read namelists and cards from environ.in
        !
        verbose_ = verbose ! set internal verbosity from input
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
    SUBROUTINE init_base_second(alat, at, comm_in, me, root, gcutm, e2_in)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: comm_in, me, root
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
        ALLOCATE(at_scaled(3, 3))
        at_scaled = at * alat
        gcutm_scaled = gcutm / alat**2
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%init_second(at_scaled, comm_in, me, root, gcutm_scaled)
        !
        DEALLOCATE(at_scaled)
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_base_second
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
    !------------------------------------------------------------------------------------
    !
    !                               COMPUTATION METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_potential(update, nnr, dvtot, local_verbose)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: update
        INTEGER, INTENT(IN) :: nnr
        INTEGER, INTENT(IN), OPTIONAL :: local_verbose
        !
        REAL(DP), INTENT(OUT) :: dvtot(nnr)
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%potential(update, local_verbose)
        !
        dvtot = env%dvtot%of_r
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
    !                                 ACCESS METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_restart(environ_restart)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: environ_restart
        !
        !--------------------------------------------------------------------------------
        !
        env%environ_restart = environ_restart
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_restart
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION get_threshold()
        !--------------------------------------------------------------------------------
        !
        get_threshold = env%environ_thr
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_threshold
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION get_nskip()
        !--------------------------------------------------------------------------------
        !
        get_nskip = env%environ_nskip
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_nskip
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    FUNCTION get_dvtot_of_r(nnr) RESULT(dvtot_of_r)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nnr
        !
        REAL(DP) :: dvtot_of_r(nnr)
        !
        !--------------------------------------------------------------------------------
        !
        dvtot_of_r = env%dvtot%of_r
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_dvtot_of_r
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION get_verbose()
        !--------------------------------------------------------------------------------
        !
        get_verbose = verbose_
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_verbose
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION is_restart()
        !--------------------------------------------------------------------------------
        !
        is_restart = env%environ_restart
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION is_restart
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION is_tddfpt()
        !--------------------------------------------------------------------------------
        !
        is_tddfpt = env%ltddfpt
        !
        RETURN
        !
        !--------------------------------------------------------------------------------
    END FUNCTION is_tddfpt
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
            CASE ('CP')
                !
                IF (env%lelectrostatic) WRITE (program_unit, 1006) env%eelectrostatic
                !
                IF (env%lsurface) WRITE (program_unit, 1007) env%esurface
                !
                IF (env%lvolume) WRITE (program_unit, 1008) env%evolume
                !
                IF (env%lconfine) WRITE (program_unit, 1009) env%econfine
                !
                IF (env%lelectrolyte) WRITE (program_unit, 1010) env%eelectrolyte
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
1006    FORMAT('     electrostatic embedding = ', F14.5, ' Hartree a.u.')
1007    FORMAT('           cavitation energy = ', F14.5, ' Hartree a.u.')
1008    FORMAT('                   PV energy = ', F14.5, ' Hartree a.u.')
1009    FORMAT('     electrolyte free energy = ', F14.5, ' Hartree a.u.')
1010    FORMAT('          confinement energy = ', F14.5, ' Hartree a.u.')
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
        RETURN
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
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_output_program_unit
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                               CLEANING METHODS
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
        RETURN
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
        RETURN
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
        RETURN
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
        RETURN
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE electrostatic_clean
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE pyec_interface
!----------------------------------------------------------------------------------------
