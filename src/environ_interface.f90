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
MODULE environ_interface
    !------------------------------------------------------------------------------------
    !
    USE class_io, ONLY: io
    !
    USE env_base_scatter, ONLY: env_scatter_grid, env_gather_grid
    USE env_types_fft, ONLY: env_fft_type_descriptor
    !
    USE environ_param, ONLY: DP
    !
    USE environ_input, ONLY: read_environ_input
    USE env_base_input
    !
    USE env_global_objects, ONLY: env, setup
    USE class_destructor, ONLY: clean
    USE class_calculator, ONLY: calc
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: init_io, read_input, init_environ, clean_environ
    !
    PUBLIC :: update_cell, update_ions, update_electrons, update_response

    PUBLIC :: add_mbx_charges
    !
    PUBLIC :: set_restart, is_restart, get_threshold, get_nnt, get_nrx
    !
    PUBLIC :: calc_potential, calc_energy, calc_force, calc_denergy
    !
    PUBLIC :: print_energies, print_potential_shift, print_potential_warning, &
              print_summary, print_clocks
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   SETUP METHODS
    !
    !------------------------------------------------------------------------------------
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
                            use_internal_pbc_corr, at, gcutm)
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
        REAL(DP), INTENT(IN) :: at(3, 3)
        REAL(DP), INTENT(IN) :: gcutm
        !
        CHARACTER(LEN=80) :: sub_name = 'init_environ'
        !
        !--------------------------------------------------------------------------------
        !
        CALL setup%init() ! flags and solvers
        !
        CALL env_allocate_mp_buffers()
        !
        CALL setup%init_cell(gcutm, comm, at)
        !
        CALL setup%init_cores(gcutm, use_internal_pbc_corr)
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
    !------------------------------------------------------------------------------------
    !
    !                                   UPDATE METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_cell(at)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: at(3, 3)
        !
        !--------------------------------------------------------------------------------
        !
        CALL setup%update_cell(at)
        !
        CALL env%update_cell_dependent_quantities()
        !
        CALL setup%end_cell_update()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_cell
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_ions(nat, tau)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: nat
        REAL(DP), INTENT(IN) :: tau(3, nat)
        !
        !--------------------------------------------------------------------------------
        !
        CALL env%update_ions(nat, tau)
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
    !------------------------------------------------------------------------------------
    !
    !                                   ACCESS METHODS
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE set_restart(flag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, INTENT(IN) :: flag
        !
        !--------------------------------------------------------------------------------
        !
        setup%restart = flag
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE set_restart
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION is_restart()
        !--------------------------------------------------------------------------------
        !
        is_restart = setup%restart
        !
        !--------------------------------------------------------------------------------
    END FUNCTION is_restart
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    REAL(DP) FUNCTION get_threshold()
        !--------------------------------------------------------------------------------
        !
        get_threshold = setup%threshold
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_threshold
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION get_nnt()
        !--------------------------------------------------------------------------------
        !
        get_nnt = setup%system_cell%dfft%nnt
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_nnt
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    INTEGER FUNCTION get_nrx(i)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: i
        !
        !--------------------------------------------------------------------------------
        !
        SELECT CASE (i)
            !
        CASE (0)
            get_nrx = setup%system_cell%dfft%nr1x
            !
        CASE (1)
            get_nrx = setup%system_cell%dfft%nr2x
            !
        CASE (2)
            get_nrx = setup%system_cell%dfft%nr3x
            !
        END SELECT
        !
        !--------------------------------------------------------------------------------
    END FUNCTION get_nrx
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
        CALL env%print_energies('PW', .FALSE.)
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
END MODULE environ_interface
!----------------------------------------------------------------------------------------
