program MCNP_builder
    use MCNPconfigOTF
    use DoplerTrearmentNJOY
    use Dopler_MCNP_OTF
    implicit none

    type(otf_config_type) :: cfg

    character(len=256) :: filename
    character(len=64)  :: material_name

    integer :: nreac, npoints
    integer, allocatable :: mt(:)

    real, allocatable :: E(:)
    real, allocatable :: sigma(:,:)

    real, allocatable :: Tunion(:)
    real, allocatable :: Tfit(:)
    real, allocatable :: Egrid(:)

    real, allocatable :: coeff(:,:,:)
    real, allocatable :: max_error(:,:)

    ! Cross sections reconstructed ONLY from fitted coefficients
    ! on the final union energy grid.
    real, allocatable :: sigma_verify(:,:)

    integer :: i, g, cidx, unit_coeff, unit_all, unit_verify
    integer, parameter :: fit_order = 8
    real, parameter :: verify_temperature = 1200.0
    character(len=256) :: coeff_filename

    ! ------------------------------------------------------------
    ! Input file name:
    !   .\MCNP_builder.exe Cross-section-data-Carbon-294.txt
    ! ------------------------------------------------------------
    if (command_argument_count() < 1) then
        print *, "Usage:"
        print *, "  MCNP_builder.exe cross_section_file.txt"
        stop
    end if

    call get_command_argument(1, filename)

    ! ============================================================
    ! THE ONLY NEW PROCEDURE IN THIS PROGRAM:
    ! read the cross-section file and fill E and sigma arrays.
    ! ============================================================
    call load_cross_sections(trim(filename), material_name, mt, E, sigma, nreac, npoints)

    print *, "Material      : ", trim(material_name)
    print *, "Reactions     : ", nreac
    print *, "Energy points : ", npoints
    print *, "MT numbers    : ", mt

    ! ------------------------------------------------------------
    ! Initial energy grid is the original energy grid from the file.
    ! build_union_grid() from Dopler_MCNP_OTF will refine it.
    ! ------------------------------------------------------------
    allocate(Egrid(npoints))
    Egrid = E

    ! ------------------------------------------------------------
    ! Existing function from Dopler_MCNP_OTF:
    ! temperature grid for union-grid construction.
    ! ------------------------------------------------------------
    call build_temerature_grid(cfg%t_min, cfg%t_max, cfg%dt_union, Tunion)

    ! ------------------------------------------------------------
    ! Existing function from Dopler_MCNP_OTF:
    ! fitting-temperature grid.
    ! ------------------------------------------------------------
    call build_temerature_grid(cfg%t_min, cfg%t_max, cfg%dt_fit, Tfit)

    ! ------------------------------------------------------------
    ! Existing build_union_grid() currently accepts one sigma(:)
    ! at a time. Apply it successively to each reaction.
    !
    ! The same Egrid is refined in-place, so after all reactions
    ! it contains the union of the points required by every MT.
    ! ------------------------------------------------------------
    do i = 1, nreac
        print *, "Building union grid for MT = ", mt(i)

        call build_union_grid( &
            Egrid,        &
            Tunion,       &
            1,            &
            cfg%ft,       &
            E,            &
            sigma(:,i),   &
            npoints       &
        )
    end do

    print *, "Final union energy-grid size = ", size(Egrid)

    ! Storage for cross sections reconstructed from the fitted
    ! coefficients at verify_temperature.
    allocate(sigma_verify(size(Egrid), nreac))
    sigma_verify = 0.0

    ! ------------------------------------------------------------
    ! Existing build_coefficients() operates on sigma(:), therefore
    ! build coefficients reaction-by-reaction.
    !
    ! If you later want all reactions in one coeff array, the old
    ! build_coefficients interface itself should be generalized.
    ! No additional helper function is introduced here.
    ! ------------------------------------------------------------
    open(newunit=unit_all, file="coefficients_all.txt", status="replace", action="write")
    write(unit_all,'(A)') "# MT index energy C1 C2 C3 C4 C5 max_error"

    do i = 1, nreac
        print *, "Building coefficients for MT = ", mt(i)

        if (allocated(coeff)) deallocate(coeff)
        if (allocated(max_error)) deallocate(max_error)

        call build_coefficients( &
            Egrid,        &
            Tfit,         &
            1,            &
            fit_order,    &
            coeff,        &
            max_error,    &
            E,            &
            sigma(:,i)    &
        )

        print *, "MT ", mt(i), " max fit error = ", maxval(max_error)

        ! --------------------------------------------------------
        ! VERIFICATION:
        ! Reconstruct sigma(E,T) using ONLY the fitted coefficients
        ! and the existing eval_expansion() function.
        !
        ! No direct doplerBroadr() call is used here.
        ! --------------------------------------------------------
        do g = 1, size(Egrid)
            sigma_verify(g,i) = eval_expansion( &
                verify_temperature, &
                minval(Tfit),       &
                maxval(Tfit),       &
                fit_order,          &
                coeff(:,1,g)        &
            )
        end do

        ! --------------------------------------------------------
        ! Output fitted coefficients for every energy-grid point.
        ! coeff(:,1,g) because build_coefficients is called here
        ! for one reaction at a time (nreac = 1).
        ! --------------------------------------------------------
        write(coeff_filename,'("coefficients_MT_",I0,".txt")') mt(i)
        open(newunit=unit_coeff, file=trim(coeff_filename), status="replace", action="write")

        write(unit_coeff,'(A)',advance='no') "# index energy"
        do cidx = 1, size(coeff,1)
            write(unit_coeff,'(A,I0)',advance='no') " C", cidx
        end do
        write(unit_coeff,'(A)') " max_error"

        do g = 1, size(Egrid)
            write(*,'(A,I0,A,I0,A,ES14.6)',advance='no') &
                "MT=", mt(i), "  point=", g, "  E=", Egrid(g)

            do cidx = 1, size(coeff,1)
                write(*,'(A,I0,A,ES14.6)',advance='no') &
                    "  C", cidx, "=", coeff(cidx,1,g)
            end do
            write(*,'(A,ES14.6)') "  err=", max_error(1,g)

            write(unit_coeff,'(I8,1X,ES24.16)',advance='no') g, Egrid(g)
            do cidx = 1, size(coeff,1)
                write(unit_coeff,'(1X,ES24.16)',advance='no') coeff(cidx,1,g)
            end do
            write(unit_coeff,'(1X,ES24.16)') max_error(1,g)

            write(unit_all,'(I6,1X,I8,1X,ES24.16)',advance='no') mt(i), g, Egrid(g)
            do cidx = 1, size(coeff,1)
                write(unit_all,'(1X,ES24.16)',advance='no') coeff(cidx,1,g)
            end do
            write(unit_all,'(1X,ES24.16)') max_error(1,g)
        end do

        close(unit_coeff)
        print *, "Saved: ", trim(coeff_filename)
    end do

    close(unit_all)
    print *, "Saved: coefficients_all.txt"

    ! ------------------------------------------------------------
    ! Save verification cross sections reconstructed from
    ! coefficients on the UNIQUE energy grid.
    !
    ! File format:
    ! energy  sigma_MT_1  sigma_MT_2  sigma_MT_102 ...
    ! ------------------------------------------------------------
    open(newunit=unit_verify, file="cross_sections_from_coefficients_T1200.txt", &
         status="replace", action="write")

    write(unit_verify,'(A)',advance='no') "# energy"
    do i = 1, nreac
        write(unit_verify,'(A,I0)',advance='no') " sigma_MT_", mt(i)
    end do
    write(unit_verify,*)

    do g = 1, size(Egrid)
        write(unit_verify,'(ES24.16)',advance='no') Egrid(g)
        do i = 1, nreac
            write(unit_verify,'(1X,ES24.16)',advance='no') sigma_verify(g,i)
        end do
        write(unit_verify,*)
    end do

    close(unit_verify)

    print *, "Verification temperature = ", verify_temperature, " K"
    print *, "Saved: cross_sections_from_coefficients_T1200.txt"

    ! ------------------------------------------------------------
    ! Save final union grid.
    ! This is ordinary main-program I/O, not another new procedure.
    ! ------------------------------------------------------------
    open(unit=20, file="unique_energy_grid.txt", status="replace", action="write")
    do i = 1, size(Egrid)
        write(20,'(ES24.16)') Egrid(i)
    end do
    close(20)

    print *, "Saved: unique_energy_grid.txt"

    print *, "=============================================="
    print *, "BOTH OUTPUTS CREATED IN THIS RUN:"
    print *, "  coefficients_all.txt"
    print *, "  cross_sections_from_coefficients_T1200.txt"
    print *, "=============================================="

contains

    ! ============================================================
    ! ONLY NEW PROCEDURE
    ! ============================================================
    subroutine load_cross_sections(filename, material_name, mt, E, sigma, nreac, npoints)
        implicit none

        character(len=*), intent(in)  :: filename
        character(len=*), intent(out) :: material_name

        integer, allocatable, intent(out) :: mt(:)
        real,    allocatable, intent(out) :: E(:)
        real,    allocatable, intent(out) :: sigma(:,:)

        integer, intent(out) :: nreac
        integer, intent(out) :: npoints

        integer :: unit_id, ios
        integer :: i
        integer :: header_value_1
        integer :: mass_number
        real    :: density

        open(newunit=unit_id, file=filename, status="old", action="read", iostat=ios)

        if (ios /= 0) then
            print *, "ERROR: cannot open file: ", trim(filename)
            error stop
        end if

        ! Expected input format, e.g.:
        !
        ! Carbon
        ! 2
        ! 12
        ! 0.08
        ! 3
        ! 185
        ! 1 2 102
        ! E sigma_MT1 sigma_MT2 sigma_MT102
        !
        read(unit_id, *, iostat=ios) material_name
        if (ios /= 0) error stop "ERROR reading material name"

        read(unit_id, *, iostat=ios) header_value_1
        if (ios /= 0) error stop "ERROR reading header value 1"

        read(unit_id, *, iostat=ios) mass_number
        if (ios /= 0) error stop "ERROR reading mass number"

        read(unit_id, *, iostat=ios) density
        if (ios /= 0) error stop "ERROR reading density"

        read(unit_id, *, iostat=ios) nreac
        if (ios /= 0 .or. nreac <= 0) error stop "ERROR reading number of reactions"

        read(unit_id, *, iostat=ios) npoints
        if (ios /= 0 .or. npoints <= 1) error stop "ERROR reading number of energy points"

        allocate(mt(nreac))
        allocate(E(npoints))
        allocate(sigma(npoints,nreac))

        read(unit_id, *, iostat=ios) mt
        if (ios /= 0) error stop "ERROR reading MT numbers"

        do i = 1, npoints
            read(unit_id, *, iostat=ios) E(i), sigma(i,:)

            if (ios /= 0) then
                print *, "ERROR reading cross-section row ", i
                error stop
            end if
        end do

        close(unit_id)

    end subroutine load_cross_sections

end program MCNP_builder
