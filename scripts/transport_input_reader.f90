!==============================================================================
! Simple reader for transport.txt
!==============================================================================
module transport_input_reader
    implicit none

    integer, parameter :: MAX_NUCLIDES = 32
    integer, parameter :: MAX_XS = 32
    integer, parameter :: MAX_OTF = 32
    integer, parameter :: PATH_LEN = 100

    integer, parameter :: XS_TABLE = 0
    integer, parameter :: XS_NJOY = 1
    integer, parameter :: XS_OTF = 2
    integer, parameter :: XS_TMS = 3

    type :: nuclide_input
        character(len=64) :: name = ''

        integer :: n_xs = 0
        real :: xs_temperature(MAX_XS) = 0.0
        character(len=PATH_LEN) :: xs_file(MAX_XS) = ''

        character(len=PATH_LEN) :: otf_grid_file = ''
        integer :: n_otf = 0
        integer :: otf_mt(MAX_OTF) = 0
        character(len=PATH_LEN) :: otf_coeff_file(MAX_OTF) = ''
    end type nuclide_input

    type :: transport_input_t
        integer :: run_mode = 0
        integer :: xs_method = XS_TABLE

        integer :: histories = 100000
        integer :: maxwell_grid_points = 500

        real :: source_energy_ev = 100.0
        real :: energy_cutoff_ev = 1.0
        real :: doppler_table_temperature_k = 3200.0

        integer :: n_nuclides = 0
        type(nuclide_input) :: nuclide(MAX_NUCLIDES)
    end type transport_input_t

contains

    subroutine read_transport_input(filename, input)
        character(len=*), intent(in) :: filename
        type(transport_input_t), intent(out) :: input

        integer :: unit_num, ios, p
        integer :: current_nuclide
        integer :: mt
        real :: temperature

        character(len=1024) :: line
        character(len=512) :: key, value

        current_nuclide = 0

        open(newunit=unit_num, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'ERROR: cannot open input file: ', trim(filename)
            error stop
        end if

        do
            read(unit_num, '(A)', iostat=ios) line
            

            line = adjustl(trim(line))
            ! Empty line or comment.
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle
            if(trim(line) == 'END_PROGRAM')exit
            !----------------------------------------------------------------------
            ! Start of a new nuclide card.
            !----------------------------------------------------------------------
            if (index(line, 'BEGIN_NUCLIDE') == 1) then
                input%n_nuclides = input%n_nuclides + 1
                current_nuclide = input%n_nuclides

                input%nuclide(current_nuclide)%name = &
                    adjustl(trim(line(len('BEGIN_NUCLIDE') + 1:)))

                cycle
            end if 

            ! End of current nuclide card.
            if (trim(line) == 'END_NUCLIDE') then
                current_nuclide = 0
                cycle
            end if

            ! All remaining input lines have KEY = VALUE form.
            p = index(line, '=')
            if (p == 0) cycle

            key   = adjustl(trim(line(:p-1)))
            value = adjustl(trim(line(p+1:)))

            !======================================================================
            ! Global input parameters.
            !======================================================================
            if (current_nuclide == 0) then

                select case (trim(key))

                case ('RUN_MODE')
                    select case (trim(value))
                    case ('TRANSPORT')
                        input%run_mode = 0
                    case ('DOPPLER_TABLE')
                        input%run_mode = 1
                    case ('OTF_PREPROCESS')
                        input%run_mode = 2
                    case ('SPATIAL_TRANSPORT')
                        input%run_mode = 3
                    end select

                case ('XS_METHOD')
                    select case (trim(value))
                    case ('TABLE')
                        input%xs_method = XS_TABLE
                    case ('NJOY')
                        input%xs_method = XS_NJOY
                    case ('OTF')
                        input%xs_method = XS_OTF
                    case ('TMS')
                        input%xs_method = XS_TMS
                    end select

                case ('HISTORIES')
                    read(value, *) input%histories

                case ('SOURCE_ENERGY_EV')
                    read(value, *) input%source_energy_ev

                case ('ENERGY_CUTOFF_EV')
                    read(value, *) input%energy_cutoff_ev

                case ('DOPPLER_TABLE_TEMPERATURE_K')
                    read(value, *) input%doppler_table_temperature_k

                case ('MAXWELL_GRID_POINTS')
                    read(value, *) input%maxwell_grid_points

                end select

            !======================================================================
            ! Current nuclide card.
            !======================================================================
            else

                if (index(key, 'XS ') == 1) then
                    input%nuclide(current_nuclide)%n_xs = input%nuclide(current_nuclide)%n_xs + 1

                    read(key(3:), *) temperature

                    p = input%nuclide(current_nuclide)%n_xs

                    input%nuclide(current_nuclide)%xs_temperature(p) = temperature
                    input%nuclide(current_nuclide)%xs_file(p) = trim(value)

                else if (trim(key) == 'OTF_GRID') then
                    input%nuclide(current_nuclide)%otf_grid_file = trim(value)

                else if (index(key, 'OTF_MT ') == 1) then
                    input%nuclide(current_nuclide)%n_otf = &
                        input%nuclide(current_nuclide)%n_otf + 1

                    read(key(8:), *) mt

                    p = input%nuclide(current_nuclide)%n_otf

                    input%nuclide(current_nuclide)%otf_mt(p) = mt
                    input%nuclide(current_nuclide)%otf_coeff_file(p) = trim(value)
                end if

            end if

        end do

        close(unit_num)
    end subroutine read_transport_input

end module transport_input_reader
