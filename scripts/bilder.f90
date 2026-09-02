!==============================================================================
! Bilder Program: bilder
!
! Purpose:
!   Initializes nuclear data and OTF data, selects the requested operating mode,
!   launches neutron histories or preprocessing tasks, and reports simulation
!   statistics.
!
! Operating modes (typeWork):
!   0 - Monte Carlo transport using the base transport path.
!   1 - Generate Doppler-broadened cross-section tables.
!   2 - Build adaptive OTF union energy grids.
!   3 - Run the spatially varying-temperature transport study.
!
!==============================================================================

program bilder
    use data_type
    use math_operation
    use operation_with_data
    use process_manager
    use DoplerTrearmentNJOY
    use MCNPconfigOTF
    use Dopler_MCNP_OTF
    use transport_input_reader

    implicit none

    integer :: N = 100000                                               ! Number of source neutron histories.
    integer i,r,j                                                       ! General-purpose history/grid/reaction loop indices.
    integer :: tem = 3200                                               ! Temperature used by mode 1 table generation [K].
    real integ,alpha                                                    ! Maxwell-CDF integral and Doppler alpha parameter.
    ! Main executable mode selector. This remains a source-level setting to preserve
    ! the original workflow; change the value before compilation when another mode is required.
    integer :: typeWork = 0                                             ! 0=transport, 1=Doppler table generation, 2=OTF grid generation, 3=spatial-temperature study.
    real :: start_time, end_time, elapsed_time                          ! CPU timing values [s].
    real :: dispersion = 0                                              ! Variance-like estimator used in FOM calculation.
    real :: age_of_netron = 0                                           ! Neutron-age estimator
    real :: livedDistance = 0                                           ! Accumulated/mean radius of surviving neutron histories.
    integer :: averageCollision = 0                                     ! Accumulated collision count
    integer :: countLivedNetron = 0                                     ! Number of histories not absorbed at transport termination.
    real :: lifeTime = 0                                                ! Accumulated/mean neutron lifetime.
    integer :: count_point = 500                                        ! Number of target-speed lookup points.
    type(netron_data),allocatable :: netrons(:)                         ! Source/history array.
    type(enviroment) env                                                ! Material and nuclear-data environment.
    type(otf_config_type) conf                                          ! OTF preprocessing configuration.

    type(transport_input_t) :: input
    character(len=100), allocatable :: directory_of_cross_section_data(:, :)
    character(len=100), allocatable :: directory_of_coeff_data(:)

    integer :: n_nuclides
    integer :: n_temperatures
    integer :: n_reactions
    integer :: current_index_nuclear
    integer :: typeDopler

    real, allocatable :: Tgrid(:),TfitGrid(:)                           ! Legacy work grid plus union/fitting temperature grids.
    character(len=100) :: directory_of_coeff_data_Carbon(3)             ! Carbon OTF coefficient files by reaction.
    character(len=100) :: directory_of_coeff_data_Uranium(3)            ! Uranium OTF coefficient files by reaction.
    character(len=100) :: file_name                                     ! Generated output filename.
    real::  local_speed_estimator(8) = 0                                ! Eight spatially resolved local estimator accumulators.
    real:: local_count_lived(8) = 0                                     ! Counts used to normalize the local estimators.
    real :: localPos = 0.0                                              ! Current neutron radial distance from the origin.

    call read_transport_input('test.txt', input)   
    N = input%histories

    tem = nint(input%doppler_table_temperature_k)

    typeWork = input%run_mode
    typeDopler = input%xs_method

    count_point = input%maxwell_grid_points



    allocate(netrons(N))
    
    n_nuclides = input%n_nuclides

    n_temperatures = input%nuclide(1)%n_xs

    allocate(directory_of_cross_section_data(n_nuclides, n_temperatures))
    print*,n_nuclides,n_temperatures
    do i = 1, n_nuclides
        do j = 1, n_temperatures
            directory_of_cross_section_data(i,j) = adjustl(trim(input%nuclide(i)%xs_file(j)))
        end do
    end do
    print*,directory_of_cross_section_data(2,2)
    

    call load_cross_section_fron_file(n_nuclides,n_temperatures,directory_of_cross_section_data,env)
    
    do i = 1, n_nuclides
        n_reactions = input%nuclide(i)%n_otf
        print*,n_reactions
        allocate(directory_of_coeff_data(n_reactions))
        print*,size(directory_of_coeff_data)
        do j = 1, n_reactions
            if(env%different_tipe_of_nuclear(i)%name_of_nuclie == input%nuclide(j)%name) then
                current_index_nuclear = j
            end if
        end do

        do j = 1, n_reactions
            directory_of_coeff_data(j) = input%nuclide(i)%otf_coeff_file(j)
        end do
        
        call load_OTF_coefficients(directory_of_coeff_data,env%different_tipe_of_nuclear(current_index_nuclear))
        call load_unique_energy_grid(input%nuclide(i)%otf_grid_file,env%different_tipe_of_nuclear(current_index_nuclear))
        deallocate(directory_of_coeff_data)
    end do



    env%tem_grid = [294.0,400.0,800.0,1200.0,1800.0,2200.0,2800.0, 3200.0]
    
    select case(typeWork)
        case(0)
        !--------------------------------------------------------------------------
        ! Mode 0: Monte Carlo transport using the base transport path.
        !--------------------------------------------------------------------------
   
            ! Build Maxwell target-speed CDF tables for both loaded nuclides.
            do i = 1,2
            ! Выделяем место под массивы
                call get_mass(env%different_tipe_of_nuclear(i)%mass_of_nuclear * 1.660539e-27)
                allocate(env%different_tipe_of_nuclear(i)%coordinate_distribution_grid(count_point,8))
                allocate(env%different_tipe_of_nuclear(i)%coordinate_velocity_grid(count_point))
                env%different_tipe_of_nuclear(i)%count_point_vel_distr=count_point
                do j = 1,8
                    call get_temp(env%tem_grid(j))
                    do r = 1, count_point
                        integ = integrate_function(maxwell_speed_distribution,0.0,real(r),10000)
                        env%different_tipe_of_nuclear(i)%coordinate_distribution_grid(r,j) = integ
                        env%different_tipe_of_nuclear(i)%coordinate_velocity_grid(r) = r
                    end do
                end do 
            end do
            ! Initialize the source population at E = 100 eV with random directions.
            netrons = give_random_direction(N,100.0)
            call cpu_time(start_time)
            ! Transport each neutron independently until absorption or the 1 eV cutoff.
            env%tem = 1200
            do i = 1, N
                do while (.not. netrons(i)%is_died .and. netrons(i)%energy > 1)
                    netrons(i) = collision_controller(netrons(i), env,typeDopler)
                ! Check the position of the neutron and assign the appropriate temperature
                    localPos = sqrt(netrons(i)%pos(1)**2+netrons(i)%pos(2)**2+netrons(i)%pos(3)**2)
                    if(localPos > 2.0 .and. localPos < 4.0)then
                        local_speed_estimator(2) = local_speed_estimator(2)+netrons(i)%speed_estimator
                    else if(localPos > 4.0 .and. localPos < 6.0)then
                        local_speed_estimator(3) = local_speed_estimator(3)+netrons(i)%speed_estimator
                    else if(localPos > 6.0 .and. localPos < 8.0)then
                        local_speed_estimator(4) = local_speed_estimator(4)+netrons(i)%speed_estimator
                    else if(localPos > 8.0 .and. localPos < 10.0)then
                        local_speed_estimator(5) = local_speed_estimator(5)+netrons(i)%speed_estimator
                    else if(localPos > 10.0 .and. localPos < 12.0)then
                        local_speed_estimator(6) = local_speed_estimator(6)+netrons(i)%speed_estimator
                    else if(localPos > 12.0 .and. localPos < 15.0)then
                        local_speed_estimator(7) = local_speed_estimator(7)+netrons(i)%speed_estimator
                    else if(localPos > 15.0)then
                        local_speed_estimator(8) = local_speed_estimator(8)+netrons(i)%speed_estimator
                    else
                        local_speed_estimator(1) = local_speed_estimator(1)+netrons(i)%speed_estimator
                    end if
                !
                end do 
                age_of_netron = netrons(i)%pos(1)**2+netrons(i)%pos(2)**2+netrons(i)%pos(3)**2 + age_of_netron
                averageCollision = averageCollision + netrons(i)%count_collision
                lifeTime = lifeTime + netrons(i)%life_time
                if(.not. netrons(i)%is_died)then
                    livedDistance = sqrt(netrons(i)%pos(1)**2+netrons(i)%pos(2)**2+netrons(i)%pos(3)**2) + livedDistance
                    countLivedNetron = countLivedNetron + 1
                end if
            end do 
            call cpu_time(end_time)
            do i = 1, N
                
                localPos = sqrt(netrons(i)%pos(1)**2+netrons(i)%pos(2)**2+netrons(i)%pos(3)**2)
                if(localPos < 2.0)then
                    local_count_lived(1) = local_count_lived(1) + 1
                end if
                if(localPos < 4.0)then
                    local_count_lived(2) = local_count_lived(2) + 1
                end if
                if(localPos < 6.0)then
                    local_count_lived(3) = local_count_lived(3) + 1
                end if
                if(localPos < 8.0)then
                    local_count_lived(4) = local_count_lived(4) + 1
                end if
                if(localPos < 10.0)then
                    local_count_lived(5) = local_count_lived(5) + 1
                end if
                if(localPos < 12.0)then
                    local_count_lived(6) = local_count_lived(6) + 1
                end if
                if(localPos < 15.0)then
                    local_count_lived(7) = local_count_lived(7) + 1
                end if
                if(localPos < 100.0)then
                    local_count_lived(8) = local_count_lived(8) + 1
                end if

            end do
            livedDistance = livedDistance/countLivedNetron
            age_of_netron = age_of_netron/real(N)
            averageCollision = averageCollision/real(N)
            lifeTime = lifeTime/real(N)

            do i = 1,N 
                if(.not. netrons(i)%is_died)then
                    dispersion = (sqrt(netrons(i)%pos(1)**2+netrons(i)%pos(2)**2+netrons(i)%pos(3)**2)-livedDistance)**2+dispersion
                end if
            end do 
            dispersion = dispersion/countLivedNetron
            elapsed_time = end_time - start_time

            print *, "Растояние пройденное без поглащения: ", livedDistance
            print *, "Время работы программы: ", elapsed_time
            print *, "FOM: ", 1/(elapsed_time*dispersion)

            
            print*, "Speen of adsorption on the distans = 2 cm: ", local_speed_estimator(1)/local_count_lived(1)
            print*, "Speen of adsorption on the distans = 4 cm: ", local_speed_estimator(2)/local_count_lived(2)
            print*, "Speen of adsorption on the distans = 6 cm: ", local_speed_estimator(3)/local_count_lived(3)
            print*, "Speen of adsorption on the distans = 8 cm: ", local_speed_estimator(4)/local_count_lived(4)
            print*, "Speen of adsorption on the distans = 10 cm: ", local_speed_estimator(5)/local_count_lived(5)
            print*, "Speen of adsorption on the distans = 12 cm: ", local_speed_estimator(6)/local_count_lived(6)
            print*, "Speen of adsorption on the distans = 15 cm: ", local_speed_estimator(7)/local_count_lived(7)
            print*, "Speen of adsorption on the distans > 15 cm: ", local_speed_estimator(8)/local_count_lived(8)
    
        case(1)
        !--------------------------------------------------------------------------
        ! Mode 1: Generate Doppler-broadened cross-section tables at temperature tem.
        !--------------------------------------------------------------------------
            do i = 1,2
                file_name = "output/Cross-section-data-"//trim(env%different_tipe_of_nuclear(i)%name_of_nuclie)//"-"//int_to_str(tem)//".txt"
                open(i, file=file_name, status="replace", action="write")
                do j = 1,env%different_tipe_of_nuclear(i)%count_point
                    alpha = env%different_tipe_of_nuclear(i)%mass_of_nuclear/(k*tem)
                    write(i, '(F20.6)', advance='no') env%different_tipe_of_nuclear(i)%energy_point_in_table(j)
                    do r = 1,3
                        write(i, '(1X, F20.6)', advance='no') doplerBroadr(env%different_tipe_of_nuclear(i)%energy_point_in_table(j),&
                                                                        env%different_tipe_of_nuclear(i)%energy_point_in_table,&
                                                                        env%different_tipe_of_nuclear(i)%cross_data(r)%cross_section_point_in_table,&
                                                                        alpha,&
                                                                        size(env%different_tipe_of_nuclear(i)%energy_point_in_table))
                    end do
                    write(i,*) " "
                end do
            end do
        case(2)
        !--------------------------------------------------------------------------
        ! Mode 2: build adaptive OTF union energy grids.
        !--------------------------------------------------------------------------
            ! Build the OTF union and fitting temperature grids.
            call build_temerature_grid(conf%t_min,conf%t_max, conf%dt_union,Tgrid)
            call build_temerature_grid(conf%t_min,conf%t_max, conf%dt_fit, TfitGrid)
            ! Initialize the adaptive OTF energy grids from the original tables.
            allocate(env%different_tipe_of_nuclear(1)%e_uniq_grid(env%different_tipe_of_nuclear(1)%count_point))
            allocate(env%different_tipe_of_nuclear(2)%e_uniq_grid(env%different_tipe_of_nuclear(2)%count_point))
            env%different_tipe_of_nuclear(1)%e_uniq_grid = env%different_tipe_of_nuclear(1)%energy_point_in_table
            env%different_tipe_of_nuclear(2)%e_uniq_grid = env%different_tipe_of_nuclear(2)%energy_point_in_table
            call build_union_grid(env%different_tipe_of_nuclear(2)%e_uniq_grid,Tgrid,3,conf%ft,env%different_tipe_of_nuclear(2)%energy_point_in_table,env,1,size(env%different_tipe_of_nuclear(1)%energy_point_in_table))
            call build_union_grid(env%different_tipe_of_nuclear(1)%e_uniq_grid,Tgrid,3,conf%ft,env%different_tipe_of_nuclear(1)%energy_point_in_table,env,2,size(env%different_tipe_of_nuclear(2)%energy_point_in_table))
            print*,  size(env%different_tipe_of_nuclear(1)%e_uniq_grid)," - " , size(env%different_tipe_of_nuclear(1)%energy_point_in_table)
            print*,  size(env%different_tipe_of_nuclear(2)%e_uniq_grid)," - " , size(env%different_tipe_of_nuclear(2)%energy_point_in_table)
        case(3)
        !--------------------------------------------------------------------------
        ! Mode 3: spatial temperature-profile transport study.
        ! Temperature is updated from the neutron radial position after collisions.
        !--------------------------------------------------------------------------
            do i = 1,2
                call get_mass(env%different_tipe_of_nuclear(i)%mass_of_nuclear * 1.660539e-27)
                allocate(env%different_tipe_of_nuclear(i)%coordinate_distribution_grid(count_point,8))
                allocate(env%different_tipe_of_nuclear(i)%coordinate_velocity_grid(count_point))
                env%different_tipe_of_nuclear(i)%count_point_vel_distr=count_point
                do j = 1,8
                    call get_temp(env%tem_grid(j))
                    do r = 1, count_point
                        integ = integrate_function(maxwell_speed_distribution,0.0,real(r),10000)
                        env%different_tipe_of_nuclear(i)%coordinate_distribution_grid(r,j) = integ
                        env%different_tipe_of_nuclear(i)%coordinate_velocity_grid(r) = r
                    end do
                end do 
            end do

            netrons = give_random_direction(N,100.0)
            call cpu_time(start_time)
            do i = 1, N
                localPos = 0.0
                env%tem = 3200.0
                do while (.not. netrons(i)%is_died .and. netrons(i)%energy > 1)
                    netrons(i) = collision_controller(netrons(i), env,typeDopler)
                ! Check the position of the neutron and assign the appropriate temperature
                    localPos = sqrt(netrons(i)%pos(1)**2+netrons(i)%pos(2)**2+netrons(i)%pos(3)**2)
                    if(localPos > 2.0 .and. localPos < 4.0)then
                        env%tem = 2800.0
                        local_speed_estimator(2) = local_speed_estimator(2)+netrons(i)%speed_estimator
                    else if(localPos > 4.0 .and. localPos < 6.0)then
                        env%tem = 2200.0
                        local_speed_estimator(3) = local_speed_estimator(3)+netrons(i)%speed_estimator
                    else if(localPos > 6.0 .and. localPos < 8.0)then
                        env%tem = 1800.0
                        local_speed_estimator(4) = local_speed_estimator(4)+netrons(i)%speed_estimator
                    else if(localPos > 8.0 .and. localPos < 10.0)then
                        env%tem = 1200.0
                        local_speed_estimator(5) = local_speed_estimator(5)+netrons(i)%speed_estimator
                    else if(localPos > 10.0 .and. localPos < 12.0)then
                        env%tem = 800.0
                        local_speed_estimator(6) = local_speed_estimator(6)+netrons(i)%speed_estimator
                    else if(localPos > 12.0 .and. localPos < 15.0)then
                        env%tem = 400.0
                        local_speed_estimator(7) = local_speed_estimator(7)+netrons(i)%speed_estimator
                    else if(localPos > 15.0)then
                        env%tem = 294.0
                        local_speed_estimator(8) = local_speed_estimator(8)+netrons(i)%speed_estimator
                    else
                        env%tem = 3200.0
                        local_speed_estimator(1) = local_speed_estimator(1)+netrons(i)%speed_estimator
                    end if
                !
                end do 
                age_of_netron = netrons(i)%pos(1)**2+netrons(i)%pos(2)**2+netrons(i)%pos(3)**2 + age_of_netron
                averageCollision = averageCollision + netrons(i)%count_collision
                lifeTime = lifeTime + netrons(i)%life_time
                if(.not. netrons(i)%is_died)then
                    livedDistance = sqrt(netrons(i)%pos(1)**2+netrons(i)%pos(2)**2+netrons(i)%pos(3)**2) + livedDistance
                    countLivedNetron = countLivedNetron + 1
                end if
            end do 
            call cpu_time(end_time)
            livedDistance = livedDistance/countLivedNetron
            age_of_netron = age_of_netron/real(N)
            averageCollision = averageCollision/real(N)
            lifeTime = lifeTime/real(N)

            do i = 1,N 
                if(.not. netrons(i)%is_died)then
                    dispersion = (sqrt(netrons(i)%pos(1)**2+netrons(i)%pos(2)**2+netrons(i)%pos(3)**2)-livedDistance)**2+dispersion
                end if
            end do 
            dispersion = dispersion/countLivedNetron
            elapsed_time = end_time - start_time
            do i = 1, N
                
                localPos = sqrt(netrons(i)%pos(1)**2+netrons(i)%pos(2)**2+netrons(i)%pos(3)**2)
                if(localPos < 2.0)then
                    local_count_lived(1) = local_count_lived(1) + 1
                end if
                if(localPos < 4.0)then
                    local_count_lived(2) = local_count_lived(2) + 1
                end if
                if(localPos < 6.0)then
                    local_count_lived(3) = local_count_lived(3) + 1
                end if
                if(localPos < 8.0)then
                    local_count_lived(4) = local_count_lived(4) + 1
                end if
                if(localPos < 10.0)then
                    local_count_lived(5) = local_count_lived(5) + 1
                end if
                if(localPos < 12.0)then
                    local_count_lived(6) = local_count_lived(6) + 1
                end if
                if(localPos < 15.0)then
                    local_count_lived(7) = local_count_lived(7) + 1
                end if
                if(localPos < 100.0)then
                    local_count_lived(8) = local_count_lived(8) + 1
                end if

            end do

            print *, "Растояние пройденное без поглащения: ", livedDistance
            print *, "Время работы программы: ", elapsed_time
            print *, "FOM: ", 1/(elapsed_time*dispersion)
            print*, countLivedNetron
            
            print*, "Speen of adsorption on the distans = 2 cm: ", local_speed_estimator(1)/local_count_lived(1)
            print*, "Speen of adsorption on the distans = 4 cm: ", local_speed_estimator(2)/local_count_lived(2)
            print*, "Speen of adsorption on the distans = 6 cm: ", local_speed_estimator(3)/local_count_lived(3)
            print*, "Speen of adsorption on the distans = 8 cm: ", local_speed_estimator(4)/local_count_lived(4)
            print*, "Speen of adsorption on the distans = 10 cm: ", local_speed_estimator(5)/local_count_lived(5)
            print*, "Speen of adsorption on the distans = 12 cm: ", local_speed_estimator(6)/local_count_lived(6)
            print*, "Speen of adsorption on the distans = 15 cm: ", local_speed_estimator(7)/local_count_lived(7)
            print*, "Speen of adsorption on the distans > 15 cm: ", local_speed_estimator(8)/local_count_lived(8)
    
    end select

    contains
    !==============================================================================
    ! function: int_to_str
    !
    ! Purpose:
    !   Convert an integer to a trimmed decimal string using Fortran format I0.
    !
    ! Parametr IN:
    !   val - integer number
    !
    ! Parametr OUT:  
    !   res - string number
    !
    !==============================================================================

    function int_to_str(val) result(res)
        integer, intent(in) :: val
        character(:), allocatable :: res
        character(len=30) :: buffer  ! Буфер с запасом

        write(buffer, '(I0)') val
        res = trim(buffer)
    end function int_to_str

end program bilder  