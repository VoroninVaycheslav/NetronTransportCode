program bilder
    use data_type
    use math_operation
    use operation_with_data
    use process_manager
    use DoplerTrearmentNJOY
    use MCNPconfigOTF
    use Dopler_MCNP_OTF

    implicit none

    integer :: N = 100000
    integer i,r,j
    integer :: tem = 3200
    real integ,alpha
    integer :: typeWork = 3 ! 0 - моделирование, 1 - генерация таблиц сечений, 2 - генерация таблиц сечений для OTF
    real :: start_time, end_time, elapsed_time
    real :: dispersion = 0
    real :: age_of_netron = 0
    real :: livedDistance = 0
    integer :: averageCollision = 0
    integer :: countLivedNetron = 0
    real :: lifeTime = 0
    integer :: count_point = 500  
    type(netron_data),allocatable :: netrons(:)
    type(enviroment) env
    type(otf_config_type) conf

    real, allocatable :: e_uniq_grid(:,:),Tgrid(:),TfitGrid(:)
    character(len=100) :: directory_of_cross_section_data(2,8)
    character(len=100) :: directory_of_coeff_data_Carbon(3)
    character(len=100) :: directory_of_coeff_data_Uranium(3)
    character(len=100) :: file_name
    real::  local_speed_estimator(8) = 0
    real:: local_count_lived(8) = 0
    real :: localPos = 0.0

    allocate(netrons(N))
    !Reader
        directory_of_cross_section_data(1,1) = "dataBaseOfCrossSection/Cross-section-data-Carbon-294.txt"
        directory_of_cross_section_data(1,2) = "dataBaseOfCrossSection/Cross-section-data-Carbon-400.txt"
        directory_of_cross_section_data(1,3) = "dataBaseOfCrossSection/Cross-section-data-Carbon-800.txt"
        directory_of_cross_section_data(1,4) = "dataBaseOfCrossSection/Cross-section-data-Carbon-1200.txt"
        directory_of_cross_section_data(1,5) = "dataBaseOfCrossSection/Cross-section-data-Carbon-1800.txt"
        directory_of_cross_section_data(1,6) = "dataBaseOfCrossSection/Cross-section-data-Carbon-2200.txt"
        directory_of_cross_section_data(1,7) = "dataBaseOfCrossSection/Cross-section-data-Carbon-2800.txt"
        directory_of_cross_section_data(1,8) = "dataBaseOfCrossSection/Cross-section-data-Carbon-3200.txt"

        directory_of_cross_section_data(2,1) = "dataBaseOfCrossSection/Cross-section-data-Uranium-294.txt"
        directory_of_cross_section_data(2,2) = "dataBaseOfCrossSection/Cross-section-data-Uranium-400.txt"
        directory_of_cross_section_data(2,3) = "dataBaseOfCrossSection/Cross-section-data-Uranium-800.txt"
        directory_of_cross_section_data(2,4) = "dataBaseOfCrossSection/Cross-section-data-Uranium-1200.txt"
        directory_of_cross_section_data(2,5) = "dataBaseOfCrossSection/Cross-section-data-Uranium-1800.txt"
        directory_of_cross_section_data(2,6) = "dataBaseOfCrossSection/Cross-section-data-Uranium-2200.txt"
        directory_of_cross_section_data(2,7) = "dataBaseOfCrossSection/Cross-section-data-Uranium-2800.txt"
        directory_of_cross_section_data(2,8) = "dataBaseOfCrossSection/Cross-section-data-Uranium-3200.txt"

        directory_of_coeff_data_Carbon(1) = "CoefData/coefficients-MT-1-Carbon.txt"
        directory_of_coeff_data_Carbon(2) = "CoefData/coefficients-MT-2-Carbon.txt"
        directory_of_coeff_data_Carbon(3) = "CoefData/coefficients-MT-102-Carbon.txt"
        
        directory_of_coeff_data_Uranium(1) = "CoefData/coefficients-MT-1-Uranium.txt"
        directory_of_coeff_data_Uranium(2) = "CoefData/coefficients-MT-2-Uranium.txt"
        directory_of_coeff_data_Uranium(3) = "CoefData/coefficients-MT-102-Uranium.txt"
    !
    call load_cross_section_fron_file(2,8,directory_of_cross_section_data,env)
    
    call load_OTF_coefficients(directory_of_coeff_data_Carbon,env%different_tipe_of_nuclear(1))
    call load_OTF_coefficients(directory_of_coeff_data_Uranium,env%different_tipe_of_nuclear(2))

    call load_unique_energy_grid("CoefData/unique-energy-grid-Carbon.txt",env%different_tipe_of_nuclear(1))
    call load_unique_energy_grid("CoefData/unique-energy-grid-Uranium.txt",env%different_tipe_of_nuclear(2))
    env%tem_grid = [294.0,400.0,800.0,1200.0,1800.0,2200.0,2800.0, 3200.0]
    
    select case(typeWork)
        case(0)
            do i = 1,2
            ! Выделяем место под массивы
                call get_mass(env%different_tipe_of_nuclear(i)%mass_of_nuclear * 1.660539e-27)
                allocate(env%different_tipe_of_nuclear(i)%coordinate_distribution_grid(count_point,8))
                allocate(env%different_tipe_of_nuclear(i)%coordinate_velocity_grid(count_point))
                env%different_tipe_of_nuclear(i)%count_point_vel_distr=count_point
                do j = 1,8
                    call get_temp(env%tem_grid(j))
                    do r = 0, count_point
                        integ = integrate_function(maxwell_speed_distribution,0.0,real(r),10000)
                        env%different_tipe_of_nuclear(i)%coordinate_distribution_grid(r,j) = integ
                        env%different_tipe_of_nuclear(i)%coordinate_velocity_grid(r) = r
                    end do
                end do 
            end do

            netrons = give_random_direction(N,100.0)
            call cpu_time(start_time)
            do i = 1, N
                do while (.not. netrons(i)%is_died .and. netrons(i)%energy > 1)
                    netrons(i) = collision_controller(netrons(i), env)

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

            print *, "Растояние пройденное без поглащения: ", livedDistance
            print *, "Время работы программы: ", elapsed_time
            print *, "FOM: ", 1/(elapsed_time*dispersion)
        case(1)
            do i = 1,2

                file_name = "output/Cross-section-data-"//trim(env%different_tipe_of_nuclear(i)%name_of_nuclie)//"-"//int_to_str(tem)//".txt"
                open(i, file=file_name, status="replace", action="write")
                ! === ЗАГОЛОВОК (с переводом строки) ===
                !!write(i, '(A)') trim(env%different_tipe_of_nuclear(i)%name_of_nuclie)
                !write(i, '(I5)') env%different_tipe_of_nuclear(i)%index_of_nuclie
                !write(i, '(F10.6)') env%different_tipe_of_nuclear(i)%mass_of_nuclear
                !write(i, '(F10.6)') env%different_tipe_of_nuclear(i)%nuclear_dencity
                !write(i, '(I5)') env%different_tipe_of_nuclear(i)%count_process
                !write(i, '(I5)') env%different_tipe_of_nuclear(i)%count_point
                
                ! === ИНДЕКСЫ ПРОЦЕССОВ (в одну строку через пробел) ===
                !do j = 1, env%different_tipe_of_nuclear(i)%count_process
                !    write(i, '(I5, 1X)', advance='no') env%different_tipe_of_nuclear(i)%cross_data(j)%index_of_process
                !end do
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
            ! Делаем темпиратурыне сетки
            call build_temerature_grid(conf%t_min,conf%t_max, conf%dt_union,Tgrid)
            call build_temerature_grid(conf%t_min,conf%t_max, conf%dt_fit, TfitGrid)
            ! Делаем уникальную сетку по энергиям
            allocate(env%different_tipe_of_nuclear(1)%e_uniq_grid(env%different_tipe_of_nuclear(1)%count_point))
            allocate(env%different_tipe_of_nuclear(2)%e_uniq_grid(env%different_tipe_of_nuclear(2)%count_point))
            env%different_tipe_of_nuclear(1)%e_uniq_grid = env%different_tipe_of_nuclear(1)%energy_point_in_table
            env%different_tipe_of_nuclear(2)%e_uniq_grid = env%different_tipe_of_nuclear(2)%energy_point_in_table
            call build_union_grid(env%different_tipe_of_nuclear(2)%e_uniq_grid,Tgrid,3,conf%ft,env%different_tipe_of_nuclear(2)%energy_point_in_table,env,1,size(env%different_tipe_of_nuclear(1)%energy_point_in_table))
            call build_union_grid(env%different_tipe_of_nuclear(1)%e_uniq_grid,Tgrid,3,conf%ft,env%different_tipe_of_nuclear(1)%energy_point_in_table,env,2,size(env%different_tipe_of_nuclear(2)%energy_point_in_table))
            print*,  size(env%different_tipe_of_nuclear(1)%e_uniq_grid)," - " , size(env%different_tipe_of_nuclear(1)%energy_point_in_table)
            print*,  size(env%different_tipe_of_nuclear(2)%e_uniq_grid)," - " , size(env%different_tipe_of_nuclear(2)%energy_point_in_table)
         case(3)
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

            netrons = give_random_direction(N,100.0)
            call cpu_time(start_time)
            do i = 1, N
                localPos = 0.0
                env%tem = 3200.0
                do while (.not. netrons(i)%is_died .and. netrons(i)%energy > 1)
                    netrons(i) = collision_controller(netrons(i), env)

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

    function int_to_str(val) result(res)
        integer, intent(in) :: val
        character(:), allocatable :: res
        character(len=30) :: buffer  ! Буфер с запасом

        ! Формат I0 автоматически подбирает нужную ширину без лишних пробелов
        write(buffer, '(I0)') val
        res = trim(buffer)
    end function int_to_str

end program bilder  