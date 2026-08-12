program bilder
    use data_type
    use math_operation
    use operation_with_data
    use process_manager
    use DoplerTrearmentNJOY
    implicit none

    integer :: N = 100000
    integer i,r,j
    integer :: tem = 1200
    real integ,alpha
    integer :: typeWork = 0  ! 0 - моделирование, 1 - генерация таблиц сечений
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
    real :: k = 8.617333262e-5

    character(len=100) :: directory_of_cross_section_data(2)
    character(len=100) :: file_name

    allocate(netrons(N))

    directory_of_cross_section_data(1) = "dataBaseOfCrossSection/Cross-section-data-Carbon-294.txt"
    directory_of_cross_section_data(2) = "dataBaseOfCrossSection/Cross-section-data-Uranium-294.txt"
    
    call load_cross_section_fron_file(2,directory_of_cross_section_data,env)
    !do i = 1,size(env%different_tipe_of_nuclear(1)%energy_point_in_table)
    !    print*,env%different_tipe_of_nuclear(1)%energy_point_in_table(i),found_cross_section_from_energy(env%different_tipe_of_nuclear(1)%energy_point_in_table(i),env%different_tipe_of_nuclear(1),1)
    !end do
    select case(typeWork)
        case(0)
            do i = 1,2
        ! Выделяем место под массивы
                call get_mass(env%different_tipe_of_nuclear(i)%mass_of_nuclear * 1.660539e-27)
                allocate(env%different_tipe_of_nuclear(i)%coordinate_distribution_grid(count_point))
                allocate(env%different_tipe_of_nuclear(i)%coordinate_velocity_grid(count_point))
                env%different_tipe_of_nuclear(i)%count_point_vel_distr=count_point
                do r = 0, count_point
                    integ = integrate_function(maxwell_speed_distribution,0.0,real(r),10000)
                    env%different_tipe_of_nuclear(i)%coordinate_distribution_grid(r) = integ
                    env%different_tipe_of_nuclear(i)%coordinate_velocity_grid(r) = r
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