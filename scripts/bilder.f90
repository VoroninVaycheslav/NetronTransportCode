program bilder
    use data_type
    use math_operation
    use operation_with_data
    use process_manager

    implicit none
    integer :: N = 100000
    integer i
    real :: start_time, end_time, elapsed_time
    real :: dispersion = 0
    real :: age_of_netron = 0
    real :: livedDistance = 0
    integer :: averageCollision = 0
    integer :: countLivedNetron = 0
    real :: lifeTime = 0
    type(netron_data),allocatable :: netrons(:)
    type(enviroment) env

    character(len=100) :: directory_of_cross_section_data(2)

    allocate(netrons(N))

    directory_of_cross_section_data(1) = "dataBaseOfCrossSection/Cross-section-data-Carbon.txt"
    directory_of_cross_section_data(2) = "dataBaseOfCrossSection/Cross-section-data-Uranium.txt"

    
    call load_cross_section_fron_file(2,directory_of_cross_section_data,env)

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

end program bilder  