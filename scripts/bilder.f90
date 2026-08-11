program bilder
    use data_type
    use math_operation
    use operation_with_data
    use process_manager

    implicit none
    integer :: N = 100000
    integer i
    real :: age_of_netron = 0
    integer :: averageCollision = 0
    real :: lifeTime = 0
    type(netron_data),allocatable :: netrons(:)
    type(enviroment) env

    character(len=100) :: directory_of_cross_section_data(1)

    allocate(netrons(N))

    directory_of_cross_section_data = "dataBaseOfCrossSection/Cross-section-data-Carbon.txt"
    netrons = give_random_direction(N,2e6)
    
    call load_cross_section_fron_file(1,directory_of_cross_section_data,env)

    do i = 1, N,1
        do while (netrons(i)%energy > 1)
            netrons(i) = collision_controller(netrons(i), env)
        end do 
        age_of_netron = netrons(i)%pos(1)**2+netrons(i)%pos(2)**2+netrons(i)%pos(3)**2 + age_of_netron
        averageCollision = averageCollision + netrons(i)%count_collision
        lifeTime = lifeTime + netrons(i)%life_time
    end do 
    age_of_netron = age_of_netron/real(N)
    averageCollision = averageCollision/real(N)
    lifeTime = lifeTime/real(N)

    print*, "Возраст: ", age_of_netron/6
    print*, "Среднее количество столкновений: ", averageCollision
    print*, "Время замедления: ", lifeTime
    

end program bilder  