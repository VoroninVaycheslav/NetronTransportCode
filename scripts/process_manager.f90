module process_manager

    use data_type
    use scattering_process
    use adsobtion_process
    use operation_with_data
    use math_operation
    implicit none

    contains
    ! Моделирование процесса взимодействия нейтрона с ядром
    function collision_controller(cur_netron_data, env) result(new_netron_data)
        
        type(netron_data), intent(in) :: cur_netron_data                        !Нейтрон до взаимодействия                  IN
        type(enviroment), intent(in) :: env                                     !Среда, в которой находится нейтрон         IN
        type(netron_data) :: new_netron_data                                    !Нейтрон после взаимодействия               OUT
        
        real, allocatable :: mic_current_cross_section_mas(:)                   !Микроскопические сечения, возможные при данной энергии и при данном ядре
        real, allocatable :: total_mac_current_cross_section_mas(:)             !Макроскопические сечения для всех ядер
        real :: total_mac_cross_section                                         !Полноне Макроскопические сечение среды 
        real :: total_micro_cross_section                                       !Полное микроскопическое сечение выбранного ядра
        integer :: type_of_nuclie                                               !Тип выбранного ядра
        integer :: type_of_process                                              !Тип выбранного процесса
        real :: sum_val
        real :: p1, p2
        integer :: i
        
        
        ! Выделяем память для массива макроскопического сечения
        allocate(total_mac_current_cross_section_mas(env%count_nuclear))
        
        ! Расчет макросечений для каждого типа ядра
        do i = 1, env%count_nuclear
            total_mac_current_cross_section_mas(i) = found_cross_section_from_energy(cur_netron_data%energy, &
                env%different_tipe_of_nuclear(i),1)*env%different_tipe_of_nuclear(i)%nuclear_dencity
        end do
        
        ! Выбор типа ядра
        total_mac_cross_section = sum(total_mac_current_cross_section_mas)
        
       ! Генириуем случайное число для выбора ядра 
        p1 = get_random_in_range(0.0000000000001, total_mac_cross_section)
        type_of_nuclie = 0
        sum_val = 0.0
        ! Выбираем ядро с помощью координаты на отрезке
        do while (sum_val < p1 .and. type_of_nuclie < env%count_nuclear)
            type_of_nuclie = type_of_nuclie + 1
            sum_val = sum_val + total_mac_current_cross_section_mas(type_of_nuclie)
        end do
        
        if (type_of_nuclie < 1 .or. type_of_nuclie > env%count_nuclear) then
            print *, "ОШИБКА: Неверный тип ядра после выбора: ", type_of_nuclie, p1
            new_netron_data = cur_netron_data
            deallocate(total_mac_current_cross_section_mas)
            return
        end if

        ! Выделяем память для массива микроскопического сечения
        allocate(mic_current_cross_section_mas(env%different_tipe_of_nuclear(type_of_nuclie)%count_process-1))
        
        ! Расчет микросечений для процессов (индексы: 2 - рассеяние, 3 - поглощение)
        total_micro_cross_section = 0.0
        do i = 2, env%different_tipe_of_nuclear(type_of_nuclie)%count_process
            mic_current_cross_section_mas(i-1) = found_cross_section_from_energy(cur_netron_data%energy, &
            env%different_tipe_of_nuclear(type_of_nuclie),i)
        end do
        total_micro_cross_section = sum(mic_current_cross_section_mas)
        
        p2 = get_random_in_range(0.00000000000001, total_micro_cross_section)
        type_of_process = 0
        sum_val = 0.0
        
        do while (sum_val < p2 .and. type_of_process < env%different_tipe_of_nuclear(type_of_nuclie)%count_process)
            type_of_process = type_of_process + 1
            sum_val = sum_val + mic_current_cross_section_mas(type_of_process)
        end do
        ! Выполнение выбранного процесса
        select case(type_of_process)
        case(1)
            ! Рассеяние с дыигвющимися ядрами
             new_netron_data = get_one_bump_netron_termalization(cur_netron_data, env%different_tipe_of_nuclear(type_of_nuclie)%mass_of_nuclear, total_mac_cross_section,env%different_tipe_of_nuclear(type_of_nuclie))
            
            ! Рассеяние с покоящимися ядрами
             !new_netron_data = get_one_bump_netron_slow_down(cur_netron_data, env%different_tipe_of_nuclear(type_of_nuclie)%mass_of_nuclear, total_mac_cross_section)
        case(2)
            ! Поглощение
            new_netron_data = get_absorption(cur_netron_data)
        case default
            print *, "ВНИМАНИЕ: Неизвестный тип процесса: ", type_of_process, p2
            new_netron_data = cur_netron_data
        end select
        
        ! Освобождаем память
        deallocate(mic_current_cross_section_mas, total_mac_current_cross_section_mas)
        
    end function collision_controller
    ! Создаем нейтроны со случайными направлениями движения и одинаковой энергией

    function give_random_direction(N, E) result(new_netron_d)
        integer, intent(in) :: N                        !Количество нейтронов       IN
        real, intent(in):: E
        type(netron_data) :: new_netron_d(N)            !Созданный нейтрон          OUT
        real :: T, Fi
        integer :: i

        do i = 1, N
            !Задаем случайное напрваление движения
            T = get_random_in_range(0.0, 3.14159)
            Fi = get_random_in_range(0.0, 6.28318)
            new_netron_d(i)%dir(1) = sin(T) * cos(Fi)
            new_netron_d(i)%dir(2) = sin(T) * sin(Fi)
            new_netron_d(i)%dir(3) = cos(T)
            new_netron_d(i)%energy = E
  
            !Обнуляем остальные характеристики нейтрона
            new_netron_d(i)%speed = 1.38e4 * sqrt(new_netron_d(i)%energy)
            new_netron_d(i)%pos = [0.0, 0.0, 0.0]
            new_netron_d(i)%life_time = 0.0
            new_netron_d(i)%count_collision = 0
            new_netron_d(i)%is_died = .False.
        end do
    end function give_random_direction

end module process_manager
