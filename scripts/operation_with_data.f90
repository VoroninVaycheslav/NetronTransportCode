module operation_with_data

    use data_type
    implicit none

    
    contains
    !Поиск сечения по заданной энергии
    function found_cross_section_from_energy(current_energy, tem, env, nuclear_d, index_process) result(current_section)
        type(nuclear_data), intent(in) :: nuclear_d         !Информация о текущем ядре (его сечения)        IN
        real, intent(in) :: current_energy                  !Текущая энергия                                IN
        real, intent(in) :: tem                             !Индекс температуры                             IN
        type(enviroment), intent(in) :: env                 !Исследуемая среды                              IN
        integer, intent(in) :: index_process                !Индекс процесса                                IN
        real :: current_section                             !Какое сечение имеет нейтрон сейчас             OUT
        real :: l_lim_cs,r_lim_cs,l_lim_e,r_lim_e           !Левые и правые границы по энергиям и сечению
        integer :: left, right, mid, tem_index              
        
        do tem_index = 1, size(env%tem_grid)
            if (tem == env%tem_grid(tem_index)) exit
        end do

        left = 1
        right = nuclear_d%count_point
        
        ! Бинарный поиск
        do while (right - left > 1)
            mid = (left + right) / 2
            if (current_energy <= nuclear_d%energy_point_in_table(mid)) then
                right = mid
            else
                left = mid
            end if
        end do

        l_lim_cs = nuclear_d%cross_data(index_process)%cross_section_point_in_table(left,tem_index)
        r_lim_cs = nuclear_d%cross_data(index_process)%cross_section_point_in_table(right,tem_index)
        l_lim_e = nuclear_d%energy_point_in_table(left)
        r_lim_e = nuclear_d%energy_point_in_table(right)
        ! Вычисление сечения на отрезке
        current_section = l_lim_cs + (r_lim_cs - l_lim_cs) * (current_energy - l_lim_e) / (r_lim_e - l_lim_e)

    end function found_cross_section_from_energy

    !Загружаем все данные о ядрах для серды из файлов
    subroutine load_cross_section_fron_file(size_data_list,count_tem,directory_of_cross_section_data, env)
        integer, intent(in)::size_data_list                                                     !Количество считываемых файлов      IN
        integer, intent(in) :: count_tem                                                        !Количество температур              IN
        character(len=100),intent(in) :: directory_of_cross_section_data(size_data_list,count_tem)        !Сами файлы (путь)                  IN
        type(enviroment), intent(inout)::env                                                    !Исследуемая среды                  IN-OUT
        
        integer err, i, unit_num, j, k, m, trust  

        !Устанавливаем количество ядер в среде
        env%count_nuclear = size_data_list

        !Выделяем место под ядра в среде с проверкой
        allocate(env%different_tipe_of_nuclear(size_data_list), stat=err)
        if (err /= 0) then
            print *, 'Ошибка выделения памяти! Код ошибки:', err
            stop
        else
            print *, 'Добавлено ', size_data_list,'ядер в среду'
        end if

        allocate(env%tem_grid(count_tem), stat=err)
        !Обработка файлов 
        do i = 1 ,size_data_list

            !Открываем файл и проверяем открытие
            open(1, file=directory_of_cross_section_data(i,1), status="old", action="read", iostat=err)
            if (err /= 0) then
                print *, "ОШИБКА: Не удалось открыть файл: ", directory_of_cross_section_data(i,1)
                stop
            end if

            !считываем все данные о ядрах кроме сечений
            read(1, *, iostat=err) env%different_tipe_of_nuclear(i)%name_of_nuclie
            read(1, *, iostat=err) env%different_tipe_of_nuclear(i)%index_of_nuclie
            read(1, *, iostat=err) env%different_tipe_of_nuclear(i)%mass_of_nuclear
            read(1, *, iostat=err) env%different_tipe_of_nuclear(i)%nuclear_dencity
            read(1, *, iostat=err) env%different_tipe_of_nuclear(i)%count_process
            read(1, *, iostat=err) env%different_tipe_of_nuclear(i)%count_point

            !выделяем пмять под столбец энергий и проверяем успех выделения
            allocate(env%different_tipe_of_nuclear(i)%energy_point_in_table(env%different_tipe_of_nuclear(i)%count_point), stat=err)
            if (err /= 0) then
                print *, 'Ошибка выделения памяти! Код ошибки:', err
                stop
            else
                print *, 'Добавлено ', env%different_tipe_of_nuclear(i)%count_point,'точек в таблицу, для ядра ',env%different_tipe_of_nuclear(i)%name_of_nuclie
            end if

            !Выделяем память для количества столбцов по сечениям и проверяем успех
            allocate(env%different_tipe_of_nuclear(i)%cross_data(env%different_tipe_of_nuclear(i)%count_process), stat=err)
            if (err /= 0) then
                print *, 'Ошибка выделения памяти! Код ошибки:', err
                stop
            else
                print *, 'Добавлено ', env%different_tipe_of_nuclear(i)%count_process,'типов процессов для ядра ',env%different_tipe_of_nuclear(i)%name_of_nuclie
            end if

            !Считываем индекс процесса
            read(1, *, iostat=err) env%different_tipe_of_nuclear(i)%cross_data%index_of_process

            !выделяем место под каждый столбец с нужным количеством точек для сечений
            do j = 1, env%different_tipe_of_nuclear(i)%count_process
                allocate(env%different_tipe_of_nuclear(i)%cross_data(j)%cross_section_point_in_table(env%different_tipe_of_nuclear(i)%count_point,count_tem), stat=err)
                if (err /= 0) then
                    print *, 'Ошибка выделения памяти! Код ошибки:', err
                    stop
                else
                    print *, 'Добавлено ', env%different_tipe_of_nuclear(i)%count_point,'точек для ядра ',env%different_tipe_of_nuclear(i)%name_of_nuclie, 'для процесса с индексом ',env%different_tipe_of_nuclear(i)%cross_data(j)%index_of_process
                end if
            end do

            !Заполняем столбцы с сечениями
            do j = 1, env%different_tipe_of_nuclear(i)%count_point
                read(1, *, iostat=err) env%different_tipe_of_nuclear(i)%energy_point_in_table(j),(env%different_tipe_of_nuclear(i)%cross_data(k)%cross_section_point_in_table(j,1), k=1, env%different_tipe_of_nuclear(i)%count_process)
            end do

            do m = 2, count_tem
                print *, "Загружаем данные о ядре ", m
                open(newunit=unit_num, file=directory_of_cross_section_data(i,m), status="old", action="read", iostat=err)
                if (err /= 0) then
                    print *, "ОШИБКА: Не удалось открыть файл: ", directory_of_cross_section_data(i,m)
                    stop
                end if
                do j = 1, env%different_tipe_of_nuclear(i)%count_point
                    read(unit_num, *, iostat=err) env%different_tipe_of_nuclear(i)%energy_point_in_table(j), (env%different_tipe_of_nuclear(i)%cross_data(k)%cross_section_point_in_table(j,m), k=1, env%different_tipe_of_nuclear(i)%count_process)
                end do
                close(unit_num)
            end do
        end do

    end subroutine load_cross_section_fron_file

end module operation_with_data
