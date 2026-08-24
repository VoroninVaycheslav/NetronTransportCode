!==============================================================================
! Module: operation_with_data
!
! Purpose:
!   Loads tabulated cross sections and OTF coefficient data, stores them in the
!   shared nuclear-data structures, and evaluates tabulated cross sections at a
!   requested neutron energy and temperature-grid index.
!
!==============================================================================
module operation_with_data

    use data_type
    implicit none

    
    contains

    !==============================================================================
    ! function: found_cross_section_from_energy
    !
    ! Purpose:
    !   Interpolate a microscopic cross section at a requested energy.
    !
    ! Parametr IN:
    !   current_energy            -   Neutron energy [eV].
    !   tem                       -   Requested material temperature [K].
    !   env                       -   Environment containing the discrete temperat
    !   nuclear_d                 -   Nuclide data and cross-section tables.
    !   index_process             -   Internal reaction-column index.
    !
    ! Parametr OUT:  
    !   current_section         -   Linearly interpolated microscopic cross section.
    !
    !==============================================================================

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

    !==============================================================================
    ! subroutine: load_cross_section_fron_file
    !
    ! Purpose:
    !   Load multi-temperature nuclide cross-section tables from text files.
    !
    ! Parametr IN:
    !   size_data_list                          -   Number of nuclides/files in the first array dimension.
    !   count_tem                               -   Number of temperatures supplied for each nuclide.
    !   directory_of_cross_section_data         -   Matrix of input file paths.
    !
    ! Parametr IN/OUT:  
    !   env                                     -   Environment populated with nuclide metadata and cross sections.
    !
    !==============================================================================

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

    !==============================================================================
    ! subroutine: load_OTF_coefficients
    !
    ! Purpose:
    !   Load precomputed OTF expansion coefficients for one nuclide.
    !
    ! Parametr IN:
    !   filenames    -   One coefficient file per reaction.
    !
    ! Parametr IN/OUT:  
    !   nuclear_d    -   Nuclide object receiving K_OTF
    !
    !==============================================================================

    subroutine load_OTF_coefficients(filenames, nuclear_d)

        character(len=*), intent(in) :: filenames(:)
        type(nuclear_data), intent(inout) :: nuclear_d

        integer :: unit_num, ios
        integer :: n_rows, n_reactions
        integer :: i, j, r
        integer :: point_index

        real :: energy_dummy
        real :: error_dummy
        real :: coeff(17)

        character(len=2048) :: line


        !============================================================
        ! Количество реакций = количество файлов
        !============================================================

        n_reactions = size(filenames)

        if (n_reactions <= 0) then
            error stop "No OTF coefficient files"
        end if


        !============================================================
        ! Первый файл:
        ! определяем количество энергетических точек
        !============================================================

        open(newunit=unit_num, &
            file=trim(filenames(1)), &
            status="old", &
            action="read", &
            iostat=ios)

        if (ios /= 0) then
            print *, "Cannot open file: ", trim(filenames(1))
            error stop
        end if


        n_rows = 0

        do

            read(unit_num, '(A)', iostat=ios) line

            if (ios /= 0) exit

            if (len_trim(line) == 0) cycle

            n_rows = n_rows + 1

        end do

        close(unit_num)


        if (n_rows <= 0) then
            error stop "OTF coefficient file is empty"
        end if


        !============================================================
        ! Выделяем трехмерный массив:
        !
        ! K_OTF(energy, reaction, coefficient)
        !
        ! 1 dimension -> energy point
        ! 2 dimension -> reaction
        ! 3 dimension -> C1...C17
        !============================================================

        if (allocated(nuclear_d%K_OTF)) then
            deallocate(nuclear_d%K_OTF)
        end if


        allocate(nuclear_d%K_OTF(n_rows, n_reactions, 17))

        nuclear_d%K_OTF = 0.0


        !============================================================
        ! Читаем ВСЕ файлы
        !============================================================

        do r = 1, n_reactions


            open(newunit=unit_num, &
                file=trim(filenames(r)), &
                status="old", &
                action="read", &
                iostat=ios)


            if (ios /= 0) then

                print *, "Cannot open file: ", trim(filenames(r))
                error stop

            end if


            !--------------------------------------------------------
            ! Формат строки:
            !
            ! index energy C1 C2 ... C17 error
            !
            ! Сохраняем только C1...C17
            !--------------------------------------------------------

            do i = 1, n_rows


                read(unit_num, *, iostat=ios) &
                    point_index,              &
                    energy_dummy,             &
                    (coeff(j), j = 1,17),     &
                    error_dummy


                if (ios /= 0) then

                    print *, "Error reading file:"
                    print *, trim(filenames(r))
                    print *, "Reaction =", r
                    print *, "Line     =", i

                    close(unit_num)

                    error stop

                end if


                !====================================================
                ! energy, reaction, coefficient
                !====================================================

                nuclear_d%K_OTF(i, r, 1:17) = coeff(1:17)


            end do


            !--------------------------------------------------------
            ! Проверяем, что в файле нет дополнительных строк
            !--------------------------------------------------------

            do

                read(unit_num, '(A)', iostat=ios) line

                if (ios /= 0) exit

                if (len_trim(line) /= 0) then

                    print *, "Different number of energy points:"
                    print *, trim(filenames(r))

                    close(unit_num)

                    error stop

                end if

            end do


            close(unit_num)


            print *, "Loaded reaction ", r, ": ", trim(filenames(r))


        end do


        print *, "========================================"
        print *, "OTF coefficients loaded"
        print *, "Energy points :", size(nuclear_d%K_OTF,1)
        print *, "Reactions     :", size(nuclear_d%K_OTF,2)
        print *, "Coefficients  :", size(nuclear_d%K_OTF,3)
        print *, "========================================"


    end subroutine load_OTF_coefficients

    !==============================================================================
    ! subroutine: load_unique_energy_grid
    !
    ! Purpose:
    !   Load the unique energy grid associated with OTF coefficients.
    !
    ! Parametr IN:
    !   filename      -   Text file containing one energy value per row.
    !
    ! Parametr IN/OUT:  
    !   nuclear_d     -   Nuclide object receiving K_OTF
    !
    !==============================================================================

    subroutine load_unique_energy_grid(filename, nuclear_d)

        character(len=*), intent(in) :: filename
        type(nuclear_data), intent(inout) :: nuclear_d

        integer :: unit_num
        integer :: ios
        integer :: n_rows
        integer :: i
        real :: energy_value

        ! ============================================================
        ! Первый проход:
        ! считаем количество энергетических точек
        ! ============================================================

        open(newunit=unit_num, &
            file=trim(filename), &
            status="old", &
            action="read", &
            iostat=ios)

        if (ios /= 0) then
            print *, "Cannot open unique energy grid file:"
            print *, trim(filename)
            error stop
        end if

        n_rows = 0

        do
            read(unit_num, *, iostat=ios) energy_value

            if (ios < 0) exit

            if (ios > 0) then
                print *, "Error while reading energy grid"
                error stop
            end if

            n_rows = n_rows + 1
        end do

        close(unit_num)


        ! ============================================================
        ! Выделяем массив e_uniq_grid
        ! ============================================================

        if (allocated(nuclear_d%e_uniq_grid)) then
            deallocate(nuclear_d%e_uniq_grid)
        end if

        allocate(nuclear_d%e_uniq_grid(n_rows))

        
        ! ============================================================
        ! Второй проход:
        ! загружаем энергии
        ! ============================================================

        open(newunit=unit_num, &
            file=trim(filename), &
            status="old", &
            action="read", &
            iostat=ios)

        if (ios /= 0) then
            print *, "Cannot open unique energy grid file:"
            print *, trim(filename)
            error stop
        end if

        do i = 1, n_rows

            read(unit_num, *, iostat=ios) nuclear_d%e_uniq_grid(i)

            if (ios /= 0) then
                print *, "Error reading energy point:", i
                error stop
            end if

        end do

        close(unit_num)


        ! ============================================================
        ! Информация
        ! ============================================================

        print *, "Unique energy grid loaded"
        print *, "Number of energy points:", size(nuclear_d%e_uniq_grid)

    end subroutine load_unique_energy_grid
end module operation_with_data
