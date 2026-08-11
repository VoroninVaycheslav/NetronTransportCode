module math_operation
    
    implicit none

    integer :: T = 300
    real :: m = 12.0 * 1.660539e-27
    contains
    !Генерация случайного числа в диапазоне
    function get_random_in_range(xmin, xmax) result(res)
        real, intent(in) :: xmin, xmax
        real :: res
        real :: t
        
        call random_number(t)
        res = xmin + (xmax - xmin) * t
    end function get_random_in_range
  
    function integrate_function(f, left_lim, right_lim, n) result(res)
        real, external :: f                             ! Интегрируемая функция
        real, intent(in) :: left_lim, right_lim         ! Пределы интегрирования
        integer, intent(in) :: n                        ! Точность разбиения (количество точек)
        real :: res                                     ! результат интегрирования
        real :: h, x_mid, area
        integer :: i
        
        h = (right_lim - left_lim) / n
        res = 0.0
        
        do i = 0, n-1
            x_mid = left_lim + (i + 0.5) * h   ! середина столбика
            area = f(x_mid) * h
            res = res + area
        end do
    end function integrate_function

    ! Функция распределения максвела
    function maxwell_distribution(E) result(res)
        real, intent(in) :: E
        real :: pi = 3.14159265
        real :: T_M = 1.333
        real :: res

        res = 2.0 * sqrt(E / (pi * T_M**3)) / exp(E / T_M)
    end function maxwell_distribution
    function maxwell_speed_distribution(v) result(res)

        real, intent(in) :: v                           ! Скорсоть ядра
        real :: res                                     ! Итог плотность вероятности
        real :: factor                                  
        real, parameter :: pi = 3.141592653589793       ! Число пи
        real, parameter :: k_B = 1.380649e-23           ! Постоянная больцмана
        if (v < 0.0) then
            res = 0.0
            return
        end if
        
        ! f(v) = 4π (m/(2πkT))^(3/2) v^2 exp(-mv^2/(2kT))
        factor = (m / (2.0 * pi * k_B * T)) ** 1.5
        res = 4.0 * pi * factor * v**2 * exp(-m * v**2 / (2.0 * k_B * T))
    
    end function maxwell_speed_distribution

    subroutine get_mass(mass_env)
        real, intent(in) :: mass_env
        m = mass_env
    end subroutine get_mass
    function found_number_with_table(num, table_search, table_res,count_point) result(res)
        real, intent(in) :: num                                             ! Вероятность подаваяемяа для поиска скорости
        integer, intent(in) :: count_point                                  ! Количество точек в сетке
        real, dimension(count_point), intent(in) :: table_search            ! Столбец скорсотей
        real, dimension(count_point), intent(in) :: table_res               ! Столбец вероятностей     
        real :: res     

        real :: l_lim_v,r_lim_v,l_lim_d,r_lim_d           
        integer :: left, right, mid                                               ! Найденная скорость
        left = 1
        right = count_point
        
        ! Бинарный поиск
        do while (right - left > 1)
            mid = (left + right) / 2
            if (num <= table_search(mid)) then
                right = mid
            else
                left = mid
            end if
        end do

        ! Линейная интерполяция
        l_lim_v = table_res(left)
        r_lim_v = table_res(right)
        l_lim_d = table_search(left)
        r_lim_d = table_search(right)
        res = l_lim_v + (r_lim_v - l_lim_v) * (num - l_lim_d) / (r_lim_d - l_lim_d)
            
    end function found_number_with_table
end module math_operation
