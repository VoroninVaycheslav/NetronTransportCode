!==============================================================================
! Module: math_operation
!
! Purpose:
!   Provides random sampling, midpoint numerical integration, Maxwell
!   distributions, and table lookup/interpolation utilities used by the
!   transport physics modules.
!
!==============================================================================

module math_operation
    use data_type
    implicit none

    integer :: T = 300                  ! Module-level Maxwell temperature [K].
    real :: m = 12.0 * 1.660539e-27     ! Module-level target mass [kg].
    contains
    
    !==============================================================================
    ! function: get_random_in_range
    !
    ! Purpose:
    !   Draw a uniformly distributed pseudo-random real value.
    !
    ! Parametr IN:
    !   xmin      -   Lower bound of the sampling interval.
    !   xmax      -   Upper bound of the sampling interval.
    !
    ! Parametr OUT:  
    !   res       -   Random value xmin + (xmax-xmin)*U, with U supplied by random_number.
    !
    !==============================================================================

    function get_random_in_range(xmin, xmax) result(res)
        real, intent(in) :: xmin, xmax
        real :: res
        real :: t
        
        call random_number(t)
        res = xmin + (xmax - xmin) * t
    end function get_random_in_range
  
    !==============================================================================
    ! function: integrate_function
    !
    ! Purpose:
    !   Integrate a scalar function by the midpoint rectangle rule.
    !
    ! Parametr IN:
    !   f           -   Scalar external function to integrate.
    !   left_lim    -   Lower integration limit.
    !   right_lim   -   Upper integration limit.
    !   n           -   Number of equal subintervals used by the midpoint rule.
    !
    ! Parametr OUT:  
    !   res         -  Approximate integral over [left_lim, right_lim].
    !
    !==============================================================================

    function integrate_function(f, left_lim, right_lim, n) result(res)
        real, external :: f                             
        real, intent(in) :: left_lim, right_lim         
        integer, intent(in) :: n                        
        real :: res               ! Numerical integral result.                                  
        real :: h, x_mid, area    ! Step size, subinterval midpoint, and local area.        
        integer :: i              ! Midpoint-rule loop index.
        
        h = (right_lim - left_lim) / n
        res = 0.0
        
        do i = 0, n-1
            x_mid = left_lim + (i + 0.5) * h   ! Midpoint of the current integration subinterval.
            area = f(x_mid) * h
            res = res + area
        end do
    end function integrate_function

    !==============================================================================
    ! function: maxwell_distribution
    !
    ! Purpose:
    !   Evaluate the legacy energy-form Maxwell distribution.
    !
    ! Parametr IN:
    !   E   -  Energy-like independent variable used by the legacy model.
    !
    ! Parametr OUT:  
    !   res -  Distribution value evaluated with the fixed T_M parameter.
    !
    !==============================================================================

    function maxwell_distribution(E) result(res)
        real, intent(in) :: E
        real :: pi = 3.14159265
        real :: T_M = 1.333
        real :: res

        res = 2.0 * sqrt(E / (pi * T_M**3)) / exp(E / T_M)
    end function maxwell_distribution
    
    !==============================================================================
    ! function: maxwell_speed_distribution
    !
    ! Purpose:
    !   Evaluate the Maxwell speed probability-density function.
    !
    ! Parametr IN:
    !   v   -  Energy-like independent variable used by the legacy model.
    !
    ! Parametr OUT:  
    !   res -  Distribution value evaluated with the fixed T_M parameter.
    !
    !==============================================================================

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

    !==============================================================================
    ! subroutine: get_mass
    !
    ! Purpose:
    !   Set the target mass used by maxwell_speed_distribution.
    !
    !==============================================================================

    subroutine get_mass(mass_env)
        real, intent(in) :: mass_env
        m = mass_env
    end subroutine get_mass

    !==============================================================================
    ! subroutine: get_temp
    !
    ! Purpose:
    !   Set the temperature used by maxwell_speed_distribution.
    !
    !==============================================================================

    subroutine get_temp(tem_env)
        real, intent(in) :: tem_env
        T = tem_env
    end subroutine get_temp

!==============================================================================
    ! function: found_number_with_table
    !
    ! Purpose:
    !   Interpolate an inverse tabulated relation at the requested temperature.
    !
    ! Parametr IN:
    !   num             -   Value searched on the table_search axis.
    !   table_search    -   Search-coordinate table, indexed by point and temperatu
    !   table_res       -   Result values associated with the search coordinate.
    !   count_point     -   Number of tabulated points.
    !   tem             -   Requested temperature [K]; expected to match env%tem_grid exactl
    !   env             -   Material environment containing the temperature grid.
    !
    ! Parametr OUT:  
    !   res             -   Linearly interpolated table_res value with endpoint clamping.
    !
    !==============================================================================

    function found_number_with_table(num, table_search, table_res, count_point,tem,env) result(res)
        real, intent(in) :: num
        type(enviroment), intent(in) :: env
        integer, intent(in) :: count_point
        real, intent(in) :: tem
        real, dimension(count_point,size(env%tem_grid)), intent(in) :: table_search
        real, dimension(count_point), intent(in) :: table_res
        real :: res
        
        real :: l_lim_v, r_lim_v, l_lim_d, r_lim_d    ! Bracketing result values and search coordinates.    
        integer :: left, right, mid,tem_index         ! Binary-search bounds/midpoint and temperature index.    

        do tem_index = 1, size(env%tem_grid)
                if (tem == env%tem_grid(tem_index)) exit
        end do

        ! Guard against empty or single-point tables.
        if (count_point < 2) then
            if (count_point == 1) then
                res = table_res(1)
                return
            else
                print *, "ERROR: count_point =", count_point
                stop
            end if
        end if

        ! Clamp requests outside the tabulated coordinate range.
        if (num <= table_search(1,tem_index)) then
            res = table_res(1)
            return
        end if

        if (num >= table_search(count_point,tem_index)) then
            res = table_res(count_point)
            return
        end if

        ! Locate the bracketing interval by binary search.
        left = 1
        right = count_point

        do while (right - left > 1)
            mid = (left + right) / 2
            if (num <= table_search(mid,tem_index)) then
                right = mid
            else
                left = mid
            end if
        end do

        ! Perform linear interpolation inside the bracketing interval.
        l_lim_v = table_res(left)
        r_lim_v = table_res(right)
        l_lim_d = table_search(left,tem_index)
        r_lim_d = table_search(right,tem_index)

        if (abs(r_lim_d - l_lim_d) < 1.0e-30) then
            res = (l_lim_v + r_lim_v) / 2.0
        else
            res = l_lim_v + (r_lim_v - l_lim_v) * (num - l_lim_d) / (r_lim_d - l_lim_d)
        end if

    end function found_number_with_table

end module math_operation
