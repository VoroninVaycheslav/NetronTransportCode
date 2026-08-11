module scattering_process

    use data_type
    use math_operation
    implicit none

    contains
    !Моделирование рассеяния частиц
    function get_one_bump_netron_slow_down(cur_netron_data, mass_enviroment, section_d) result(new_netron_data)
        type(netron_data), intent(in) :: cur_netron_data        !Нейтрон до взаимодействия                                              IN
        real, intent(in) :: mass_enviroment                     !Масса ядра, на котором происходит рассеивание                          IN
        real, intent(in) :: section_d                           !Сечение рассеивания, отвечающее данной энергии нейтрона и ядру         IN
        type(netron_data) :: new_netron_data                    !Нейтрон после взаимодействия                                           OUT
       
        real :: T, newT, fi, P
        real :: norm_factor, eps
        
        eps = 1.0e-10
        T = get_random_in_range(-1.0, 1.0)
        newT = (1.0 + mass_enviroment * T) / sqrt(1.0 + mass_enviroment**2 + 2.0 * mass_enviroment * T)

        !Задаем случайные числа для рассивания (два угла рассеяния)
        fi = get_random_in_range(0.0, 6.283185)
        P = get_random_in_range(0.001, 1.0)
        
        ! Расчет нового направления
        if (abs(1.0 - cur_netron_data%dir%z**2) > eps) then
            norm_factor = sqrt(1.0 - newT**2) / sqrt(1.0 - cur_netron_data%dir%z**2)
            
            new_netron_data%dir%x = norm_factor * (cur_netron_data%dir%x * cur_netron_data%dir%z * cos(fi) - cur_netron_data%dir%y * sin(fi)) + cur_netron_data%dir%x * newT
            new_netron_data%dir%y = norm_factor * (cur_netron_data%dir%y * cur_netron_data%dir%z * cos(fi) + cur_netron_data%dir%x * sin(fi)) + cur_netron_data%dir%y * newT
            new_netron_data%dir%z = -sqrt(1.0 - newT**2) * sqrt(1.0 - cur_netron_data%dir%z**2) * cos(fi) + cur_netron_data%dir%z * newT
        else
            new_netron_data%dir%x = sqrt(1.0 - newT**2) * cos(fi)
            new_netron_data%dir%y = sqrt(1.0 - newT**2) * sin(fi)
            new_netron_data%dir%z = newT * sign(1.0, cur_netron_data%dir%z)
        end if
        ! Расчет новой энергии
        new_netron_data%energy = cur_netron_data%energy * (mass_enviroment**2 + 1.0 + 2.0 * mass_enviroment * T) / ((mass_enviroment + 1.0)**2)
        
        ! Расчет нового положения
        new_netron_data%pos%x = cur_netron_data%pos%x + new_netron_data%dir%x * (-log(P)) / section_d
        new_netron_data%pos%y = cur_netron_data%pos%y + new_netron_data%dir%y * (-log(P)) / section_d
        new_netron_data%pos%z = cur_netron_data%pos%z + new_netron_data%dir%z * (-log(P)) / section_d
        
        ! Расчет скорости и времени жизни
        new_netron_data%speed = 1.38e4 * sqrt(new_netron_data%energy)
        new_netron_data%life_time = cur_netron_data%life_time + (-log(P)) / (section_d * new_netron_data%speed) * 0.01
        new_netron_data%count_collision = cur_netron_data%count_collision + 1
        
        
        
    end function get_one_bump_netron_slow_down

    function get_one_bump_netron_termalization(cur_netron_data, mass_enviroment, section_d,nuc) result(new_netron_data)
        type(netron_data), intent(in) :: cur_netron_data        !Нейтрон до взаимодействия                                              IN
        real, intent(in) :: mass_enviroment                     !Масса ядра, на котором происходит рассеивание                          IN
        real, intent(in) :: section_d                           !Сечение рассеивания, отвечающее данной энергии нейтрона и ядру         IN
        type(nuclear_data), intent(in) :: nuc
        type(netron_data) :: new_netron_data                    !Нейтрон после взаимодействия                                           OUT

        ! Начальные скорости нейтрона и ядра в лабораторной системе
        real :: vn0(3), va0(3)
        real :: rel_dir(3), v_rel(3), u_n(3)
        real :: abs_v_rel, abs_va
        real :: P,T_a,fi_a
        real :: random_velocity
        real :: rundom_probability_speed_of_nuclie
        real :: vn_star_new(3), vn_lab_new(3), abs_vn_lab_new
        new_netron_data = cur_netron_data
        rundom_probability_speed_of_nuclie = get_random_in_range(0.0,0.999)
        random_velocity = found_number_with_table(rundom_probability_speed_of_nuclie,nuc%coordinate_distribution_grid,nuc%coordinate_velocity_grid,nuc%count_point_vel_distr)

        vn0(1) = cur_netron_data%dir%x*cur_netron_data%speed
        vn0(2) = cur_netron_data%dir%y*cur_netron_data%speed
        vn0(3) = cur_netron_data%dir%z*cur_netron_data%speed
        T_a = get_random_in_range(-1.0, 1.0)
        fi_a = get_random_in_range(0.0, 6.283185)
        P = get_random_in_range(0.001, 1.0)

        va0(1) = random_velocity*cos(fi_a)*T_a
        va0(2) = random_velocity*sin(fi_a)*T_a
        va0(3) = random_velocity*sqrt(1-T_a**2)
        abs_va = sqrt(va0(1)**2 + va0(2)**2 + va0(3)**2)
        v_rel = vn0 - va0
        abs_v_rel = sqrt(v_rel(1)**2 + v_rel(2)**2 + v_rel(3)**2)
        rel_dir = v_rel/abs_v_rel
        new_netron_data%dir%x = rel_dir(1)
        new_netron_data%dir%y = rel_dir(2)
        new_netron_data%dir%z = rel_dir(3)

        new_netron_data%speed = abs_v_rel
        new_netron_data%energy = 0.5 * (v_rel(1)**2 + v_rel(2)**2 + v_rel(3)**2)*10.35/10e8

        new_netron_data = get_one_bump_netron_slow_down(new_netron_data,mass_enviroment,section_d)


        vn_star_new(1) = new_netron_data%speed * new_netron_data%dir%x
        vn_star_new(2) = new_netron_data%speed * new_netron_data%dir%y
        vn_star_new(3) = new_netron_data%speed * new_netron_data%dir%z

        vn_lab_new(1) = vn_star_new(1) + va0(1)
        vn_lab_new(2) = vn_star_new(2) + va0(2)
        vn_lab_new(3) = vn_star_new(3) + va0(3)

        abs_vn_lab_new = sqrt(vn_lab_new(1)**2 + vn_lab_new(2)**2 + vn_lab_new(3)**2)

        new_netron_data%speed = abs_vn_lab_new

        new_netron_data%dir%x = vn_lab_new(1) / abs_vn_lab_new
        new_netron_data%dir%y = vn_lab_new(2) / abs_vn_lab_new
        new_netron_data%dir%z = vn_lab_new(3) / abs_vn_lab_new

        new_netron_data%energy = 0.5 * abs_vn_lab_new**2 * 10.35 / 10e8

    end function get_one_bump_netron_termalization

end module scattering_process
