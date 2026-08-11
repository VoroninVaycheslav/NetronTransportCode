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
        if (abs(1.0 - cur_netron_data%dir(3)**2) > eps) then
            norm_factor = sqrt(1.0 - newT**2) / sqrt(1.0 - cur_netron_data%dir(3)**2)
            
            new_netron_data%dir(1) = norm_factor * (cur_netron_data%dir(1) * cur_netron_data%dir(3) * cos(fi) - cur_netron_data%dir(2) * sin(fi)) + cur_netron_data%dir(1) * newT
            new_netron_data%dir(2) = norm_factor * (cur_netron_data%dir(2) * cur_netron_data%dir(3) * cos(fi) + cur_netron_data%dir(1) * sin(fi)) + cur_netron_data%dir(2) * newT
            new_netron_data%dir(3) = -sqrt(1.0 - newT**2) * sqrt(1.0 - cur_netron_data%dir(3)**2) * cos(fi) + cur_netron_data%dir(3) * newT
        else
            new_netron_data%dir(1) = sqrt(1.0 - newT**2) * cos(fi)
            new_netron_data%dir(2) = sqrt(1.0 - newT**2) * sin(fi)
            new_netron_data%dir(3) = newT * sign(1.0, cur_netron_data%dir(3))
        end if
        ! Расчет новой энергии
        new_netron_data%energy = cur_netron_data%energy * (mass_enviroment**2 + 1.0 + 2.0 * mass_enviroment * T) / ((mass_enviroment + 1.0)**2)
        
        ! Расчет нового положения
        new_netron_data%pos(1) = cur_netron_data%pos(1) + new_netron_data%dir(1) * (-log(P)) / section_d
        new_netron_data%pos(2) = cur_netron_data%pos(2) + new_netron_data%dir(2) * (-log(P)) / section_d
        new_netron_data%pos(3) = cur_netron_data%pos(3) + new_netron_data%dir(3) * (-log(P)) / section_d
        
        ! Расчет скорости и времени жизни
        new_netron_data%speed = 1.38e4 * sqrt(new_netron_data%energy)
        new_netron_data%life_time = cur_netron_data%life_time + (-log(P)) / (section_d * new_netron_data%speed) * 0.01
        new_netron_data%count_collision = cur_netron_data%count_collision + 1
        
        
        
    end function get_one_bump_netron_slow_down

    

end module scattering_process
