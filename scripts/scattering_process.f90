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

    function get_one_bump_netron_termalization(cur_netron_data, mass_enviroment, section_d,nuc,tem,env) result(new_netron_data)
        type(netron_data), intent(in) :: cur_netron_data        !Нейтрон до взаимодействия                                              IN
        real, intent(in) :: mass_enviroment                     !Масса ядра, на котором происходит рассеивание                          IN
        real, intent(in) :: section_d                           !Сечение рассеивания, отвечающее данной энергии нейтрона и ядру         IN
        type(nuclear_data), intent(in) :: nuc
        real, intent(in) :: tem
        type(enviroment), intent(in) :: env
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
        random_velocity = found_number_with_table(rundom_probability_speed_of_nuclie,nuc%coordinate_distribution_grid,nuc%coordinate_velocity_grid,nuc%count_point_vel_distr, tem, env)

        vn0(1) = cur_netron_data%dir(1)*cur_netron_data%speed
        vn0(2) = cur_netron_data%dir(2)*cur_netron_data%speed
        vn0(3) = cur_netron_data%dir(3)*cur_netron_data%speed
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
        new_netron_data%dir(1) = rel_dir(1)
        new_netron_data%dir(2) = rel_dir(2)
        new_netron_data%dir(3) = rel_dir(3)

        new_netron_data%speed = abs_v_rel
        new_netron_data%energy = 0.5 * (v_rel(1)**2 + v_rel(2)**2 + v_rel(3)**2)*10.35/10e8

        new_netron_data = get_one_bump_netron_slow_down(new_netron_data,mass_enviroment,section_d)


        vn_star_new(1) = new_netron_data%speed * new_netron_data%dir(1)
        vn_star_new(2) = new_netron_data%speed * new_netron_data%dir(2)
        vn_star_new(3) = new_netron_data%speed * new_netron_data%dir(3)

        vn_lab_new(1) = vn_star_new(1) + va0(1)
        vn_lab_new(2) = vn_star_new(2) + va0(2)
        vn_lab_new(3) = vn_star_new(3) + va0(3)

        abs_vn_lab_new = sqrt(vn_lab_new(1)**2 + vn_lab_new(2)**2 + vn_lab_new(3)**2)

        new_netron_data%speed = abs_vn_lab_new

        new_netron_data%dir(1) = vn_lab_new(1) / abs_vn_lab_new
        new_netron_data%dir(2) = vn_lab_new(2) / abs_vn_lab_new
        new_netron_data%dir(3) = vn_lab_new(3) / abs_vn_lab_new

        new_netron_data%energy = 0.5 * abs_vn_lab_new**2 * 10.35 / 10e8

    end function get_one_bump_netron_termalization

    function move_neutron_TMS_candidate(neutron, env, Tbase, Tmax,Sigma_maj) result(new_netron_d)
        
        type(netron_data), intent(in) :: neutron
        type(enviroment), intent(in) :: env

        real, intent(in) :: Tbase
        real, intent(in) :: Tmax

        real, intent(in) :: Sigma_maj
        type(netron_data)::new_netron_d
        real :: path
        real :: xi

       new_netron_d = neutron
        call random_number(xi)

        if (xi <= tiny(1.0)) then
            xi = tiny(1.0)
        end if

        path = -log(xi) / Sigma_maj

        new_netron_d%pos(1) = neutron%pos(1) + &
                        neutron%dir(1) * path

        new_netron_d%pos(2) = neutron%pos(2) + &
                        neutron%dir(2) * path

        new_netron_d%pos(3) = neutron%pos(3) + &
                        neutron%dir(3) * path


        if (neutron%speed > tiny(1.0)) then

            new_netron_d%life_time = neutron%life_time + &
                                path / neutron%speed * 0.01

        end if

    end function move_neutron_TMS_candidate
 
    function change_dir_TMS(cur_netron_data, mass_enviroment, target_velocity) result(new_netron_d)
        implicit none

        type(netron_data), intent(in) :: cur_netron_data
        real, intent(in) :: mass_enviroment
        real, intent(in) :: target_velocity(3)

        type(netron_data) :: new_netron_d

        real :: vn_lab(3)
        real :: v_rel(3)
        real :: v_cm(3)

        real :: v_rel_abs
        real :: v_rel_new(3)
        real :: vn_lab_new(3)
        real :: abs_vn_lab_new

        real :: e1(3), e2(3), e3(3)
        real :: norm_e1

        real :: mu
        real :: phi
        real :: sin_theta
        real :: u

        new_netron_d = cur_netron_data
        vn_lab = cur_netron_data%speed * cur_netron_data%dir
        v_rel = vn_lab - target_velocity
        v_rel_abs = sqrt(sum(v_rel**2))
        if (v_rel_abs <= tiny(1.0)) then
            return
        end if
        v_cm = (vn_lab + mass_enviroment * target_velocity) / &
            (1.0 + mass_enviroment)
        e3 = v_rel / v_rel_abs
        if (abs(e3(3)) < 0.9) then
            e1 = [-e3(2), e3(1), 0.0]
        else
            e1 = [0.0, -e3(3), e3(2)]
        end if

        norm_e1 = sqrt(sum(e1**2))

        if (norm_e1 <= tiny(1.0)) then
            return
        end if

        e1 = e1 / norm_e1

        e2(1) = e3(2)*e1(3) - e3(3)*e1(2)
        e2(2) = e3(3)*e1(1) - e3(1)*e1(3)
        e2(3) = e3(1)*e1(2) - e3(2)*e1(1)


        call random_number(u)
        mu = 2.0*u - 1.0

        call random_number(u)
        phi = 2.0 * acos(-1.0) * u

        sin_theta = sqrt(max(0.0, 1.0 - mu*mu))

        v_rel_new = v_rel_abs * ( &
            sin_theta*cos(phi)*e1 &
            + sin_theta*sin(phi)*e2 &
            + mu*e3 )


        vn_lab_new = v_cm + &
            (mass_enviroment / (mass_enviroment + 1.0)) * v_rel_new


        abs_vn_lab_new = sqrt(sum(vn_lab_new**2))

        if (abs_vn_lab_new <= tiny(1.0)) then
            return
        end if


        new_netron_d%dir = vn_lab_new / abs_vn_lab_new


        new_netron_d%speed = abs_vn_lab_new

        new_netron_d%energy = &
            0.5 * abs_vn_lab_new**2 * 10.35 / 10e8


        new_netron_d%count_collision = &
            cur_netron_data%count_collision + 1

    end function change_dir_TMS
end module scattering_process
