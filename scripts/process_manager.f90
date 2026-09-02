!==============================================================================
! module: process_manager
!
! Purpose:
!   Selects the collision nuclide and reaction from macroscopic/microscopic
!   cross sections, dispatches scattering or absorption physics, and provides
!   the cross-section treatment selector used by the collision kernel.
!
!==============================================================================
module process_manager

    use data_type
    use serpent_dopler
    use scattering_process
    use adsobtion_process
    use operation_with_data
    use math_operation
    use DoplerTrearmentNJOY
    use Dopler_MCNP_OTF
    implicit none

    contains
    !==============================================================================
    ! function: collision_controller
    !
    ! Purpose:
    !   Execute one collision-selection and reaction-dispatch step.
    !
    ! Selection sequence:
    !   1. Build per-nuclide macroscopic total cross sections.
    !   2. Sample the interacting nuclide by cumulative probability.
    !   3. Optionally perform the TMS acceptance/rejection branch.
    !   4. Build reaction microscopic cross sections for the selected nuclide.
    !   5. Sample and dispatch scattering or absorption.
    !
    ! Parametr IN:
    !   cur_netron_data     -   Neutron state at entry t
    !   env                 -   Material/nuclear-data environment.
    !
    ! Parametr OUT:  
    !   new_netron_d        -   Updated neutron state after scattering, absorption, or a TMS rejection cycle.
    !
    !==============================================================================

    function collision_controller(cur_netron_data, env,typeW) result(new_netron_data)
        
        type(netron_data), intent(in) :: cur_netron_data                        
        type(enviroment), intent(in) :: env 
        integer, intent(in) :: typeW                                    
        type(netron_data) :: new_netron_data                                    
        
        real, allocatable :: mic_current_cross_section_mas(:)                   ! Reaction microscopic cross sections for the selected nuclide.
        real, allocatable :: total_mac_current_cross_section_mas(:)             ! Per-nuclide macroscopic total cross sections.
        real :: total_mac_cross_section                                         ! Total macroscopic cross section of the material.
        real :: total_micro_cross_section                                       ! Sum of sampled reaction microscopic cross sections.
        integer :: type_of_nuclie                                               ! Index of the sampled interacting nuclide.
        integer :: type_of_process                                              ! Index of the sampled reaction within the local reaction array.
        real :: sum_val,pacc,dTlocal,Eprime,Stot_prime,g_local                  ! Cumulative XS, TMS acceptance, dT, relative energy, XS, and g factor.
        real :: p1, p2,u, path                                                  ! Sampling coordinates, uniform deviate, and traveled path length.
        integer :: i                                                            ! Nuclide/reaction loop index.
        ! Cross-section treatment used internally by collision_controller.
        ! 0 = direct tabulated interpolation
        ! 1 = NJOY-style analytical broadening
        ! 2 = OTF fitted expansion  [current preserved default]
        ! 3 = TMS majorant / rejection branch
        ! NOTE: this remains local and hard-coded to preserve validated behavior.
       
        real target_vel(3)
        new_netron_data = cur_netron_data
        do
            total_mac_cross_section = 0
            
            allocate(total_mac_current_cross_section_mas(env%count_nuclear))
            
            ! Build macroscopic total cross sections for all nuclides.
            do i = 1, env%count_nuclear
                
                total_mac_current_cross_section_mas(i) = get_sigma(new_netron_data%energy, env,i,1,typeW)*env%different_tipe_of_nuclear(i)%nuclear_dencity
            end do
            
            ! Sum macroscopic cross sections before probabilistic nuclide selection.
            total_mac_cross_section = sum(total_mac_current_cross_section_mas)
            
            ! Draw a random coordinate on the total macroscopic-cross-section interval.
            p1 = get_random_in_range(0.0000000000001, total_mac_cross_section)
            type_of_nuclie = 0
            sum_val = 0.0
            ! Select the nuclide by cumulative cross-section weight.
            do while (sum_val < p1 .and. type_of_nuclie < env%count_nuclear)
                type_of_nuclie = type_of_nuclie + 1
                sum_val = sum_val + total_mac_current_cross_section_mas(type_of_nuclie)
            end do
            
            if (type_of_nuclie < 1 .or. type_of_nuclie > env%count_nuclear) then
                print *, "ОШИБКА: Неверный тип ядра после выбора: ", type_of_nuclie, p1,total_mac_current_cross_section_mas
                new_netron_data = cur_netron_data
                deallocate(total_mac_current_cross_section_mas)
                return
            end if



            ! TMS-only candidate-collision acceptance/rejection branch.
            if(typeW == 3)then
                dTlocal = env%tem-294
                if (dTlocal > 0.0) then
                    call tms_sample_target_energy( &
                    new_netron_data%energy, dTlocal, env%different_tipe_of_nuclear(type_of_nuclie)%mass_of_nuclear, Eprime,target_vel)
                 else
                    Eprime = new_netron_data%energy
                end if
                
                Stot_prime = found_cross_section_from_energy( &
                        Eprime, &
                        294.0, &
                        env, &
                        env%different_tipe_of_nuclear(type_of_nuclie), &
                        1) * &
                        env%different_tipe_of_nuclear(type_of_nuclie)%nuclear_dencity

                g_local = tms_g(new_netron_data%energy, dTlocal, env%different_tipe_of_nuclear(type_of_nuclie)%mass_of_nuclear)

                pacc = g_local * Stot_prime / total_mac_current_cross_section_mas(type_of_nuclie)
                
                pacc = min(1.0, max(0.0, pacc))
                call random_number(u)
                new_netron_data =  move_neutron_TMS_candidate(new_netron_data,env,294.0,3200.0,total_mac_cross_section)
                    
                if(u>pacc)then
                    deallocate(total_mac_current_cross_section_mas)
                    cycle
                
                end if
            end if


            ! Evaluate reaction cross sections (table columns 2..count_process).
            allocate(mic_current_cross_section_mas(env%different_tipe_of_nuclear(type_of_nuclie)%count_process-1))
            
            total_micro_cross_section = 0.0
            do i = 2, env%different_tipe_of_nuclear(type_of_nuclie)%count_process
                if(typeW == 3)then
                    
                     mic_current_cross_section_mas(i-1) = &
                        found_cross_section_from_energy( &
                            Eprime, &
                            294.0, &
                            env, &
                            env%different_tipe_of_nuclear(type_of_nuclie), &
                            i)
                else
                    mic_current_cross_section_mas(i-1) = get_sigma(new_netron_data%energy, env,type_of_nuclie,i,typeW)
                end if
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
                if(typeW == 3)then

                    new_netron_data%speed_estimator = g_local*found_cross_section_from_energy(Eprime,294.0,env,env%different_tipe_of_nuclear(type_of_nuclie),2)*env%different_tipe_of_nuclear(type_of_nuclie)%nuclear_dencity/total_mac_cross_section
                    new_netron_data = change_dir_TMS(new_netron_data, env%different_tipe_of_nuclear(type_of_nuclie)%mass_of_nuclear,target_vel)
                else
                    ! Рассеяние с дыигвющимися ядрами
                    new_netron_data = get_one_bump_netron_termalization(new_netron_data, env%different_tipe_of_nuclear(type_of_nuclie)%mass_of_nuclear, total_mac_cross_section,env%different_tipe_of_nuclear(type_of_nuclie),env%tem,env)
                    
                    ! Рассеяние с покоящимися ядрами
                    !new_netron_data = get_one_bump_netron_slow_down(new_netron_data, env%different_tipe_of_nuclear(type_of_nuclie)%mass_of_nuclear, total_mac_cross_section)
                    path = sqrt((new_netron_data%pos(1)-cur_netron_data%pos(1))**2+(new_netron_data%pos(2)-cur_netron_data%pos(2))**2+(new_netron_data%pos(3)-cur_netron_data%pos(3))**2)
                    new_netron_data%speed_estimator = path * mic_current_cross_section_mas(2)
                end if
                case(2)
                    new_netron_data = get_absorption(new_netron_data)
            case default
                print *, "ВНИМАНИЕ: Неизвестный тип процесса: ", type_of_process, p2, cur_netron_data%energy
            end select
            
            deallocate(mic_current_cross_section_mas, total_mac_current_cross_section_mas)
            exit
        end do
    end function collision_controller
    
    !==============================================================================
    ! function: give_random_direction
    !
    ! Purpose:
    !   Initialize N neutron histories with random directions and common energy E.
    !
    ! Parametr IN:
    !   N            -   Number of neutron histories.
    !   E            -   Initial neutron energy [eV].
    !
    ! Parametr OUT:  
    !   new_netron_d -   Array of initialized neutron states at the origin.
    !
    !==============================================================================

    function give_random_direction(N, E) result(new_netron_d)
        integer, intent(in) :: N                        
        real, intent(in):: E
        type(netron_data) :: new_netron_d(N)            
        real :: T, Fi
        integer :: i

        do i = 1, N
            T = get_random_in_range(0.0, 3.14159)
            Fi = get_random_in_range(0.0, 6.28318)
            new_netron_d(i)%dir(1) = sin(T) * cos(Fi)
            new_netron_d(i)%dir(2) = sin(T) * sin(Fi)
            new_netron_d(i)%dir(3) = cos(T)
            new_netron_d(i)%energy = E
  
            new_netron_d(i)%speed = 1.38e4 * sqrt(new_netron_d(i)%energy)
            new_netron_d(i)%pos = [0.0, 0.0, 0.0]
            new_netron_d(i)%life_time = 0.0
            new_netron_d(i)%count_collision = 0
            new_netron_d(i)%is_died = .False.
        end do
    end function give_random_direction

    !==============================================================================
    ! function: get_sigma
    !
    ! Purpose:
    !   Collision-selection and reaction-dispatch logic.
    !
    !
    ! Parametr IN:
    !   energy          -   Neutron energy [eV].
    !   env             -   Material/nuclear-data environment.
    !   type_of_nuclie  -   Nuclide index.
    !   type_of_process -   Internal reaction-column index
    !   typeWorking     -   Cross-section treatment selector.
    !
    ! Parametr OUT:  
    !   sigma           -   Microscopic cross section.
    !
    !==============================================================================

    function get_sigma(energy, env, type_of_nuclie, type_of_process, typeWorking) result(sigma)
        real, intent(in) :: energy
        type(enviroment), intent(in) :: env
        real:: tem
        integer, intent(in) :: type_of_nuclie
        integer, intent(in) :: type_of_process
        integer, intent(in) :: typeWorking
        integer index_enrgy
        real :: sigma, alpha

        real :: kl(17)
        tem = env%tem
        index_enrgy = 1
        select case(typeWorking)
            case(0)
                sigma = found_cross_section_from_energy(energy,tem,env, env%different_tipe_of_nuclear(type_of_nuclie),type_of_process)
            case(1)
                alpha = env%different_tipe_of_nuclear(type_of_nuclie)%mass_of_nuclear / (k * tem)
                sigma = doplerBroadr(energy, env%different_tipe_of_nuclear(type_of_nuclie)%energy_point_in_table,&
                env%different_tipe_of_nuclear(type_of_nuclie)%cross_data(type_of_process)%cross_section_point_in_table, &
                alpha, size(env%different_tipe_of_nuclear(type_of_nuclie)%energy_point_in_table))
            case(2) 
                do while(energy>env%different_tipe_of_nuclear(type_of_nuclie)%e_uniq_grid(index_enrgy))
                    index_enrgy = index_enrgy + 1
                end do
                
                kl =env%different_tipe_of_nuclear(type_of_nuclie)%K_OTF(index_enrgy,type_of_process,1:17)
                sigma = eval_expansion(tem,300.0,3200.0,8,kl)
                !sigma = doplerBroadrOTF_MCNP(env%different_tipe_of_nuclear(type_of_nuclie), tem, energy)
                
            case(3)
                sigma = tms_nuclide_majorant(env%different_tipe_of_nuclear(type_of_nuclie),energy,294.0,3200.0,env,type_of_process)
        end select
    end function get_sigma


end module process_manager
