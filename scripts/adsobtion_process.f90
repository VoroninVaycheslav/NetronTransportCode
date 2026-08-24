!==============================================================================
! Module: adsobtion_process
!
! Purpose:
!   Implements neutron absorption by terminating the current neutron history.
!
!==============================================================================

module adsobtion_process

    use data_type
    implicit none

    contains
    !==============================================================================
    ! function: get_absorption
    !
    ! Purpose:
    !   terminate the current neutron history.
    !
    ! Parametr IN:
    !   cur_netron_d - Neutron state immediately before absorption.
    !
    ! Parametr OUT:  
    !   new_netron_d - Neutron state after absorption.
    !
    !==============================================================================

    function get_absorption(cur_netron_d) result(new_netron_d)
        type(netron_data), intent(in) :: cur_netron_d     
        type(netron_data) :: new_netron_d                 

        new_netron_d = cur_netron_d
        new_netron_d%is_died = .True. ! Mark the neutron history as terminated by absorpt
    end function get_absorption

end module adsobtion_process
