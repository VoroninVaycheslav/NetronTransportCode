module adsobtion_process

    use data_type
    implicit none

    contains
    !Моделирование процесса поглащений нейтрона ядрном
    function get_absorption(cur_netron_d) result(new_netron_d)
        type(netron_data), intent(in) :: cur_netron_d       !Текущий нейтрон                    IN
        type(netron_data) :: new_netron_d                   !Нейтрон после взаимодействия       OUT

        new_netron_d = cur_netron_d
        new_netron_d%is_died = .True.   !Делаем нейтрон поглащенным
    end function get_absorption

end module adsobtion_process
