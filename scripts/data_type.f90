module data_type

    implicit none
    
    !Вектор трехмерного измерения
    type :: vector3
        real :: x
        real :: y
        real :: z
    end type vector3

    !Данные среды
    type :: enviroment
        integer :: count_nuclear = 0                                        ! Количество ядер
        type(nuclear_data),allocatable :: different_tipe_of_nuclear(:)      ! Массив ядер
    end type enviroment

    !Данные отдельного ядра в среде 
    type nuclear_data 
        character(len=40) name_of_nuclie                                    ! Имя нуклида
        integer index_of_nuclie                                             ! Индекс ядра
        real mass_of_nuclear                                                ! Масса ядра 
        real nuclear_dencity                                                ! Ядерная плотность
        integer count_process                                               ! Количество процессов
        integer count_point                                                 ! Количество точек в таблице
        integer :: count_point_vel_distr = 0
        real, allocatable :: energy_point_in_table(:)                       ! Столбец энергии (Общий для всех)
        type(cross_section_data),allocatable:: cross_data(:)                ! Столбы различных сечений
        real, dimension(:), allocatable ::coordinate_distribution_grid      ! Сетка по плотности вероятности ядер
        real, dimension(:), allocatable ::coordinate_velocity_grid          ! Сетка по скоростям ядер
        real, dimension(:), allocatable ::e_uniq_grid                       ! Сетка по уникальным энергиям (для OTF)
    end type nuclear_data

    !Таблица сечений для ядра 
    type cross_section_data
        integer :: index_of_process = -1                                    ! Индекс процесса
        real, allocatable :: cross_section_point_in_table(:)                ! Стобец процесса
    end type cross_section_data

    !информация о нейтроне
    type :: netron_data
        real :: dir(3)                                                      ! Направление движения
        real :: pos(3)                                                      ! Расположение в пространстве
        real :: energy                                                      ! Энергия нейтрона, eV
        real :: speed                                                       ! Скорость нейтрона
        real :: life_time = 0                                               ! Время рассеяния
        integer :: count_collision = 0                                      ! Количество испытанных столкновений
        logical :: is_died = .False.                                        ! Флаг поглощенности нейтрона 
    end type netron_data

end module data_type