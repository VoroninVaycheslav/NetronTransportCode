!==============================================================================
! Module: data_type
!
! Purpose:
!   Defines the derived types shared by the transport, nuclear-data, collision,
!   scattering, absorption, and Doppler-treatment modules.
!
!==============================================================================

module data_type

    !==============================================================================
    ! Derived Types: enviroment
    !
    ! Purpose:
    !   Material/environment state shared by the transport kernel.
    !   Current implementation stores all nuclides in a homogeneous material object.
    !   Temperature-dependent tables are selected by exact equality against tem_grid.
    !
    !==============================================================================
    
    type :: enviroment
        integer :: count_nuclear = 0                                        ! Number of nuclide species in the material.
        type(nuclear_data),allocatable :: different_tipe_of_nuclear(:)      ! Allocatable array of nuclide data objects.
        real, allocatable :: tem_grid(:)                                    ! Material temperature [K].
        real :: tem = 294.0                                                 ! Material temperature grid [K].
    end type enviroment 


    !==============================================================================
    ! Derived Types: nuclear_data
    !
    ! Purpose:
    !   Nuclear data associated with one nuclide in the material.
    !   
    ! Cross-section storage convention:
    !   energy_point_in_table(i)                        -> energy coordinate
    !   cross_data(r)%cross_section_point_in_table(i,t) -> reaction r at temperature t
    ! 
    ! OTF storage convention:
    !   K_OTF(i,r,c) -> coefficient c for energy point i and reaction r.
    !
    !==============================================================================
    
    type nuclear_data 
        character(len=40) name_of_nuclie                                   ! Human-readable nuclide name from the input file.
        integer index_of_nuclie                                            ! Nuclide identifier read from the input file.
        real mass_of_nuclear                                               ! Target mass parameter.
        real nuclear_dencity                                               ! Number-density factor.
        integer count_process                                              ! Number of cross-section/reaction columns.
        integer count_point                                                ! Number of energy points in the tabulated data.
        integer :: count_point_vel_distr = 0                               ! Number of points in the target-speed lookup grid.
        real, allocatable :: energy_point_in_table(:)                      ! Energy grid shared by all reactions for this nuclide [eV].
        type(cross_section_data),allocatable:: cross_data(:)               ! Reaction-specific cross-section tables.
        real, dimension(:,:), allocatable ::coordinate_distribution_grid   ! Maxwell CDF lookup table versus speed and temperature.
        real, dimension(:), allocatable ::coordinate_velocity_grid         ! Speed coordinate associated with the Maxwell lookup table.
        real, dimension(:), allocatable ::e_uniq_grid                      ! Unique/adaptive energy grid used by the OTF representation.
        real, dimension(:,:,:), allocatable ::K_OTF                        ! OTF coefficients: energy x reaction x coefficient.
        
    end type nuclear_data

    !==============================================================================
    ! Derived Types: cross_section_data
    !
    ! Purpose:
    !   One reaction cross-section table for a nuclide.
    !
    !==============================================================================
    
    type cross_section_data
        integer :: index_of_process = -1                                   ! Reaction/process identifier from the input data.   
        real, allocatable :: cross_section_point_in_table(:,:)             ! Cross sections indexed by energy point and temperature.   
    end type cross_section_data

    !==============================================================================
    ! Derived Types: netron_data
    !
    ! Purpose:
    !   Complete state carried by one neutron history.
    !   Particle-history state is updated in place by successive collision kernels.
    !   is_died is the history-termination flag used by the main transport loops.
    !
    !==============================================================================
    
    type :: netron_data
        real :: dir(3)                                                     ! Unit direction-cosine vector.                
        real :: pos(3)                                                     ! Cartesian position; current transport convention uses cm.               
        real :: energy                                                     ! Neutron kinetic energy [eV].               
        real :: speed                                                      ! Neutron speed; current implementation uses m/s.                
        real :: life_time = 0                                              ! Accumulated neutron lifetime [s in the current convention].                
        integer :: count_collision = 0                                     ! Number of accepted collision events.                
        logical :: is_died = .False.                                       ! True after absorption terminates the history.                
        real :: speed_estimator = 0                                        ! Per-history local estimator accumulated by transport modes.
    end type netron_data

end module data_type