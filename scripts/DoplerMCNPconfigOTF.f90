module MCNPconfigOTF
    implicit none

    type, public :: otf_config_type
      real :: awr = 238
      real :: t_base = 293
      real :: t_min = 300
      real :: t_max = 3200.0
      real :: dt_union = 50.0
      real :: dt_fit = 10.0
      real :: ft = 1.0e-3
      real :: fit_ft = 1.0e-3
      real :: abs_xs_tol = 1.0e-12
      real :: energy_to_ev = 1.0
      real :: e_min = 0.0
      real :: e_max = 100.0
   end type otf_config_type
contains
    
end module MCNPconfigOTF