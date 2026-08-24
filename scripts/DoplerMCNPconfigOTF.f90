!==============================================================================
! Module: MCNPconfigOTF
!
! Purpose:
!   Stores parameters used to construct temperature/energy grids and fit the
!   temperature-dependent OTF cross-section representation.
!==============================================================================

module MCNPconfigOTF
   implicit none
   !==============================================================================  
   ! Derived Types: otf_config_type
   !
   ! Purpose:
   !   Saves the configuration with settings for modeling coefficients 
   !   using the MCNP OTF method
   !
   !==============================================================================
    
   type, public :: otf_config_type
      real :: awr = 238                      ! Atomic-weight-ratio parameter used by OTF preprocessing.
      real :: t_base = 293                   ! Reference/base temperature [K].
      real :: t_min = 300                    ! Minimum fitting/refinement temperature [K].
      real :: t_max = 3200.0                 ! Maximum fitting/refinement temperature [K].
      real :: dt_union = 50.0                ! Temperature spacing for union-grid refinement [K].
      real :: dt_fit = 10.0                  ! Temperature spacing for coefficient fitting [K].
      real :: ft = 1.0e-3                    ! Fractional tolerance for energy-grid refinement.
      real :: fit_ft = 1.0e-3                ! Reserved legacy fit tolerance parameter.
      real :: abs_xs_tol = 1.0e-12           ! Absolute cross-section floor/tolerance parameter.
      real :: energy_to_ev = 1.0             ! Energy conversion factor to eV.
      real :: e_min = 0.0                    ! Lower energy bound [eV].
      real :: e_max = 100.0                  ! Upper energy bound [eV].
   end type otf_config_type
contains
    
end module MCNPconfigOTF