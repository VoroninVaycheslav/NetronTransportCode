!==============================================================================
! Neutron Transport Kernel - NJOY-style Doppler broadening kernel
!
! Purpose:
!   Evaluates the Doppler-broadened cross section using the existing analytical
!   formulation based on sigma_star, H_n functions, and recursive F_n moments.
!
! Important:
!   The implementation is documented here; its equations and numerical order
!   are not modified.
!==============================================================================
module DoplerTrearmentNJOY
   implicit none

contains
    
   !==============================================================================
   ! function: doplerBroadr
   !
   ! Purpose:
   !   Evaluate the legacy Doppler-broadened cross section at currentEnergy.
   !
   ! Parametr IN:
   !   currentEnergy    -   Incident neutron energy [eV].
   !   E                -   Original energy grid [eV].
   !   sigma            -   Original cross-section values on E.
   !   alpha            -   Doppler parameter used by the analytical transformation.
   !   n                -   Number of tabulated energy points.
   !
   ! Parametr OUT:  
   !   res              -   Broadened cross section.
   !
   !==============================================================================

   function doplerBroadr(currentEnergy, E, sigma, alpha, n) result(res)
      integer, intent(in) :: n
      real, intent(in) :: alpha
      real, intent(in) :: currentEnergy, E(n), sigma(n)
      real :: res
      real y                     ! Transformed incident-energy coordinate sqrt(alpha*E).
      
      y = sqrt(alpha * currentEnergy)

      res = sigma_star(y, E, sigma, alpha, n) - sigma_star(-y, E, sigma, alpha, n)

   end function doplerBroadr

   !==============================================================================
   ! function: sigma_star
   !
   ! Purpose:
   !   Evaluate the accumulated sigma-star contribution for a signed y value.
   !
   ! Parametr IN:
   !   y        -   Signed transformed-energy coordinate.
   !   E        -   Energy grid [eV].
   !   sigma    -   Cross-section values.
   !   alpha    -   Doppler transformation parameter.
   !   n        -   Number of grid points.
   !
   ! Parametr OUT:  
   !   res      -   Accumulated analytical contribution.
   !
   !==============================================================================

   function sigma_star(y, E, sigma,alpha, n) result(res)
      integer, intent(in) :: n
      real, intent(in) :: alpha
      real, intent(in) :: y, E(n), sigma(n)
      real :: res
   
      integer :: i                     ! Energy-interval loop index.
      real :: x_i, x_ip1               ! Transformed left/right interval coordinates.
      real :: s_i, slope               ! Left cross section and linear slope versus x^2.
      real :: H0,H1,H2,H3,H4           ! Analytical H_n interval moments.
      real :: A_i, B_i                 ! Interval coefficients in the broadened-cross-section sum.
   
      res = 0.0
   
      do i = 1, n-1
   
         ! Transform the tabulated energy coordinate to x = sqrt(alpha*E).
         x_i   = sqrt(alpha * E(i))
         x_ip1 = sqrt(alpha * E(i+1))
   
         ! The legacy derivation assumes linearity in x^2 on each interval.
         s_i = sigma(i)
         if (abs(x_ip1**2 - x_i**2) < 1.0e-12) then
            slope = 0.0
         else
            slope = (sigma(i+1)-sigma(i)) / (x_ip1**2 - x_i**2)
         end if
   
         ! Evaluate H0..H4 over the transformed interval.
         call compute_H_all(x_i - y, x_ip1 - y, H0,H1,H2,H3,H4)
   
         ! Assemble interval coefficients used by the analytical expression.
         A_i = (1.0/y**2)*H2 + (2.0/y)*H1 + H0
         B_i = (1.0/y**2)*H4 + (4.0/y)*H3 + 6.0*H2 + 4.0*y*H1 + y*y*H0
   
         ! Accumulate the legacy analytical expression (formula 88 in the source derivation).
         res = res + A_i*(s_i - slope*x_i**2) + B_i*slope
   
      end do
   
   end function

   !==============================================================================
   ! subroutine: compute_H_all
   !
   ! Purpose:
   !   Compute H0..H4 as differences of F_n values at interval endpoints.
   !
   !==============================================================================

  subroutine compute_H_all(a,b,H0,H1,H2,H3,H4)
      real, intent(in) :: a,b
      real, intent(out) :: H0,H1,H2,H3,H4
  
      H0 = Fn(0,a) - Fn(0,b)
      H1 = Fn(1,a) - Fn(1,b)
      H2 = Fn(2,a) - Fn(2,b)
      H3 = Fn(3,a) - Fn(3,b)
      H4 = Fn(4,a) - Fn(4,b)
  end subroutine
  
   !==============================================================================
   ! function: Fn
   !
   ! Purpose:
   !   Evaluate the recursive F_n moment used by the broadening formula.
   !
   !==============================================================================
  
   recursive function Fn(n,a) result(res)
      integer, intent(in) :: n
      real, intent(in) :: a
      real :: res

      if (n == 0) then
         res = 0.5 * erfc(a)
      else if (n == 1) then
         res = exp(-a*a)/(2.0*sqrt(3.141592653589793))
      else
         res = (real(n-1)/2.0)*Fn(n-2,a) + a**(n-1)*Fn(1,a)
      end if
  
   end function

end module DoplerTrearmentNJOY