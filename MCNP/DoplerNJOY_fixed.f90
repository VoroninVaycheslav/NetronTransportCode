module DoplerTrearmentNJOY
    implicit none

contains
    
! ====================================================
! σ*(y)
! ====================================================
function doplerBroadr(currentEnergy, E, sigma, alpha, n) result(res)
   integer, intent(in) :: n
   real, intent(in) :: alpha
   real, intent(in) :: currentEnergy, E(n), sigma(n)
   real :: res
   real y
   
   y = sqrt(alpha * currentEnergy)

   res = sigma_star(y, E, sigma, alpha, n) - sigma_star(-y, E, sigma, alpha, n)

   end function
function sigma_star(y, E, sigma,alpha, n) result(res)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: alpha
    real, intent(in) :: y, E(n), sigma(n)
    real :: res
  
    integer :: i
    real :: x_i, x_ip1
    real :: s_i, slope
    real :: H0,H1,H2,H3,H4
    real :: A_i, B_i
  
    res = 0.0
  
    do i = 1, n-1
  
       ! переход к x
       x_i   = sqrt(alpha * E(i))
       x_ip1 = sqrt(alpha * E(i+1))
  
       ! ЛИНЕЙНОСТЬ ПО x^2 (ЭТО КЛЮЧЕВОЕ)
       s_i = sigma(i)

       ! Repeated energy points are allowed in tabulated cross sections.
       ! They represent a zero-width interval, so do not evaluate a slope
       ! across that interval (division by zero would produce NaN).
       if (abs(x_ip1**2 - x_i**2) <= 100.0*epsilon(1.0) * &
           max(1.0, abs(x_i**2), abs(x_ip1**2))) cycle

       slope = (sigma(i+1)-sigma(i)) / (x_ip1**2 - x_i**2)
  
       ! H функции
       call compute_H_all(x_i - y, x_ip1 - y, H0,H1,H2,H3,H4)
  
       ! коэффициенты
       A_i = (1.0/y**2)*H2 + (2.0/y)*H1 + H0
       B_i = (1.0/y**2)*H4 + (4.0/y)*H3 + 6.0*H2 + 4.0*y*H1 + y*y*H0
  
       ! формула (88)
       res = res + A_i*(s_i - slope*x_i**2) + B_i*slope
  
    end do
  
  end function
  
  ! ====================================================
  ! H_n = F_n(a) - F_n(b)
  ! ====================================================
  subroutine compute_H_all(a,b,H0,H1,H2,H3,H4)
    implicit none
    real, intent(in) :: a,b
    real, intent(out) :: H0,H1,H2,H3,H4
  
    H0 = Fn(0,a) - Fn(0,b)
    H1 = Fn(1,a) - Fn(1,b)
    H2 = Fn(2,a) - Fn(2,b)
    H3 = Fn(3,a) - Fn(3,b)
    H4 = Fn(4,a) - Fn(4,b)
  end subroutine
  
  ! ====================================================
  ! F_n
  ! ====================================================
  recursive function Fn(n,a) result(res)
    implicit none
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