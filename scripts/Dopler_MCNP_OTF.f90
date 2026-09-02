!==============================================================================
! Module: Dopler_MCNP_OTF
!
! Purpose:
!   Builds adaptive unionized energy grids, evaluates Doppler-broadened data on
!   temperature grids, and fits the existing temperature expansion used by the
!   transport kernel.
!
! Numerical behavior:
!   Adaptive-grid tolerances, basis functions, Jacobi-SVD solver, and expansion
!   evaluation are preserved exactly in computational content.
!==============================================================================

module Dopler_MCNP_OTF
   use MCNPconfigOTF
   use data_type
   use DoplerTrearmentNJOY
   
   implicit none

   real :: AWR = 238.0             ! Atomic-weight-ratio used by build_coefficients.
   real :: T = 1200.0              ! Module temperature parameter
   real :: k = 8.617333262e-5      ! Boltzmann constant [eV/K].

   contains
   !==============================================================================
   ! subroutine: build_temerature_grid
   !
   ! Purpose:
   !   Build an equally spaced temperature grid from t_min to t_max.
   !
   ! Parametr IN:
   !   t_min - Minimum temperature [K].
   !   t_max - Maximum temperature [K].
   !   delta - Temperature spacing [K].
   !
   ! Parametr IN/OUT:  
   !   t_grid - Allocatable output grid.
   !
   !==============================================================================

   subroutine build_temerature_grid(t_min,t_max, delta ,t_grid)
      real, allocatable, intent(inout) :: t_grid(:)
      real, intent(in) :: t_max
      real, intent(in) :: delta
      real, intent(in) :: t_min
      
      integer :: amount_point_in_temp_grid, i             ! Step of temterature and index of iteration
   
      amount_point_in_temp_grid = (t_max-t_min)/delta
      
      
      allocate(t_grid(amount_point_in_temp_grid+1))
      do i = 0, amount_point_in_temp_grid
          t_grid(i+1) = t_min + i*delta
         
        end do
   end subroutine

   !==============================================================================
   ! subroutine: build_union_grid
   !
   ! Purpose:
   !   Adaptively refine an energy grid against Doppler-broadened curvature.
   !
   ! Parametr IN:
   !   Tunion         -   Temperatures at which refinement tolerance is checked.
   !   nreac          -   Number of reaction columns considered.
   !   ft             -   Fractional midpoint interpolation tolerance.
   !   E_old          -   Original tabulated energy grid.
   !   env            -   Nuclear-data environment.
   !   typeOfNuclear  -   Nuclide index in env%different_tipe_of_nuclear.
   !   n              -   Number of original energy-grid points.
   !
   ! Parametr IN/OUT:  
   !   Egrid          -   Energy grid refined by midpoint insertion.
   !
   !==============================================================================

   subroutine build_union_grid(Egrid,Tunion,nreac,ft,E_old,env,typeOfNuclear,n)
     
      integer, intent(in) :: n
      real, intent(in) :: E_old(n)
      type(enviroment), intent(in) :: env
      real, allocatable, intent(inout) :: Egrid(:)
      real, intent(in) :: Tunion(:)
      real, intent(in) :: ft
      integer, intent(in) :: nreac
      integer, intent(in) :: typeOfNuclear
   
      integer :: i,j,r                  ! Energy-interval, temperature, and reaction indices.
   
      real :: E1,E2,Em                  ! Left, right, and midpoint energies [eV].
      real :: s1,s2,sm                  ! Broadened cross sections at E1, E2, and Em.
      real :: slin,err                  ! Linear midpoint estimate and fractional interpolation error.
      real :: alfa                      ! Doppler alpha parameter at the active temperature.
      logical :: add_point              ! True when the current interval must be refined.
   
      i = 1
      ! Refinement proceeds left-to-right. Inserting a midpoint keeps the interval
      ! index fixed so the new left subinterval is rechecked immediately.
      do while (i < size(Egrid))
         print*,i," ",size(Egrid)
         E1 = Egrid(i)
         E2 = Egrid(i+1)
   
         ! Midpoint of the current energy interval.
         Em = 0.5*(E1+E2)
   
         add_point = .false.
   
         ! Check all reaction channels included in the refinement criterion.
         do r = 1,nreac
   
            ! Check every temperature on the union-temperature grid.
            do j = 1,size(Tunion)
   
               alfa = env%different_tipe_of_nuclear(typeOfNuclear)%mass_of_nuclear/(k*Tunion(j))
               ! Evaluate broadened cross sections at both endpoints and the midpoint.
               s1 = doplerBroadr(E1,E_old,env%different_tipe_of_nuclear(typeOfNuclear)%cross_data(r)%cross_section_point_in_table,alfa,n)
               s2 = doplerBroadr(E2,E_old,env%different_tipe_of_nuclear(typeOfNuclear)%cross_data(r)%cross_section_point_in_table,alfa,n)
               sm = doplerBroadr(Em,E_old,env%different_tipe_of_nuclear(typeOfNuclear)%cross_data(r)%cross_section_point_in_table,alfa,n)
   
               ! Midpoint value predicted by linear interpolation.
               slin = 0.5*(s1+s2)
   
               ! Fractional tolerance.
               err = abs(sm-slin) / max(abs(sm),1.0e-12)
               if (err -ft > 1.0e-12 .AND. E2-E1>0.00001) then
                  print*, E1,E2
                  add_point = .true.
                  exit
               end if
   
            end do
   
            if (add_point) exit
   
         end do
   
         if (add_point) then
   
            ! Insert the midpoint when the fractional error exceeds tolerance.
            !
            ! Keep i unchanged after insertion:
            ! the newly created left interval [E1, Em] is checked first.
            call insert_point(Egrid,i+1,Em)
   
         else
   
            ! The interval satisfies the fractional-tolerance criterion.
            i = i+1
   
         end if
   
      end do
     
   end subroutine build_union_grid
     
   !==============================================================================
   ! subroutine: build_coefficients
   !
   ! Purpose:
   !   Fit temperature-expansion coefficients at each energy-grid point.
   !  
   ! Parametr IN:
   !   cur_netron_d     -   Neutron state immediately before absorption.
   !   Egrid            -   Energy points at which coefficients are generated.
   !   Tfit             -   Temperatures used for fitting.
   !   nreac            -   Number of fitted reaction channels.
   !   order            -   Expansion order; the coefficient count is 2*order+1.
   !   E                -   Original energy grid supplied to Doppler broadening.
   !   sigma            -   Cross-section vector supplied to Doppler broadening
   !
   ! Parametr OUT:  
   !   coeff            -   Fitted coefficients indexed by coefficient/reaction/energy.
   !   max_error        -   Maximum relative fit error at each reaction/energy point.
   !
   !==============================================================================

   subroutine build_coefficients(Egrid,Tfit,nreac,order,coeff,max_error,E,sigma)
      real, intent(in) :: Egrid(:),Tfit(:), E(:),sigma(:)
      integer, intent(in) :: nreac,order
      real, allocatable, intent(out) :: coeff(:,:,:),max_error(:,:)
   
      integer :: g,j,r,ncoef                     ! Energy, temperature, reaction indices and coefficient count.
      real, allocatable :: sigma_T(:),c(:)        ! Temperature samples and fitted coefficient vector.

      real :: alfa
   
      ncoef=2*order+1
      allocate(coeff(ncoef,nreac,size(Egrid)))
      allocate(max_error(nreac,size(Egrid)))
      allocate(sigma_T(size(Tfit)),c(ncoef))
   
      do g=1,size(Egrid)
         print*,g, " - ", size(Egrid)
         do r=1,nreac
            do j=1,size(Tfit)
               !(E1,E_old,sigma,alfa,n)
               alfa = AWR/(k*Tfit(j))
               sigma_T(j)=doplerBroadr(Egrid(g),E,sigma,alfa,size(E))
            end do
            call fit_one_energy(Tfit,sigma_T,order,c,max_error(r,g))
            coeff(:,r,g)=c
         end do
      end do
   end subroutine build_coefficients
     
   !==============================================================================
   ! subroutine: fit_one_energy
   !
   ! Purpose:
   !   Solve one temperature-fit problem at a fixed energy point.
   !
   ! Parametr IN:
   !   Tfit             -   Temperatures used for fitting.
   !   sigma_T          -   Cross-section vector supplied to Doppler broadening at temperature T
   !   order            -   Expansion order; the coefficient count is 2*order+1.
   !
   ! Parametr OUT:
   !  coeff             -   Resalting coeficient for OTF
   !  max_error         -   Max fluctuation from all modeling coeficient
   !
   !==============================================================================

   subroutine fit_one_energy(Tfit,sigma_T,order,coeff,max_error)
      real, intent(in) :: Tfit(:),sigma_T(:)
      integer, intent(in) :: order
      real, intent(out) :: coeff(2*order+1),max_error
   
      integer :: m,ncoef,rank,info,j          ! Fit dimensions, numerical rank/status, and loop index.
      real :: fit_value                       ! Expansion value used to evaluate fit error.
      real, allocatable :: A(:,:),b(:)        ! Least-squares design matrix and right-hand side.
   
      m=size(Tfit)
      ncoef=2*order+1
      if(m<ncoef) error stop "Too few fit temperatures"
   
      allocate(A(m,ncoef),b(m))
      call fill_matrix(Tfit,order,A)
      b=sigma_T
   
      call least_squares_svd(A,b,coeff,rank,info)
      if(info/=0) error stop "Internal Jacobi SVD did not converge"
   
      max_error=0.0
      do j=1,m
         fit_value=eval_expansion(Tfit(j),minval(Tfit),maxval(Tfit),order,coeff)
         max_error=max(max_error,abs(fit_value-sigma_T(j))/max(abs(sigma_T(j)),1.0e-12))
      end do
   end subroutine fit_one_energy
     
   !==============================================================================
   ! subroutine: fill_matrix
   !
   ! Purpose:
   !   Construct the temperature-basis design matrix used by the OTF fit.
   ! Parametr IN:
   !   Tfit             -   Temperatures used for fitting.
   !   order            -   Expansion order; the coefficient count is 2*order+1.
   !
   ! Parametr OUT:
   !   A                -   Matrix of temperature for different superation and velue
   !
   !==============================================================================

   subroutine fill_matrix(Tfit,order,A)
      real, intent(in) :: Tfit(:)
      integer, intent(in) :: order
      real, intent(out) :: A(:,:)
      integer :: i,j
      real :: q,Tmin,Tmax
   
      Tmin=minval(Tfit)
      Tmax=maxval(Tfit)
      A=0.0
   
      do j=1,size(Tfit)
         q=scaled_temperature(Tfit(j),Tmin,Tmax)
         A(j,1)=1.0
         do i=1,order
            A(j,1+i)=q**(-0.5*real(i))
            A(j,1+order+i)=q**(0.5*real(i))
         end do
      end do
   end subroutine fill_matrix
     
   !==============================================================================
   ! subroutine: least_squares_svd
   !
   ! Purpose:
   !   Solve a least-squares system with the internal one-sided Jacobi SVD.
   !   The matrix columns are scaled before cyclic Jacobi rotations; unresolved
   !   singular directions are discarded using an epsilon-based cutoff.
   !
   !==============================================================================

   subroutine least_squares_svd(Ain,bin,x,rank,info)
      ! Self-contained least-squares solver.  No LAPACK is used.
      ! One-sided cyclic Jacobi SVD is applied to a column-scaled matrix.
      real, intent(in) :: Ain(:,:),bin(:)
      real, intent(out) :: x(size(Ain,2))
      integer, intent(out) :: rank,info
   
      integer :: m,n,p,q,j,sweep
      integer, parameter :: max_sweeps=200
      real :: app,aqq,apq,tau,t,c,s
      real :: tol,threshold,smax,cutoff,projection
      real, allocatable :: B(:,:),V(:,:),scale(:),sv(:),z(:)
      real, allocatable :: bp(:),bq(:),vp(:),vq(:)
      logical :: changed
   
      m=size(Ain,1)
      n=size(Ain,2)
      x=0.0
      rank=0
      info=0
   
      if(size(bin)/=m .or. m<n) then
         info=-1
         return
      end if
   
      allocate(B(m,n),V(n,n),scale(n),sv(n),z(n))
      allocate(bp(m),bq(m),vp(n),vq(n))
   
      do j=1,n
         scale(j)=sqrt(dot_product(Ain(:,j),Ain(:,j)))
         if(scale(j)<=tiny(1.0)) then
            info=-2
            return
         end if
         B(:,j)=Ain(:,j)/scale(j)
      end do
   
      V=0.0
      do j=1,n
         V(j,j)=1.0
      end do
   
      tol=100.0*epsilon(1.0)
   
      do sweep=1,max_sweeps
         changed=.false.
   
         do p=1,n-1
            do q=p+1,n
               app=dot_product(B(:,p),B(:,p))
               aqq=dot_product(B(:,q),B(:,q))
               apq=dot_product(B(:,p),B(:,q))
   
               if(app<=tiny(1.0) .or. aqq<=tiny(1.0)) cycle
   
               threshold=tol*sqrt(app*aqq)
               if(abs(apq)<=threshold) cycle
   
               tau=(aqq-app)/(2.0*apq)
               if(tau>=0.0) then
                  t=1.0/(tau+sqrt(1.0+tau*tau))
               else
                  t=-1.0/(-tau+sqrt(1.0+tau*tau))
               end if
   
               c=1.0/sqrt(1.0+t*t)
               s=c*t
   
               bp=B(:,p)
               bq=B(:,q)
               B(:,p)=c*bp-s*bq
               B(:,q)=s*bp+c*bq
   
               vp=V(:,p)
               vq=V(:,q)
               V(:,p)=c*vp-s*vq
               V(:,q)=s*vp+c*vq
   
               changed=.true.
            end do
         end do
   
         if(.not.changed) exit
      end do
   
      if(changed) then
         info=1
         return
      end if
   
      do j=1,n
         sv(j)=sqrt(dot_product(B(:,j),B(:,j)))
      end do
   
      smax=maxval(sv)
      if(smax<=tiny(1.0)) then
         info=-3
         return
      end if
   
      ! Discard numerically unresolved singular directions.
      cutoff=epsilon(1.0)*real(max(m,n))*smax
   
      z=0.0
      do j=1,n
         if(sv(j)>cutoff) then
            rank=rank+1
            projection=dot_product(B(:,j),bin)/(sv(j)*sv(j))
            z=z+projection*V(:,j)
         end if
      end do
   
      if(rank==0) then
         info=-4
         return
      end if
   
      do j=1,n
         x(j)=z(j)/scale(j)
      end do
   end subroutine least_squares_svd
   
   !==============================================================================
   ! function: eval_expansion
   !
   ! Purpose:
   !   Evaluate the fitted OTF temperature expansion at temperature T.
   !
   ! Parametr IN:
   !   T                -   choosing temperature
   !   t_min            -    Minimum temperature [K].
   !   t_max            -    Maximum temperature [K].
   !   order            -   Expansion order; the coefficient count is 2*order+1.
   !   coeff            -   Fitted coefficients indexed by coefficient/reaction/energy.
   !
   ! Parametr OUT:  
   !   sigma            -   Cross-section vector culculated with OTF
   !
   !==============================================================================

   function eval_expansion(T,Tmin,Tmax,order,coeff) result(sigma)
      real, intent(in) :: T,Tmin,Tmax
      integer, intent(in) :: order
      real, intent(in) :: coeff(2*order+1)
      real :: sigma,q
      integer :: i
     
      q=scaled_temperature(T,Tmin,Tmax)
      sigma=coeff(1)
      do i=1,order
         sigma=sigma+coeff(1+i)*q**(-0.5*real(i)) &
                     +coeff(1+order+i)*q**(0.5*real(i))
      end do
   end function eval_expansion
    
   !==============================================================================
   ! function: scaled_temperature
   !
   ! Purpose:
   !   Map physical temperature to the legacy dimensionless fit coordinate q.
   !
   ! Parametr IN:
   !   T        -    choosing temperature
   !   t_min    -    Minimum temperature [K].
   !   t_max    -    Maximum temperature [K].
   !
   ! Parametr OUT:
   !   q        -   scaled temperature
   !
   !==============================================================================

   function scaled_temperature(T,Tmin,Tmax) result(q)
      real, intent(in) :: T,Tmin,Tmax
      real :: q,Toff
      Toff=(Tmin-Tmax)/50.0
      q=(T-Tmin-Toff)/(Tmax-Tmin-Toff)
   end function scaled_temperature
    
   !==============================================================================
   ! subroutine: insert_point
   !
   ! Purpose:
   !   Insert one value into an allocatable real vector at position pos.
   !
   ! Parametr IN:
   !   pos           -   position in energy mass
   !   value         -   value of new energy
   !
   ! Parametr IN/OUT:
   !   x             -   new mass of energy
   !
   !==============================================================================

   subroutine insert_point(x,pos,value)
      real, allocatable, intent(inout) :: x(:)
      integer, intent(in) :: pos
      real, intent(in) :: value
      real, allocatable :: tmp(:)
      integer :: n
   
      n = size(x)
   
      allocate(tmp(n+1))
   
      if (pos > 1) then
         tmp(1:pos-1) = x(1:pos-1)
      end if
   
      tmp(pos) = value
   
      if (pos <= n) then
         tmp(pos+1:n+1) = x(pos:n)
      end if
   
      call move_alloc(tmp,x)
   end subroutine insert_point
     
end module Dopler_MCNP_OTF