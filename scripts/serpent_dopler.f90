module serpent_dopler
    
    use data_type
    use operation_with_data
    implicit none
    
    real, parameter :: KB_EV = 8.617333262e-5
    real, parameter :: PI_R  = 3.14159265358979323846

contains
    real function tms_g(E, dT, A) result(g)
        real, intent(in) :: E, dT, A
        real :: kappa, x, x2
        if (dT <= 0.0) then
            g = 1.0
            return
        end if
        kappa = sqrt(A/(KB_EV*dT))
        x = kappa*sqrt(E)
        x2 = x*x

        if (x > 12.0) then
            ! erf(x) ~ 1 and exp(-x^2) ~ 0
            g = 1.0 + 0.5/x2
        else
            g = (1.0 + 0.5/x2)*erf(x) + exp(-x2)/(sqrt(PI_R)*x)
        end if
    end function tms_g

    subroutine tms_energy_bounds(E, dT, A, Emin, Emax)
        real, intent(in)  :: E, dT, A
        real, intent(out) :: Emin, Emax
        real :: kappa, rootE, droot

        kappa = sqrt(A/(KB_EV*dT))
        rootE = sqrt(max(E, 0.0))
        droot = 4.0/kappa

        Emin = max(0.0, rootE - droot)**2
        Emax = (rootE + droot)**2
    end subroutine tms_energy_bounds
    
    real function max_total_micro_xs(nuclide, Emin, Emax,tem,env,index) result(sigmax)
        type(enviroment), intent(in) :: env
        real, intent(in) ::tem
        integer, intent(in):: index
        type(nuclear_data), intent(in) :: nuclide
        real, intent(in) :: Emin, Emax

        integer :: j, n
        real :: Ej

        sigmax = max(found_cross_section_from_energy(Emin, tem, env, nuclide, index), &
                     found_cross_section_from_energy(Emax, tem, env, nuclide, index))

        n = size(nuclide%energy_point_in_table)
        do j = 1, n
            Ej = nuclide%energy_point_in_table(j)
            if (Ej >= Emin .and. Ej <= Emax) then
                sigmax = max(sigmax, &
                    max(0.0, nuclide%cross_data(index)%cross_section_point_in_table(j,1)))!!!!!!!
            end if
        end do
    end function max_total_micro_xs


    real function tms_nuclide_majorant(nuclide, E, Tbase, Tmax,env,index) result(Smaj)
        type(nuclear_data), intent(in) :: nuclide
        real, intent(in) :: E, Tbase, Tmax
        type(enviroment), intent(in) :: env
        integer, intent(in):: index
        real :: dT, Emin, Emax, sigmax, g

        dT = max(0.0, Tmax - Tbase)

        call tms_energy_bounds(E, dT, nuclide%mass_of_nuclear, Emin, Emax)

        sigmax = max_total_micro_xs(nuclide, Emin, Emax, Tbase,env,index)
        g = tms_g(E, dT, nuclide%mass_of_nuclear)

        Smaj = g * sigmax
    end function tms_nuclide_majorant

    subroutine tms_sample_target_energy(E, dT, A, Eprime)
        real, intent(in)  :: E, dT, A
        real, intent(out) :: Eprime

        real :: kappa, x, y, mean_y, mix_mb
        real :: mu, vrel, denom, u, u1, u2
        real :: gx, gy, gz, sigma_comp

        if (dT <= 0.0) then
                Eprime = E
                return
        end if
        kappa = sqrt(A/(KB_EV*dT))
        x = sqrt(E)

        ! Mean Maxwell speed in the same sqrt(eV) velocity unit.
        mean_y = 2.0/(sqrt(PI_R)*kappa)

        ! Mixture weight of the ordinary Maxwell speed distribution.
        mix_mb = x/(x + mean_y)

        do
            call random_number(u)

            if (u < mix_mb) then
                ! Maxwell speed: magnitude of 3 Gaussian components.
                sigma_comp = 1.0/(sqrt(2.0)*kappa)
                gx = gaussian01()*sigma_comp
                gy = gaussian01()*sigma_comp
                gz = gaussian01()*sigma_comp
                y = sqrt(gx*gx + gy*gy + gz*gz)
            else
                ! Speed-biased Maxwell:
                ! z = kappa^2*y^2 ~ Gamma(shape=2, scale=1)
                call positive_uniform(u1)
                call positive_uniform(u2)
                y = sqrt(-log(u1*u2))/kappa
            end if

            call random_number(u)
            mu = 2.0*u - 1.0

            vrel = sqrt(max(0.0, x*x + y*y - 2.0*x*y*mu))
            denom = x + y

            if (denom <= tiny(1.0)) then
                cycle
            end if

            call random_number(u)
            if (u <= min(1.0, vrel/denom)) exit
        end do

        Eprime = vrel*vrel
    end subroutine tms_sample_target_energy
    
    real function gaussian01() result(z)
        real :: u1, u2
        call positive_uniform(u1)
        call random_number(u2)
        z = sqrt(-2.0*log(u1))*cos(2.0*PI_R*u2)
    end function gaussian01

    subroutine positive_uniform(u)
        real, intent(out) :: u
        call random_number(u)
        if (u <= tiny(1.0)) u = tiny(1.0)
    end subroutine positive_uniform

end module serpent_dopler