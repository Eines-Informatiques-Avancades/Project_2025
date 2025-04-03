!
! thermostat.f90
! Molecular Dynamics simulation of 3D Van der Waals gas
! Oriol Miro
!
! Andersen thermostat for Molecular Dynamics.
! A 3-dimensional cubic system is assumed.
!

module thermostat
    use global_vars

    implicit none

    contains
        ! Andersen thermostat. Updates the velocity magnitude and direction of a random
        ! number of particles to one that follows the Maxwell-Boltzmann distribution.
        subroutine andersen_thermostat(velocities)
            implicit none

            real, allocatable, intent(inout) :: velocities(:)

            real :: sigma, rnumber
            integer :: i

            sigma = sqrt(temperature)

            do i = 1, part_num
                call random_number(rnumber)

                if (rnumber < collision_frequence) then
                    velocities(i) = sigma * random_gaussian()
                end if
            end do
        end subroutine andersen_thermostat

        ! Function that returns a random number following the gaussian distribution.
        function random_gaussian() result(r)
            implicit none

            real :: r, u1, u2

            call random_number(u1)
            call random_number(u2)

            r = real(sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * 3.141592653589793d0 * u2))
        end function random_gaussian
end module thermostat
