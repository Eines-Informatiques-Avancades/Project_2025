!
! thermostat.f90
! Molecular Dynamics simulation of 3D Van der Waals gas
! Ricard Rodríguez
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
            use mpi

            implicit none

            real(8), allocatable, intent(inout) :: velocities(:, :)

            integer :: i
            real(8) :: sigma, rnumber


            sigma = sqrt(temperature)

            do i = start_part, end_part
                call random_number(rnumber)

                if (rnumber < collision_frequence) then
                    velocities(i, 1) = sigma * random_gaussian()
                    velocities(i, 2) = sigma * random_gaussian()
                    velocities(i, 3) = sigma * random_gaussian()
                end if
            end do

            ! Gather the results into the root process
            do i = 1, 3
                call mpi_allgatherv( &
                    velocities(start_part : end_part, i), counts(rank), MPI_REAL8, &
                    velocities(:, i), counts, displs, MPI_REAL8, MPI_COMM_WORLD, ierr &
                )
            end do

        end subroutine andersen_thermostat

        ! Function that returns a random number following the gaussian distribution.
        function random_gaussian() result(r)
            implicit none

            real(8) :: r, u1, u2

            call random_number(u1)
            call random_number(u2)

            r = dble(sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * 3.141592653589793d0 * u2))
        end function random_gaussian
end module thermostat
