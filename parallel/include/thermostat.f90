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
            integer :: rank, size, chunk_size, start, end, ierr
            real(8) :: sigma, rnumber


            call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
            call mpi_comm_size(MPI_COMM_WORLD, size, ierr)

            ! Calculate the range of particles each process will handle.
            chunk_size = part_num / size
            start = rank * chunk_size + 1

            ! Ensure the last process handles any remaining particles.
            if (rank == size - 1) then
                end = part_num
            else
                end = (rank + 1) * chunk_size
            end if

            sigma = sqrt(temperature)

            do i = start, end
                call random_number(rnumber)

                if (rnumber < collision_frequence) then
                    velocities(i, 1) = sigma * random_gaussian()
                    velocities(i, 2) = sigma * random_gaussian()
                    velocities(i, 3) = sigma * random_gaussian()
                end if
            end do

            ! Gather results back to the root process (rank 0)
            if (rank == 0) then
                ! Gather the results into the root process
                call mpi_gather( &
                    velocities(start:end, :), end-start+1, MPI_REAL8, &
                    velocities, end-start+1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr &
                )
            else
                ! Other processes send their results to rank 0
                call mpi_gather( &
                    velocities(start:end, :), end-start+1, MPI_REAL8, &
                    velocities, end-start+1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr &
                )
            end if
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
