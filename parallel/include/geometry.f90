!
! geometry.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Ricard Rodriguez
!
! Contains all the subroutines dependent or related to the system geometry.
!

module geometry
    use global_vars

    implicit none

    contains
        ! Apply PBC conditions to the positions array.
        ! This subroutine must be run every time the positions of the particles are
        ! updated.
        subroutine apply_pbc(positions, counts, displs)
            use mpi

            implicit none

            real(8), allocatable, intent(inout) :: positions(:, :)
            integer, allocatable, intent(in) :: counts(:), displs(:) ! MPI arguments.

            integer :: i, j


            ! Apply PBC to the assigned chunk of positions
            do i = start_part, end_part
                do j = 1, 3
                    positions(i, j) = pbc(positions(i, j), system_size)
                end do
            end do

            call mpi_barrier(MPI_COMM_WORLD, ierr)

            ! Gather the results into the root process
            do i = 1, 3
                call mpi_allgatherv( &
                    positions(start_part : end_part, i), counts(rank), MPI_REAL8, &
                    positions(:, i), counts, displs, MPI_REAL8, MPI_COMM_WORLD, ierr &
                )
            end do
        end subroutine apply_pbc

        ! Apply PBC to a certain distance, given the box size.
        function pbc(distance, box_size)
            implicit none

            real(8) :: distance, box_size, pbc

            if (distance > box_size/2) then
                distance = distance - box_size
            else if (distance < -box_size/2) then
                distance = distance + box_size
            end if

            pbc = distance
        end function pbc
end module geometry
