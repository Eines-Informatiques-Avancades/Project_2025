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
        subroutine apply_pbc(positions)
            use mpi

            implicit none

            real(8), allocatable, intent(inout) :: positions(:, :)

            integer :: i, j
            integer, allocatable :: counts(:), displs(:)

            allocate(counts(0:nproc - 1))
            allocate(displs(0:nproc - 1))

            do i = 0, nproc - 1
                if (i < mod(part_num, nproc)) then
                    counts(i) = (part_num / nproc + 1)
                else
                    counts(i) = (part_num / nproc)
                end if
            end do

            displs(0) = 0
            do i = 1, nproc - 1
                displs(i) = displs(i - 1) + counts(i - 1)
            end do

            ! Apply PBC to the assigned chunk of positions
            do i = displs(rank) + 1, displs(rank) + counts(rank)
                do j = 1, 3
                    positions(i, j) = pbc(positions(i, j), system_size)
                end do
            end do

            print *, 'Rank: ', rank, 'start: ', start, 'd+1: ', displs(rank)+1, 'end: ', end, 'd+c', displs(rank)+counts(rank)

            call mpi_barrier(MPI_COMM_WORLD, ierr)

            ! Gather the results into the root process
            do i = 1, 3
                call mpi_allgatherv( &
                    positions(displs(rank) + 1 : displs(rank) + counts(rank), i), counts(rank), MPI_REAL8, &
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
