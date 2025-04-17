!
! integrators.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Oriol Miro
!
! Integrators for Molecular Dynamics.
! A 3-dimensional cubic system is assumed.
!

module integrators
    use global_vars
    use lj_forces
    use geometry
    use mpi

    implicit none

    contains

        subroutine velocity_verlet(positions, velocities, lj_potential, verlet_list, n_neighbors, counts, displs)
            implicit none

            real(8), allocatable, intent(inout) :: positions(:, :), velocities(:, :)
            real(8), intent(out) :: lj_potential
            integer, allocatable, intent(in) :: verlet_list(:, :), n_neighbors(:)
            integer, allocatable, intent(in) :: counts(:), displs(:) ! MPI arguments.

            integer :: i
            real(8), allocatable :: forces(:, :)


            call compute_forces(positions, forces, lj_potential, verlet_list, n_neighbors)
            do i = start_part, end_part
                positions(i,1) = positions(i,1) + velocities(i,1)*timestep + 0.5 * forces(i,1)*timestep*timestep
                positions(i,2) = positions(i,2) + velocities(i,2)*timestep + 0.5 * forces(i,2)*timestep*timestep
                positions(i,3) = positions(i,3) + velocities(i,3)*timestep + 0.5 * forces(i,3)*timestep*timestep

                velocities(i,1) = velocities(i,1) + 0.5 * forces(i,1)*timestep
                velocities(i,2) = velocities(i,2) + 0.5 * forces(i,2)*timestep
                velocities(i,3) = velocities(i,3) + 0.5 * forces(i,3)*timestep
            end do

            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            do i = 1, 3
                call MPI_Allgatherv(positions(start_part : end_part, i), counts(rank), MPI_REAL8, &
                                    positions(:,i), counts, displs, MPI_REAL8, MPI_COMM_WORLD, ierr)
            end do

            if (rank == 0) then
                call apply_pbc(positions, counts, displs)
            end if

            call compute_forces(positions, forces, lj_potential, verlet_list, n_neighbors)
            velocities = velocities + 0.5 * forces*timestep

            ! The simulation works without this gathering, but it's better if all the velocities are only divided inside this
            ! subroutine. to work better in the final versions of the code.
            do i = 1, 3
                call MPI_Allgatherv(velocities(start_part : end_part, i), counts(rank), MPI_REAL8, &
                                    velocities(:, i), counts, displs, MPI_REAL8, MPI_COMM_WORLD, ierr)
            end do
        end subroutine velocity_verlet
end module integrators
