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
        subroutine verlet(dt, positions, positions_old, velocities, lj_potential)
            implicit none

            real(8), intent(in) :: dt
            real(8), allocatable, intent(inout) :: positions(:, :), positions_old(:, :), velocities(:, :)
            real(8), intent(out) :: lj_potential
            
            integer :: i
            real(8), allocatable :: positions_aux(:, :), forces(:, :)

            allocate(positions_aux(part_num, 3))

            call compute_forces(positions, forces, lj_potential)
            positions_aux = positions
            positions = 2*positions - positions_old + forces*dt*dt
            positions_old = positions_aux
            velocities = (positions - positions_old)/dt

            do i = 1, 3
                call apply_pbc(positions(:,3))
            end do
        end subroutine verlet

        subroutine velocity_verlet(dt, positions, velocities, lj_potential, counts, displs, rank, nproc)
            implicit none

            real(8), intent(in) :: dt
            real(8), allocatable, intent(inout) :: positions(:, :), velocities(:, :)
            real(8), intent(out) :: lj_potential

            integer, intent(in) :: nproc, rank
            integer, intent(in) :: counts(0:nproc-1), displs(0:nproc-1)
            integer :: i, local_n, ierror
            real(8), allocatable :: positions_local(:, :), velocities_local(:, :), forces_local(:, :), forces(:,:)



            ! Número total de partícules
            local_n = counts(rank)

            allocate(positions_local(local_n, 3), velocities_local(local_n, 3), forces_local(local_n, 3))
            print*, "counts", lbound(counts), Ubound(counts)
            ! Repartir datos
            call compute_forces(positions, forces, lj_potential)
            forces_local = forces(displs(rank) + 1:displs(rank)+ local_n, :)

            do i = 1,3
                call MPI_Scatterv(positions(:, i), counts, displs, MPI_REAL8, &
                                positions_local(:, i), local_n, MPI_REAL8, 0, MPI_COMM_WORLD, ierror)

                call MPI_Scatterv(velocities(:, i), counts, displs, MPI_REAL8, &
                                velocities_local(:, i), local_n, MPI_REAL8, 0, MPI_COMM_WORLD, ierror)
            end do
        ! ==== VELOCITY VERLET LOCAL ====


            positions_local = positions_local + velocities_local*dt + 0.5 * forces_local*dt*dt
            velocities_local = velocities_local + 0.5 * forces_local*dt
            do i = 1, 3
                call apply_pbc(positions_local(:,3))
            end do
            call MPI_Allgatherv(positions_local(:, i), local_n, MPI_REAL8, &
            positions, counts, displs, MPI_REAL8, MPI_COMM_WORLD, ierror)

            call compute_forces(positions, forces, lj_potential)
            forces_local = forces(displs(rank)+1:displs(rank)+ local_n, :)

            velocities_local = velocities_local + 0.5 * forces_local*dt
            call MPI_Barrier(MPI_COMM_WORLD, ierror)
            ! ==== RECOLECCIÓ GATHERV ====
            do i = 1,3
 
                call MPI_Allgatherv(positions_local(:, i), local_n, MPI_REAL8, &
                                positions, counts, displs, MPI_REAL8, MPI_COMM_WORLD, ierror)

                call MPI_Allgatherv(velocities_local(:, i), local_n, MPI_REAL8, &
                                velocities, counts, displs, MPI_REAL8, MPI_COMM_WORLD, ierror)
            end do


            deallocate(positions_local, velocities_local, forces_local)
        end subroutine velocity_verlet


        subroutine euler(dt, positions, velocities, lj_potential)
            implicit none

            real(8), intent(in) :: dt
            real(8), allocatable, intent(inout) :: positions(:, :), velocities(:, :)
            real(8), intent(out) :: lj_potential

            real(8), allocatable :: forces(:, :)

            integer :: i

            call compute_forces(positions, forces, lj_potential)
            positions = positions + velocities * dt + 0.5 * forces*dt*dt
            velocities = velocities + forces*dt

            do i = 1, 3
                call apply_pbc(positions(:,3))
            end do
        end subroutine euler
end module integrators
