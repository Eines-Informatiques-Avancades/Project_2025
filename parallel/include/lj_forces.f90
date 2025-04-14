!
! lj_forces.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Alejandro Díaz
!
! Compute forces derived from a Lennard-Jones potential using Verlet lists.
!

module lj_forces
    use global_vars
    use geometry

    implicit none

    contains
        ! Compute Verlet lists, one for each particle in the system.
        ! The Verlet list for one particle contains a list of neighbours it interactis with.
        ! This subroutine must not be run at each step, but every 5 to 15 timesteps.
        subroutine compute_verlet_list(positions, verlet_list, n_neighbors)
            implicit none

            real(8), intent(in) :: positions(:, :)
            integer, allocatable, intent(out) :: verlet_list(:, :)
            integer, allocatable, intent(out) :: n_neighbors(:)
            integer :: i, j
            real(8) :: dx, dy, dz, r2, cutoff_verlet

            ! Must be slightly bigger than the Lennard-Jones cutoff.
            cutoff_verlet = cutoff * 1.2

            if (allocated(n_neighbors)) then
                deallocate(n_neighbors)
            endif
            if (allocated(verlet_list)) then
                deallocate(verlet_list)
            endif

            ! Set a long enough list range
            allocate(n_neighbors(part_num))
            n_neighbors = 0

            allocate(verlet_list(part_num, part_num))
            verlet_list = 0

            ! Compute pairs (i, j) with j > i.
            do i = 1, part_num - 1
                do j = i + 1, part_num
                    dx = positions(i, 1) - positions(j, 1)
                    dy = positions(i, 2) - positions(j, 2)
                    dz = positions(i, 3) - positions(j, 3)

                    ! Apply pbc
                    dx = pbc(dx, system_size)
                    dy = pbc(dy, system_size)
                    dz = pbc(dz, system_size)

                    r2 = dx*dx + dy*dy + dz*dz

                    ! Allocate interactions within the cutoff_verlet
                    if (r2 <= cutoff_verlet**2) then
                        n_neighbors(i) = n_neighbors(i) + 1
                        verlet_list(i, n_neighbors(i)) = j
                    endif
                end do
            end do
        end subroutine compute_verlet_list

        subroutine compute_forces(positions, forces, lj_potential, verlet_list, n_neighbors)
            use mpi
            use global_vars
            use geometry

            implicit none

            real(8), allocatable, intent(in)  :: positions(:, :)
            real(8), allocatable, intent(out) :: forces(:, :)
            real(8), intent(out)              :: lj_potential
            integer, intent(in)               :: verlet_list(:, :)
            integer, intent(in)               :: n_neighbors(:)

            integer :: i, k, j, h, ierr, rank, nprocs
            integer :: tot_inter, local_target, n_inter_acc
            integer, allocatable :: assign_start(:), assign_end(:)
            integer :: i_start, i_end, current_proc, target
            real(8) :: r, f, local_lj_potential
            real(8) :: r_vec(3)
            real(8), allocatable :: local_forces(:, :)

            call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
            call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)

            tot_inter = 0
            do i = 1, part_num - 1
                tot_inter = tot_inter + n_neighbors(i)
            end do

            local_target = tot_inter / nprocs

            ! Divide tasks between processors (last particle does not have a list)
            allocate(assign_start(nprocs))
            allocate(assign_end(nprocs))
            current_proc = 0
            n_inter_acc = 0
            assign_start(current_proc+1) = 1

            ! Task assigned based on the total number of interactions
            do i = 1, part_num - 1
                n_inter_acc = n_inter_acc + n_neighbors(i)

                if (current_proc < mod(tot_inter, nprocs)) then
                    target = local_target + 1
                else
                    target = local_target
                endif

                if (n_inter_acc >= target .and. current_proc < nprocs - 1) then
                    assign_end(current_proc+1) = i
                    current_proc = current_proc + 1
                    assign_start(current_proc+1) = i + 1
                    n_inter_acc = 0
                end if
            end do
            assign_end(nprocs) = part_num - 1

            i_start = assign_start(rank + 1)
            i_end   = assign_end(rank + 1)
            ! Assign local matrix
            allocate(local_forces(part_num, 3))
            local_forces = 0.0
            local_lj_potential = 0.0

            ! Each processors just works with their indexes
            do i = i_start, i_end
                do k = 1, n_neighbors(i)
                    j = verlet_list(i, k)

                    ! Compute distance
                    r_vec(:) = positions(i, :) - positions(j, :)

                    ! Apply PBC
                    do h = 1, 3
                        r_vec(h) = pbc(r_vec(h), system_size)
                    end do

                    r = sqrt(dot_product(r_vec, r_vec))

                    if (r < cutoff .and. r > 1.0e-8) then
                        ! Local Lennard-Jones potential
                        local_lj_potential = local_lj_potential + 4.0 * (1.0 / (r**12) - 1.0 / (r**6))

                        ! Compute forces
                        f = 48.0 / (r**14) - 24.0 / (r**8)

                        ! Update forces
                        local_forces(i, :) = local_forces(i, :) + f * r_vec(:)
                        local_forces(j, :) = local_forces(j, :) - f * r_vec(:)
                    end if
                end do
            end do

            ! Combine all tasks
            if (allocated(forces)) then
                deallocate(forces)
            endif
            allocate(forces(part_num, 3))
            forces = 0.0
            call mpi_barrier(MPI_COMM_WORLD, ierr)
            call mpi_allreduce(local_forces, forces, part_num*3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            call mpi_allreduce(local_lj_potential, lj_potential, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

            deallocate( &
                local_forces, assign_start, assign_end &
            )
        end subroutine compute_forces
end module lj_forces
