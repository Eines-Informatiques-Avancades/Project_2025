!
! post_trajectory_analysis
! Molecular Dynamics Simulation of a Van der Waals Gas
! Haoyu Huang, Ricard Rodriguez
!

module post_trajectory_analysis
    use global_vars

    implicit none

    contains
        ! Read data from the positions file and store it into arrays.
        ! Particle positions in time are stored in columns: time x y z
        ! Each row refers to a particle in an specific time.
    subroutine read_trajectory(positions_file, x, y, z, time)
        use mpi
        
        implicit none

        character(50), intent(in) :: positions_file
        real(8), allocatable, intent(inout) :: x(:, :), y(:, :), z(:, :), time(:)

        integer :: part, step, ios, n
        real(8) :: t


        if (rank == 0) then
            print *, 'Reading trajectory data from ', positions_file
            open(4, file = positions_file, status = 'old', action = 'read')

            ! n is the total number of particles and t is the time
            do step = 1, step_num
                read(4, *, iostat = ios) n
                read(4, *, iostat = ios) t
                do part = 1, part_num
                    read(4, *, iostat = ios) atom_type, x(part, step), y(part, step), z(part, step)

                    ! Store the different values of t (once for each different time).
                    if (part == 1) then
                        time(step) = t
                    end if
                end do
            end do

            close(4)
        end if

        ! 0 (master) broadcast arrays of information to other processes
        call MPI_Bcast(x, part_num * step_num, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)  ! (x array, (xarray component1) * (xarray component2), type, 'the thread rank that bcast', mpi_comm_world, ierr )
        call MPI_Bcast(y, part_num * step_num, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)  ! same with y, z and time arrays
        call MPI_Bcast(z, part_num * step_num, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(time, step_num, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    end subroutine read_trajectory

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                   !
! TEST DONT NEED TO PARALLELIZE                                     !
!                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Test to proove that the program has access to all the stored xyz information.
    subroutine test_access_xyz_data(x, y, z, time)
            implicit none

            real(8), allocatable, intent(in) :: x(:, :), y(:, :), z(:, :), time(:)

            integer :: i, j

            print *, 'Processing stored coordinates...'

            do j = 1, step_num
                print *, 'time:', time(j)

                do i = 1, part_num
                    print *, '(', x(i, j), y(i, j), z(i, j), ')'
                end do
            end do
    end subroutine test_access_xyz_data

        ! Compute RDF using the stored data.
        ! Must be executed after read_trajectory, as it depends on the arrays it
        ! creates.
    subroutine compute_rdf(x, y, z, rdf_file)
            use mpi
            
            implicit none

            real(8), allocatable, intent(in) :: x(:, :), y(:, :), z(:, :)
            character(50), intent(in) :: rdf_file

            integer :: i, j, k, time_index
            real(8) :: maximum_radius, volume, density
            integer :: bins, start, fin, frame_chunk
            real(8), allocatable :: h(:), rdf(:), r_values(:), h_total(:)
            real(8) :: r, r_sq, dx, dy, dz, dv, r_lo, r_hi, const, nid
            integer :: bin_index, nproc, ierr, rank

            call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
            call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

            ! Parameters
            maximum_radius = system_size / 2
            bins = int(maximum_radius /timestep)
            volume = system_size**3
            density = part_num / volume

            ! Allocate arrays
            allocate(h(bins), rdf(bins), r_values(bins), h_total(bins))

            h = 0.0
            h_total = 0.0
            frame_chunk = step_num / nproc       ! divide the work, chunk is the amount of work for each thread.
            start = rank * frame_chunk + 1      ! start and fin are index for the mpi loop, this line ensures beginning index of each thread.
            fin = start + frame_chunk - 1       ! this line ensures the ending index.
                                                ! i.e  100 steps in 4 nproc, meaning we have 4 loops for each thread,
                                                ! |rank=0, start=0*25+1=1, fin= 1+25-1=25| |rank=1, start=1*25+1=26, fin= 26+25-1=50|, ...

            if (rank == nproc - 1) then
                fin = step_num                  ! ensuring for step_num is not divisible by nproc (like 100/7) the ending index is the step_num.
            endif                               ! 100/7 = 14; 14*7 = 98, if this condition were not imposed, the loop is incomplete.

            ! Compute histogram h(k)
            do time_index = start, fin
                do i = 1, part_num - 1
                    do j = i + 1, part_num

                        ! Compute the eucledian distance
                        dx = x(j, time_index) - x(i, time_index)
                        dy = y(j, time_index) - y(i, time_index)
                        dz = z(j, time_index) - z(i, time_index)

                        ! PBC
                        dx = dx - system_size * anint(dx / system_size)
                        dy = dy - system_size * anint(dy / system_size)
                        dz = dz - system_size * anint(dz / system_size)

                        r_sq = dx**2 + dy**2 + dz**2
                        r = sqrt(r_sq)

                        if (r < maximum_radius) then
                            bin_index = floor(r / timestep) + 1

                            ! Boundary check to prevent out-of-bounds access
                            if (bin_index >= 1 .and. bin_index <= bins) then
                                h(bin_index) = h(bin_index) + 2  ! pairwise counting
                            else
                                 print *, 'Warning: bin_index out of range:', bin_index, ' (max bins:', bins, ')'
                            endif
                        endif
                    end do
                end do
            end do

            ! Here the reduction is applied to all local values of h, converting it to h_total at rank 0
            call mpi_reduce(h,h_total,bins, mpi_real8, mpi_sum, 0, mpi_comm_world, ierr)

            ! Now normalize the h stored in h_total by master
            if (rank == 0) then
                const = 4.0 * 3.14159265358979 / 3.0
                do k = 1, bins
                    r_lo = (k - 1) * timestep
                    r_hi = r_lo + timestep
                    dv = const * (r_hi**3 - r_lo**3)  ! Shell volume
                    nid = dv * density
                    rdf(k) = h_total(k) / (part_num * step_num * nid )
                    r_values(k) = (k - 0.5) * timestep  ! Bin center
                end do

            ! Save RDF results to file
                open(12, file = rdf_file, status = 'replace')
                do k = 1, bins
                    write(12, '(F10.5, F15.8)') r_values(k), rdf(k)
                end do
                close(12)

                print *, 'RDF calculation completed and saved to ', rdf_file
            endif

            ! Deallocate memory
            deallocate(h,h_total, rdf, r_values)

    end subroutine compute_rdf

        ! Compute RMSD using the stored data.
        ! Must be executed after read_trajectory, as it depends on the arrays it
        ! creates.
    subroutine compute_rmsd(x, y, z, time, rmsd_file)
        use mpi

        implicit none

        real(8), allocatable, intent(in) :: x(:, :), y(:, :), z(:, :), time(:)
        character(50), intent(in) :: rmsd_file

        integer :: i, j, p, rank, nproc, ierr
        integer :: partial_steps, frame_chunk, start, fin, actual_step
        integer, allocatable :: recvcounts(:), displs(:)
        real(8), allocatable :: total_rmsd(:), partial_rmsd(:)
        real(8) :: dx, dy, dz, sum_sq

        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

        frame_chunk = step_num / nproc       ! Division of work
        start = rank * frame_chunk + 1      ! Each division start indexation
        fin = start + frame_chunk - 1       ! Each division ending indexation
        if (rank == nproc - 1) then          ! Any residual work in the no divisible case gives to the last worker, by forcing the last indexation be the total step num
            fin = step_num
        endif
        partial_steps = fin - start + 1

        allocate(partial_rmsd(partial_steps))

        do j = 1, partial_steps
            actual_step = start + j - 1
            sum_sq = 0.0
            do i = 1, part_num
                dx = x(i, actual_step) - x(i, 1)
                dy = y(i, actual_step) - y(i, 1)
                dz = z(i, actual_step) - z(i, 1)

                sum_sq = sum_sq + (dx**2 + dy**2 + dz**2)
            end do
            partial_rmsd(j) = sqrt(sum_sq / part_num)
        end do

        ! Rank 0 allocates full arrays and prepare displacements
        if (rank == 0) then
            allocate(total_rmsd(step_num))
            allocate(recvcounts(nproc), displs(nproc))

            do p = 0, nproc - 1
                recvcounts(p + 1) = frame_chunk
            end do
            recvcounts(nproc) = step_num - (nproc - 1) * frame_chunk  ! last rank may get more
            displs(1) = 0
            do p = 2, nproc
                displs(p) = displs(p - 1) + recvcounts(p - 1)
            end do
        end if

        ! Use Gatherv to merge all rmsd
        call MPI_Gatherv(partial_rmsd, partial_steps, MPI_REAL8, &
            total_rmsd, recvcounts, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

        if (rank == 0) then
            open(12, file = rmsd_file, status = 'replace')
            do j = 1, step_num
                write(12, *) time(j), total_rmsd(j)
            end do
            close(12)
            print *, 'RMSD calculation completed and saved to ', rmsd_file
            deallocate(total_rmsd, recvcounts, displs)
        end if

        deallocate(partial_rmsd)
    end subroutine compute_rmsd
end module
