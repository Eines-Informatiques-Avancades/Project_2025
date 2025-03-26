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
        real, allocatable, intent(inout) :: x(:, :), y(:, :), z(:, :), time(:)

        integer :: part, step, ios, n, ierr
        integer :: rank, size
        real :: t

        ! Initialize MPI
        !!call MPI_Init(ierr)          ! MOVE THIS TO THE MAIN
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

        if (rank == 0) then
         print *, 'Reading trajectory data from ', positions_file
         open(4, file = positions_file, status = 'old', action = 'read')

         do step = 1, step_num
            do part = 1, part_num
                read(4, *, iostat = ios) n
                read(4, *, iostat = ios) t
                read(4, *, iostat = ios) atom_type, x(part, step), y(part, step), z(part, step)

                if (part == 1) then
                    time(step) = t
                end if
            end do
         end do

        close(4)
        end if

    ! Broadcast arrays of information to other processes
        call MPI_Bcast(x, part_num * step_num, MPI_REAL, 0, MPI_COMM_WORLD, ierr)  ! (x array, (xarray component1) * (xarray component2), type, 'the thread rank that bcast', mpi_comm_world, ierr )
        call MPI_Bcast(y, part_num * step_num, MPI_REAL, 0, MPI_COMM_WORLD, ierr)  ! same with y, z and time arrays
        call MPI_Bcast(z, part_num * step_num, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(time, step_num, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

    
    ! call MPI_Finalize(ierr)     ! MOVE THIS TO THE MAIN
    end subroutine read_trajectory

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                   !
! TEST DONT NEED TO PARALLELIZE                                     !
!                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Test to proove that the program has access to all the stored xyz information.
    subroutine test_access_xyz_data(x, y, z, time)
            implicit none

            real, allocatable, intent(in) :: x(:, :), y(:, :), z(:, :), time(:)

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
            implicit none

            real, allocatable, intent(in) :: x(:, :), y(:, :), z(:, :)
            character(50), intent(in) :: rdf_file

            integer :: i, j, k, time_index
            real(8), parameter :: dr = 0.1
            real(8) :: maximum_radius, volume, density
            integer :: bins
            real(8), allocatable :: h(:), rdf(:), r_values(:), h_total(:)
            real(8) :: r, r_sq, dx, dy, dz, dv, h_id, r_lo, r_hi, const, nid
            integer :: bin_index, size, ierr, rank

            ! Parameters
            maximum_radius = system_size
            bins = int(maximum_radius / dr)
            volume = system_size**3
            density = part_num / volume


            call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
            call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

            ! Allocate arrays
            allocate(h(bins), rdf(bins), r_values(bins), h_total(bins))

            h = 0.0
            h_total = 0.0
            frame_chunk = step_num / size               ! divide the work, chunk is the amount of work for each thread.
            start = rank * frame_chunk + 1              ! start and fin are index for the mpi loop, this line ensures beginning index of each thread .
            fin = start + frame_chunk - 1               ! this line ensures the ending index.
                                                        ! i.e  100 steps in 4 size, meaning we have 4 loops for each thread,  
                                                        !|rank=0, start=0*25+1=1, fin= 1+25-1=25| |rank=1, start=1*25+1=26, fin= 26+25-1=50| , ...

            if (rank == size - 1) then              
                fin = step_num                          ! ensuring for step_num is not divisible by size (like 100/7) the ending index is the step_num.
            endif                                       ! 100//7 = 14   14*7 = 98, if this condition were not imposed, the loop is incomplete.        

            ! Compute histogram h(k)
            do time_index = start, fin
                do i = 1, part_num - 1
                    do j = i + 1, part_num
                        ! Compute minimum image distance between particles with PBC
                        dx = x(j, time_index) - x(i, time_index)
                        dy = y(j, time_index) - y(i, time_index)
                        dz = z(j, time_index) - z(i, time_index)

                        dx = dx - system_size * anint(dx / system_size)
                        dy = dy - system_size * anint(dy / system_size)
                        dz = dz - system_size * anint(dz / system_size)

                        r_sq = dx**2 + dy**2 + dz**2
                        r = sqrt(r_sq)

                        if (r < maximum_radius) then
                            bin_index = floor(r / dr) + 1
                        endif
                            h(bin_index) = h(bin_index) + 2  ! pairwise counting
                    end do
                end do
            end do

            ! Here the reduction is applied to all local values of h, converting it to h_total at rank 0
            call mpi_reduce(h,h_total,bins, mpi_real, mpi_sum, 0, mpi_comm_world, ierr)


            ! Normalize RDF
            ! Now normalize the h stored in h_total 
            if (rank == 0) then 
                const = 4.0 * 3.14159265358979 * density / 3.0
                do k = 1, bins
                    r_lo = (k - 1) * dr
                    r_hi = r_lo + dr
                    dv = const * (r_hi**3 - r_lo**3)  ! Shell volume
                    nid = dv
                    rdf(k) = (h_total(k) / (part_num * step_num)) / nid
                    r_values(k) = (k - 0.5) * dr  ! Bin center
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
            implicit none

            real, allocatable, intent(in) :: x(:, :), y(:, :), z(:, :), time(:)
            character(50), intent(in) :: rmsd_file

            integer :: i, j, size, rank, ierr, partial_steps,frame_chunk, start, fin
            real(4), allocatable :: total_rmsd(:), partial_rmsd(:)
            real(4) :: dx, dy, dz, sum_sq

            call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
            call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

            frame_chunk = step_num / size               ! Divide the work, chunk is the amount of work for each thread.
            start = rank * frame_chunk + 1              ! Start and fin are index for the mpi loop, this line ensures beginning index of each thread .
            fin = start + frame_chunk - 1               ! This line ensures the ending index.


            if (rank == size - 1) then              
                fin = step_num                          ! Ensuring for step_num is not divisible by size (like 100/7) the ending index is the step_num.
            endif                                       ! 100//7 = 14   14*7 = 98, if this condition were not imposed, the loop is incomplete.        

        
            partial_steps = fin - start + 1             ! To divide the work and not spending so much memory, a partial step definition is then needed. 
            allocate(partial_rmsd(partial_steps))

            if (rank == 0) then                         ! When rank 0 can merge all partial rmsd happens the following
                allocate(total_rmsd(step_num))          ! Ensure that the total_rmsd has the total indexation size.
            endif

                                                                   
            do j = 1, partial_steps                             ! For each partial steps.

                actual_step = start + j - 1                     ! i.e.  start = 1 and fin = 5, do j = 1,5 , actual step = 1 + 1 -1 = 1, then do all loop of i with step = 1.
                                                                ! then repeating for others j.
                                                                
                sum_sq = 0.0                                    ! Summation using the acumulation.
                do i = 1, part_num
                    dx = x(i, actual_step) - x(i, 1)
                    dy = y(i, actual_step) - y(i, 1)                      ! Difference of positions is between the xj and x1 (as reference).
                    dz = z(i, actual_step) - z(i, 1)
                    sum_sq = sum_sq + (dx**2 + dy**2 + dz**2)
                end do

                partial_rmsd(j) = sqrt(sum_sq / part_num)               ! accumulate partials rmsd 
            end do                                                    ! rmsd=sqrt(summation(r-r')^2 / n)

            ! Merge the partial rmsd to totalrmsd
            call mpi_gather (partial_rmsd, partial_steps, mpi_real, total_rmsd, partial_steps, mpi_real, 0, mpi_comm_world, ierr )

            ! Rank 0 does the final work of output
            if (rank == 0) then
                open(12, file = rmsd_file, status = 'replace')
                do j = 1, step_num
                    write(12, *) time(j), total_rmsd(j)
                end do
                close(12)

                print *, 'RMSD calculation completed and saved to ', rmsd_file
            endif

            deallocate(partial_rmsd)
            if (rank == 0) then
                deallocate(total_rmsd)
            endif

    end subroutine compute_rmsd
end module post_trajectory_analysis
