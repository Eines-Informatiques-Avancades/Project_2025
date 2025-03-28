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
            implicit none

            character(50), intent(in) :: positions_file
            real, allocatable, intent(inout) :: x(:, :), y(:, :), z(:, :), time(:)

            integer :: part, step, ios, n
            real :: t

            print *, 'Reading trajectory data from ', positions_file

            open(4, file = positions_file, status = 'old', action = 'read')

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
        end subroutine read_trajectory

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
            real(8), parameter :: dr = 0.01
            real(8) :: maximum_radius, volume, density
            integer :: bins
            real(8), allocatable :: h(:), rdf(:), r_values(:)
            real(8) :: r, r_sq, dx, dy, dz, dv, h_id, r_lo, r_hi, const, nid
            integer :: bin_index

            ! Parameters
            maximum_radius = system_size / 2
            bins = int(maximum_radius / dr)
            volume = system_size**3
            density = part_num / volume

            ! Allocate arrays
            allocate(h(bins), rdf(bins), r_values(bins))
            h = 0.0

            ! Compute histogram h(k)
            do time_index = 1, step_num
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
                            h(bin_index) = h(bin_index) + 2   ! pairwise counting
                        endif
                             
                    end do
                end do
            end do

            ! Normalize RDF
            const = 4.0 * 3.14159265358979 / 3.0
            do k = 1, bins
                r_lo = (k - 1) * dr
                r_hi = r_lo + dr
                dv = const * (r_hi**3 - r_lo**3)  ! Shell volume
                nid = density * dv
                rdf(k) = h(k) / (part_num * step_num * nid)
                r_values(k) = (k - 0.5) * dr  ! Bin center
            end do

            ! Save RDF results to file
            open(12, file = rdf_file, status = 'replace')
             do k = 1, bins
                write(12, '(F10.5, F15.8)') r_values(k), rdf(k)
             end do
            close(12)

            print *, 'RDF calculation completed and saved to ', rdf_file

            ! Deallocate memory
            deallocate(h, rdf, r_values)

        end subroutine compute_rdf

        ! Compute RMSD using the stored data.
        ! Must be executed after read_trajectory, as it depends on the arrays it
        ! creates.
        subroutine compute_rmsd(x, y, z, time, rmsd_file)
            implicit none

            real, allocatable, intent(in) :: x(:, :), y(:, :), z(:, :), time(:)
            character(50), intent(in) :: rmsd_file

            integer :: i, j
            real(4), allocatable :: rmsd(:)
            real(4) :: dx, dy, dz, sum_sq

            allocate(rmsd(step_num))

            do j = 1, step_num                                  ! for all steps
                sum_sq = 0.0                                    ! summation using the acumulation
                do i = 1, part_num
                    dx = x(i, j) - x(i, 1)
                    dy = y(i, j) - y(i, 1)                      ! difference of positions is between the xj and x1 (as reference).
                    dz = z(i, j) - z(i, 1)
                    sum_sq = sum_sq + (dx**2 + dy**2 + dz**2)
                end do

                rmsd(j) = sqrt(sum_sq / part_num)               ! rmsd=sqrt(summation(r-r')^2 / n)
            end do

            open(12, file = rmsd_file, status = 'replace')
             do j = 1, step_num
                write(12, *) time(j), rmsd(j)
             end do
            close(12)

            print *, 'RMSD calculation completed and saved to ', rmsd_file

            deallocate(rmsd)
        end subroutine compute_rmsd
end module post_trajectory_analysis
