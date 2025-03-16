!
! post_trajectory_analysis
! Molecular Dynamics Simulation of a Van der Waals Gas
! Haoyu Huang, Ricard Rodriguez
!

! Read data from the positions file and store it into arrays.
! Particle positions in time are stored in columns: time x y z
! Each row refers to a particle in an specific time.
subroutine read_trajectory(part_num, step_num, positions_file, x, y, z, time)
    implicit none

    integer, intent(in) :: part_num, step_num
    character(50), intent(in) :: positions_file
    real, allocatable, intent(inout) :: x(:, :), y(:, :), z(:, :), time(:)

    integer :: part, step, ios
    real :: t

    print *, 'Reading trajectory data from ', positions_file

    open(4, file = positions_file, status = 'old', action = 'read')

    do step = 1, step_num
        do part = 1, part_num
            read(4, *, iostat = ios) t, x(part, step), y(part, step), z(part, step)

            ! Store the different values of t (once for each different time).
            if (part == 1) then
                time(step) = t
            end if
        end do
    end do

    close(4)
end subroutine read_trajectory

! Test to proove that the program has access to all the stored xyz information.
subroutine test_access_xyz_data(part_num, step_num, x, y, z, time)
    implicit none

    integer, intent(in) :: part_num, step_num
    real, allocatable, intent(in) :: x(:, :), y(:, :), z(:, :), time(:)

    integer :: i, j

    print *, 'Processing stored coordinates...'

    do j = 1, step_num
        print *, 'time:', time(j)

        do i = 1, part_num
            print *,'(', x(i, j), y(i, j), z(i, j), ')'
        end do
    end do
end subroutine test_access_xyz_data

! Test to access specific step xyz information.
subroutine test_access_xyz_frame(part_num, step_num, x, y, z, time)
    implicit none

    integer, intent(in) :: part_num, step_num
    real, allocatable, intent(in) :: x(:, :), y(:, :), z(:, :), time(:)

    integer :: i, j

    print *, 'Testing specific step xyz information...'

    do j = 1, step_num
        if (j == 3) then                            ! test for step = 3
            print *, 'time:', time(j)

            do i = 1, part_num
                print *, 'x', x(i, j), 'y', y(i, j), 'z', z(i, j)
            end do
        end if
    end do
end subroutine test_access_xyz_frame

! Compute RDF using the stored data.
! Must be executed after read_trajectory, as it depends on the arrays it
! creates.
subroutine compute_rdf(part_num, step_num, system_size, x, y, z, rdf_file)
    implicit none

    integer, intent(in) :: part_num, step_num
    real, intent(in) :: system_size
    real, allocatable, intent(in) :: x(:, :), y(:, :), z(:, :)
    character(50), intent(in) :: rdf_file

    integer :: i, j, k, time_index
    real(4) :: maximum_radius                      ! maximum radius
    real(4), parameter :: dr = 0.05
    integer :: bins                                ! Define a number of bins to set the size of rdf and r vectors
    real(4), allocatable :: rdf(:), r_values(:)    ! dx,dy,dz,dr and dv are respectively x,y,z,r positions variation and volume variation
    real(4) :: r, dx, dy, dz, dv, density
    integer :: bin_index                           ! bin_index correspond to the zone of sphere that this r belongs
    real(4) :: volume                              ! volume refers to the sphere's volume

    maximum_radius = 0.5 * system_size
    bins = int(maximum_radius / dr)
    volume = (4.0/3.0) * 3.1415926 * maximum_radius**3
    density = part_num / volume

    allocate(rdf(bins), r_values(bins))
    rdf = 0.0

    do time_index = 1, step_num                    ! loop for all registered time
        do i = 1, part_num - 1                     ! loop for all atoms i.e. i=1 j=2 dx=x1-x2 ...
            do j = i + 1, part_num                 ! loop for the next atom (atom j) of atom i
                dx = x(j, time_index) - x(i, time_index)
                dy = y(j, time_index) - y(i, time_index)
                dz = z(j, time_index) - z(i, time_index)
                r = sqrt(dx**2 + dy**2 + dz**2)

                if (r < maximum_radius) then
                    bin_index = int(r / dr) + 1
                    rdf(bin_index) = rdf(bin_index) + 1
                end if
            end do
        end do
    end do

    do k = 1, bins
        r_values(k) = k * dr                             ! r_values are grid points of r
        dv = 4.0 * 3.14159 * r_values(k)**2 * dr         ! volume between r to r + dr to normalize the rdf
        rdf(k) = rdf(k) / (density * part_num * dv )     ! normalize the rdf
    end do

    open(12, file = rdf_file, status = 'replace')
    do k = 1, bins
        write(12, *) r_values(k), rdf(k)
    end do
    close(12)

    print *, 'RDF calculation completed and saved to ', rdf_file

    deallocate(rdf, r_values)
end subroutine compute_rdf

! Compute RMSD using the stored data.
! Must be executed after read_trajectory, as it depends on the arrays it
! creates.
subroutine compute_rmsd(part_num, step_num, x, y, z, time, rmsd_file)
    implicit none

    integer, intent(in) :: part_num, step_num
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
