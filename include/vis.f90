!
! post_trajectory_analysis
! Molecular Dynamics Simulation of a Van der Waals Gas
! Haoyu Huang, Ricard Rodriguez
!


! This module is used for the storage of the variables
! Also it ensures that these variables can be accessed by the subroutine.
module var
    implicit none

    integer :: part_num, step_num
    real :: timestep, part_density, temperature, system_size, collision_frequence
    character(6) :: lattice_type
    character(50) :: input_file

    real, allocatable :: x(:,:), y(:,:), z(:,:), time(:)
end module var

program post_trajectory_analysis
    use var
    use subroutines, only: read_input

    implicit none

    ! read initial configuration to obtain the parameters of the simulation
    input_file = 'input_parameters.in'
    call read_input(input_file, part_num, system_size, lattice_type, timestep, step_num, temperature, collision_frequence)

    allocate( &
        x(part_num, step_num), y(part_num, step_num), &
        z(part_num, step_num), time(step_num) &
    )

    call read_trajectory()

    ! call test_access_xyz_data()
    ! call test_access_xyz_frame()

    call compute_rdf()
    call compute_rmsd()

    ! Release the stored memory
    deallocate(x, y, z, time)
end program post_trajectory_analysis


subroutine read_trajectory()
    use var

    implicit none

    integer :: i, j, ios
    real :: t

    open(1, file = 'positions.xyz', status = 'old', action = 'read')

    ! Read information line by line.
    ! i.e. at time j=1 or time=1, read n atoms xyz information, then for
    ! j = 2, 3, ... the same
    do j = 1, step_num
        do i = 1, part_num
            ! this 'xyz' format is： time x y z
            ! these positions are vectorized following the the atom number
            ! (i = atom1, atom2, ...) and the time that is registered
            ! (j = 1, 2, ...)
            read(1, *, iostat = ios) t, x(i, j), y(i, j), z(i, j)

            ! store the time information (it's enough with register first atom
            ! time)
            if (i == 1) then
                time(j) = t
            end if
        end do
    end do

    close(1)

    print *, 'Finished reading the trajectories.'
end subroutine read_trajectory

! Test to proove that the program has access to all the stored xyz information.
subroutine test_access_xyz_data()
    use var

    implicit none

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
subroutine test_access_xyz_frame()
    use var

    implicit none

    integer :: i, j

    print *, 'Testing specific step xyz information:'

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
subroutine compute_rdf()
    use var

    implicit none

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
            do j = i+1, part_num                   ! loop for the next atom (atom j) of atom i
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

    open(2, file = 'rdf_data.txt', status = 'replace')
    do k = 1, bins
        write(2, *) r_values(k), rdf(k)
    end do
    close(2)

    print *, 'RDF calculation completed and saved to rdf_data.txt'

    deallocate(rdf, r_values)
end subroutine compute_rdf

! Compute RMSD using the stored data.
subroutine compute_rmsd()
    use var

    implicit none

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

    open(2, file = 'rmsd_data.txt', status = 'replace')
    do j = 1, step_num
        write(2, *) time(j), rmsd(j)
    end do
    close(2)

    print *, 'RMSD calculation completed and saved to rmsd_data.txt'

    deallocate(rmsd)
end subroutine compute_rmsd
