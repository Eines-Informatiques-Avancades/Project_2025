!
! vdw_gas
! Molecular Dynamics Simulation of a Van der Waals Gas
! Ricard Rodriguez
!
! Main program corresponding to the "Molecular Dynamics Simulation of a Van der
! Waals Gas" project for the Advanced Computer Tools subject. Please check the
! rest of the files containing the module(s) necessary for the program
! execution.
!

program vdw_gas
    use subroutines

    implicit none

    integer :: part_num, step_num, step
    real :: part_density, system_size, volume, cutoff, time, timestep, lj_potential
    real, allocatable :: positions(:, :), forces(:, :), velocities(:, :)
    character(6) :: lattice_type
    character(50) :: positions_file

    ! System parameters.
    lattice_type = 'SC'
    part_num = 125
    system_size = 6
    timestep = 0.0001
    step_num = 10000

    volume = system_size**(3.)  ! System is a cubic box.
    cutoff = 1.5                ! Cutoff radius for molecular interactions.

    !
    ! Generate initial system configuration.
    !

    ! Generate the initial configuration from a lattice.
    call gen_initial_conf(lattice_type, system_size, part_num, part_density, positions)
    print *, 'Initial lattice particle density: ', part_density

    ! Center initial config at the origin of coordinates.
    call apply_pbc(positions, system_size)

    allocate(velocities(part_num, 3))
    velocities(:, :) = 0

    print *, 'Computing initial Lennard-Jones forces...'
    call compute_forces(part_num, positions, forces, lj_potential, system_size, cutoff)

    ! Create a new positions_file or replace the existing one.
    positions_file = 'positions.xyz'
    open(4, file = positions_file, status = 'replace')
    write(4, *) '# time, particle coordinates'
    close(4)

    open(5, file = 'lj_potential.dat')
    time = 0
    do step = 1, step_num
        time = time + timestep

        call velocity_verlet(timestep, part_num, system_size, cutoff, positions, velocities)

        call compute_forces(part_num, positions, forces, lj_potential, system_size, cutoff)

        write(5, *) time, lj_potential
        call write_positions_xyz(part_num, time, positions, positions_file)
    end do
    close(5)
end program vdw_gas
