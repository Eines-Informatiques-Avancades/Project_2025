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

    integer :: part_num, time
    real :: part_density, system_size, volume, cutoff, dt, temperature, collision_frequence
    real, allocatable :: positions(:, :), forces(:, :), velocities(:,:)
    character(6) :: lattice_type
    character(50) :: positions_file
	
    allocate(positions(system_size,3), velocities(system_size,3), forces(system_size,3))
    velocities = 0
    ! Generate initial lattice.
    lattice_type = 'SC'
    part_num = 125
    system_size = 500
    dt = 0.01
    temperature = 10.0
    collision_frequence = 0.05

    volume = system_size**(3.)  ! System is a cubic box.
    cutoff = 0.5*system_size    ! Cutoff radius for molecular interactions.

    call gen_initial_conf(lattice_type, system_size, part_num, part_density, positions)

    ! Center initial config at the origin of coordinates.
    call apply_pbc(positions, system_size)

    positions_file = 'positions.xyz'

    ! Create a new positions_file or replace the existing one.
    open(4, file = positions_file, status = 'replace')
    write(4, *) '# time, particle coordinates'
    close(4)

    do time = 1, 5
        call write_positions_xyz(part_num, time, positions, positions_file)
	
	call velocity_verlet(dt, part_num, system_size, cutoff, positions, velocities)
	call andersen_thermostat(part_num, temperature, collision_frequence, velocities)
    end do
end program vdw_gas
