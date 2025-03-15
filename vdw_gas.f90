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
    real :: part_density, system_size, volume, cutoff, time, timestep, lj_potential, &
        temperature, temperature_inst, collision_frequence, kinetic_energy, total_energy
    real, allocatable :: positions(:, :), forces(:, :), velocities(:, :)
    character(6) :: lattice_type
    character(50) :: positions_file, input_file

    ! System parameters.
    input_file = 'input_parameters.in'
    call read_input(input_file, part_num, system_size, lattice_type, timestep, step_num, temperature, collision_frequence)

    volume = system_size**(3.)  ! System is a cubic box.
    cutoff = 1.5                ! Cutoff radius for molecular interactions.

    !
    ! Generate initial system configuration.
    !

    ! Generate the initial configuration from a lattice.
    call gen_initial_conf(lattice_type, system_size, part_num, part_density, positions)
    print *, 'Initial lattice particle density: ', part_density
    print *

    ! Center initial config at the origin of coordinates.
    call apply_pbc(positions, system_size)

    allocate(velocities(part_num, 3))
    velocities(:, :) = 0

    print *, 'Computing initial Lennard-Jones forces...'
    call compute_forces(part_num, positions, forces, lj_potential, system_size, cutoff)

    print *, 'Generating initial configuration for a VdW gas from the lattice...'

    do step = 1, step_num
        call velocity_verlet(timestep, part_num, system_size, cutoff, positions, velocities, lj_potential)
        call andersen_thermostat(part_num, temperature, collision_frequence, velocities)
    end do

    !
    ! System evolution.
    !

    ! Create a new positions_file or replace the existing one.
    positions_file = 'positions.xyz'
    open(4, file = positions_file, status = 'replace')
    write(4, *) '# time, particle coordinates'
    close(4)

    open(5, file = 'temperature_inst.dat', status = 'replace')
    open(6, file = 'lj_potential.dat', status = 'replace')
    open(7, file = 'kinetic_energy.dat', status = 'replace')
    open(8, file = 'total_energy.dat', status = 'replace')
    time = 0
    do step = 1, step_num
        time = time + timestep

        call velocity_verlet(timestep, part_num, system_size, cutoff, positions, velocities, lj_potential)
        call andersen_thermostat(part_num, temperature, collision_frequence, velocities)

        call compute_total_kinetic_energy(part_num, velocities, kinetic_energy)
        total_energy = lj_potential + kinetic_energy

        temperature_inst = instantaneous_temperature(part_num, kinetic_energy)

        write(5, *) time, temperature_inst
        write(6, *) time, lj_potential
        write(7, *) time, kinetic_energy
        write(8, *) time, total_energy
        call write_positions_xyz(part_num, time, positions, positions_file)
    end do
    close(5)
    close(6)
    close(7)
    close(8)
end program vdw_gas
