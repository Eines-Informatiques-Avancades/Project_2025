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
    use global_vars
    use lj_forces
    use lj_potentials
    use geometry
    use initial_conf
    use integrators
    use io
    use post_trajectory_analysis
    use thermodynamics
    use thermostat

    implicit none

    integer :: step, seed_size, i
    integer, allocatable :: seed(:)
    real :: part_density, volume, time, lj_potential, temperature_inst, &
        kinetic_energy, total_energy, kinetic_energy_x, kinetic_energy_y, &
        kinetic_energy_z
    real, allocatable :: x(:, :), y(:, :), z(:, :), time_points(:), &
        initial_positions(:, :), positions_x(:), positions_y(:), positions_z(:), &
        forces_x(:), forces_y(:), forces_z(:), velocities_x(:), velocities_y(:), velocities_z(:)
    character(50) :: input_file, positions_file, thermodynamics_file, rdf_file, rmsd_file

    ! System parameters.
    input_file = 'input_parameters.in'
    call read_input(input_file)

    volume = system_size**(3.)  ! System is a cubic box.

    !
    ! Generate initial system configuration.
    !

    ! Generate the initial configuration from a lattice.
    call gen_initial_conf(part_density, initial_positions)
    print *, 'Initial lattice particle density: ', part_density
    print *

    allocate(positions_x(part_num), positions_y(part_num), positions_z(part_num))

    do i = 1, part_num
        positions_x(i) = initial_positions(i, 1)
        positions_y(i) = initial_positions(i, 2)
        positions_z(i) = initial_positions(i, 3)
    end do

    ! Center initial config at the origin of coordinates.
    call apply_pbc(positions_x)
    call apply_pbc(positions_y)
    call apply_pbc(positions_z)

    allocate(velocities_x(part_num), velocities_y(part_num), velocities_z(part_num))
    call gen_velocities_bimodal_distr(velocities_x)
    call gen_velocities_bimodal_distr(velocities_y)
    call gen_velocities_bimodal_distr(velocities_z)

    print *, 'Computing initial Lennard-Jones forces...'
    call compute_forces(positions_x, positions_y, positions_z, forces_x, forces_y, forces_z)

    print *, 'Generating initial configuration for a VdW gas from the lattice...'

    ! Seed initialization for the Andersen_thermsostat.
    call random_seed(size = seed_size)
    allocate(seed(seed_size))

    if (test_mode == "ON") then
        seed = 123456789    ! Set arbitrary seed to all elements.
        call random_seed(put = seed)
    endif

    do step = 1, equilibration_step_num
        call compute_forces(positions_x, positions_y, positions_z, forces_x, forces_y, forces_z)

        call velocity_verlet(timestep, positions_x, velocities_x, forces_x)
        call andersen_thermostat(velocities_x)

        call velocity_verlet(timestep, positions_y, velocities_y, forces_y)
        call andersen_thermostat(velocities_y)

        call velocity_verlet(timestep, positions_z, velocities_z, forces_z)
        call andersen_thermostat(velocities_z)
    end do

    deallocate(seed)

    !
    ! System evolution.
    !

    ! Create a new positions_file or replace the existing one.
    ! This needs to be done previously, as the write_positions_xyz subroutine
    ! appends the positions to an existing file.
    positions_file = 'positions.xyz'
    open(4, file = positions_file, status = 'replace')
    close(4)

    thermodynamics_file = 'thermodynamics.dat'
    open(12, file = thermodynamics_file, status = 'replace')
    write(12, *) '# time, lj_potential, kinetic_energy, total_energy, temperature_inst'

    time = 0
    do step = 1, step_num
        time = time + timestep

        call compute_forces(positions_x, positions_y, positions_z, forces_x, forces_y, forces_z)

        call velocity_verlet(timestep, positions_x, velocities_x, forces_x)
        call andersen_thermostat(velocities_x)

        call velocity_verlet(timestep, positions_y, velocities_y, forces_y)
        call andersen_thermostat(velocities_y)

        call velocity_verlet(timestep, positions_z, velocities_z, forces_z)
        call andersen_thermostat(velocities_z)

        call compute_lj(positions_x, positions_y, positions_z, lj_potential)
        call compute_total_kinetic_energy(velocities_x, kinetic_energy_x)
        call compute_total_kinetic_energy(velocities_y, kinetic_energy_y)
        call compute_total_kinetic_energy(velocities_z, kinetic_energy_z)

        kinetic_energy = kinetic_energy_x + kinetic_energy_y + kinetic_energy_z
        total_energy = lj_potential + kinetic_energy
        temperature_inst = instantaneous_temperature(kinetic_energy)

        write(12, *) time, lj_potential, kinetic_energy, total_energy, temperature_inst
        call write_positions_xyz(time, positions_x, positions_y, positions_z, positions_file)
    end do

    close(12)

    deallocate( &
        positions_x, positions_y, positions_z, &
        velocities_x, velocities_y, velocities_z, &
        forces_x, forces_y, forces_z &
    )

    !
    ! Post-trajectory analysis (RDF and RMSD computation).
    !

    print *
    print *, 'Performing post-trajectory analysis...'

    allocate( &
        x(part_num, step_num), y(part_num, step_num), &
        z(part_num, step_num), time_points(step_num) &
    )

    call read_trajectory(positions_file, x, y, z, time_points)

    rdf_file = 'rdf.dat'
    rmsd_file = 'rmsd.dat'
    call compute_rdf(x, y, z, rdf_file)
    call compute_rmsd(x, y, z, time_points, rmsd_file)

    deallocate(x, y, z, time_points)
end program vdw_gas
