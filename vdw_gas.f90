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

    integer :: part_num, step_num, step, seed_size
    real :: part_density, system_size, volume, cutoff, time, timestep, lj_potential, &
        temperature, temperature_inst, collision_frequence, kinetic_energy, total_energy
    real, allocatable :: positions(:, :), forces(:, :), velocities(:, :)
    real, allocatable :: x(:, :), y(:, :), z(:, :), time_points(:)
    character(6) :: lattice_type
<<<<<<< HEAD
    character(50) :: positions_file, input_file, rdf_file, rmsd_file
=======
    character(3) :: test_mode
    integer, allocatable :: seed(:)
    character(50) :: positions_file, input_file
>>>>>>> 01a9310 (seed added to control thermostat random numbers)

    ! System parameters.
    input_file = 'input_parameters.in'
    call read_input(input_file, part_num, system_size, lattice_type, timestep, step_num, temperature, &
        collision_frequence, cutoff, test_mode)

    volume = system_size**(3.)  ! System is a cubic box.

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

    
    ! seed inizialization for andersen_thermsostat
    call random_seed(size=seed_size)
    allocate(seed(seed_size))

    if (test_mode == "ON") then
        seed = 123456789    ! putting arbitrary seed to all elements
    endif

    call random_seed(put=seed)

    do step = 1, step_num
        call velocity_verlet(timestep, part_num, system_size, cutoff, positions, velocities, lj_potential)
        call andersen_thermostat(part_num, temperature, collision_frequence, velocities)
    end do

<<<<<<< HEAD
    !
    ! System evolution.
    !
=======
    deallocate(seed)
>>>>>>> 01a9310 (seed added to control thermostat random numbers)

    ! Create a new positions_file or replace the existing one.
    positions_file = 'positions.xyz'
    open(4, file = positions_file, status = 'replace')
    write(4, *) '# time, particle coordinates'
    close(4)

    open(10, file = 'temperature_inst.dat', status = 'replace')
    open(11, file = 'lj_potential.dat', status = 'replace')
    open(12, file = 'kinetic_energy.dat', status = 'replace')
    open(13, file = 'total_energy.dat', status = 'replace')
    time = 0
    do step = 1, step_num
        time = time + timestep

        call velocity_verlet(timestep, part_num, system_size, cutoff, positions, velocities, lj_potential)
        call andersen_thermostat(part_num, temperature, collision_frequence, velocities)

        call compute_total_kinetic_energy(part_num, velocities, kinetic_energy)
        total_energy = lj_potential + kinetic_energy

        temperature_inst = instantaneous_temperature(part_num, kinetic_energy)

        write(10, *) time, temperature_inst
        write(11, *) time, lj_potential
        write(12, *) time, kinetic_energy
        write(13, *) time, total_energy
        call write_positions_xyz(part_num, time, positions, positions_file)
    end do
    close(10)
    close(11)
    close(12)
    close(13)

    deallocate(positions, velocities, forces)

    !
    ! Post-trajectory analysis (RDF and RMSD computation).
    !

    print *, 'Performing post-trajectory analysis...'

    allocate( &
        x(part_num, step_num), y(part_num, step_num), &
        z(part_num, step_num), time_points(step_num) &
    )

    call read_trajectory(part_num, step_num, positions_file, x, y, z, time_points)

    rdf_file = 'rdf.dat'
    rmsd_file = 'rmsd.dat'
    call compute_rdf(part_num, step_num, system_size, x, y, z, rdf_file)
    call compute_rmsd(part_num, step_num, x, y, z, time_points, rmsd_file)

    deallocate(x, y, z, time_points)
end program vdw_gas
