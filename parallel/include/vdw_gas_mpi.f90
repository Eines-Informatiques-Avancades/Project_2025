!
! vdw_gas_mpi
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
    use geometry
    use initial_conf
    use integrators
    use io
    use post_trajectory_analysis
    use thermodynamics
    use thermostat
    use mpi

    implicit none
    
    ! Particles per worker, and worker starting index
    integer, allocatable :: counts(:), displs(:)
    integer :: ierror, rank, nproc

    integer :: step, seed_size, i, MASTER
    integer, allocatable :: seed(:)
    real(8) :: volume, time, lj_potential, temperature_inst, &
        kinetic_energy, total_energy
    real(8) :: part_density
    real(8), allocatable :: positions(:, :), forces(:, :), velocities(:, :), &
        x(:, :), y(:, :), z(:, :), time_points(:)
    character(50) :: input_file, positions_file, thermodynamics_file, rdf_file, rmsd_file

    call MPI_INIT(ierror)

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierror) 

    MASTER = 0

    ! System parameters.
    input_file = 'input_parameters.in'
    call read_input(input_file)

        ! Distribute the particles between workers for integrators
        ! N-remainder will work with Nlocal particles and remainder # of particles with Nlocal + 1

    allocate(counts(0:nproc-1))
    allocate(displs(0:nproc-1))


    do i = 0, nproc-1
        if (i < mod(part_num, nproc)) then
            counts(i) = (part_num / nproc + 1) 
        else
            counts(i) = (part_num / nproc) 
        end if
    end do
    print*, rank, counts(rank)
    displs(0) = 0

    do i = 1, nproc-1
        displs(i) = displs(i-1) + counts(i-1) 
    end do

    volume = system_size**(3.)  ! System is a cubic box.




    call MPI_Barrier(MPI_COMM_WORLD, ierror)
    print*, rank, "count", counts(rank), "displacement", displs(rank)

    !
    ! Generate initial system configuration.
    !
    if (rank == MASTER) then
        ! Generate the initial configuration from a lattice.
        call gen_initial_conf(part_density, positions)
        print *, 'Initial lattice particle density: ', part_density
        print *

        ! Center initial config at the origin of coordinates.
        do i = 1, 3
            call apply_pbc(positions(:,3))
        end do
    end if
    allocate(velocities(part_num, 3))
    call gen_velocities_bimodal_distr(velocities)

    print *, 'Computing initial Lennard-Jones forces...'
    call compute_forces(positions, forces, lj_potential)
    print *, 'Generating initial configuration for a VdW gas from the lattice...'

    ! Seed initialization for the Andersen_thermsostat.
    call random_seed(size = seed_size)
    allocate(seed(seed_size))

    if (test_mode == "ON") then
        seed = 123456789    ! Set arbitrary seed to all elements.
        call random_seed(put = seed)
    endif

    do step = 1, equilibration_step_num
        call velocity_verlet(timestep, positions, velocities, lj_potential, counts, displs, rank, nproc)
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
        call velocity_verlet(timestep, positions, velocities, lj_potential, counts, displs, rank, nproc)
        call andersen_thermostat(velocities)

        call compute_total_kinetic_energy(velocities, kinetic_energy)

        total_energy = lj_potential + kinetic_energy
        temperature_inst = instantaneous_temperature(kinetic_energy)

        write(12, *) time, lj_potential, kinetic_energy, total_energy, temperature_inst
        call write_positions_xyz(time, positions, positions_file)
    end do

    close(12)

    deallocate(positions, velocities, forces)

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

    call MPI_FINALIZE(ierror);
end program vdw_gas
