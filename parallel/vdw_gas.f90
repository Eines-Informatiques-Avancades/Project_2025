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
    use geometry
    use initial_conf
    use integrators
    use io
    use post_trajectory_analysis
    use thermodynamics
    use thermostat

    implicit none

    integer :: i, step, seed_size, tclockstart, tclockend, clock_rate
    integer, allocatable :: seed(:), verlet_list(:, :), n_neighbors(:)
    real(8) :: part_density, volume, time, lj_potential, temperature_inst, &
        kinetic_energy, total_energy
    real(8), allocatable :: positions(:, :), forces(:, :), velocities(:, :), &
        x(:, :), y(:, :), z(:, :), time_points(:)
    character(50) :: input_file, positions_file, thermodynamics_file, rdf_file, rmsd_file

    ! CPU Time measurement.
    real ::  tcpustart, tcpuend


    ! Initialize MPI
    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)

    if (rank == 0) then
        ! Initialize simulation time measurement.
        call system_clock(count_rate = clock_rate)  ! Set reference clock rate (tick/s).
        call system_clock(tclockstart)              ! Start actual tick couting.
        call cpu_time(tcpustart)                    ! Start counting cpu time.

        ! Read simulation parameters.
        input_file = 'input_parameters.in'
        call read_input(input_file)

        volume = system_size**(3.)  ! System is a cubic box.

        ! Initialize a flag to let the user know if at a certain time of the
        ! simulation a particle has escaped the box with infinite velocity.
        infinite_distance = .false.
    end if

    ! Broadcast simulation parameters.
    call mpi_bcast(part_num                 , 1                 , MPI_INTEGER   , 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(atom_type                , len(atom_type)    , MPI_CHARACTER , 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(system_size              , 1                 , MPI_REAL8     , 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(lattice_type             , len(lattice_type) , MPI_CHARACTER , 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(timestep                 , 1                 , MPI_REAL8     , 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(step_num                 , 1                 , MPI_INTEGER   , 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(equilibration_step_num   , 1                 , MPI_INTEGER   , 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(temperature              , 1                 , MPI_REAL8     , 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(collision_frequence      , 1                 , MPI_REAL8     , 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(cutoff                   , 1                 , MPI_REAL8     , 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(test_mode                , len(test_mode)    , MPI_CHARACTER , 0, MPI_COMM_WORLD, ierr)

    !
    ! MPI variable definition.
    !

    allocate( &
        counts(0:nproc - 1), &
        displs(0:nproc - 1)  &
    )

    do i = 0, nproc - 1
        if (i < mod(part_num, nproc)) then
            counts(i) = (part_num / nproc + 1)
        else
            counts(i) = (part_num / nproc)
        end if
    end do

    displs(0) = 0
    do i = 1, nproc - 1
        displs(i) = displs(i - 1) + counts(i - 1)
    end do

    start_part = displs(rank) + 1
    end_part = displs(rank) + counts(rank)

    !
    ! Generate initial system configuration.
    !

    allocate( &
        velocities(part_num, 3), positions(part_num, 3), forces(part_num, 3), &
        n_neighbors(part_num), verlet_list(part_num, part_num) &
    )

    if (rank == 0) then
        ! Generate the initial configuration from a lattice.
        call gen_initial_conf(part_density, positions)
        print *, 'Initial lattice particle density: ', part_density
        print *
    end if
    call mpi_bcast(positions, size(positions), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    ! Center initial config at the origin of coordinates.
    call apply_pbc(positions, counts, displs)

    if (rank == 0) call gen_velocities_bimodal_distr(velocities)
    call mpi_bcast(velocities, size(velocities), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
        print *, 'Computing initial Verlet lists...'
    end if
    call compute_verlet_list(positions, verlet_list, n_neighbors)

    if (rank == 0) print *, 'Computing initial Lennard-Jones forces...'
    call compute_forces(positions, forces, lj_potential, verlet_list, n_neighbors)

    if (rank == 0) print *, 'Generating initial configuration for a VdW gas from the lattice...'

    ! Seed initialization for the Andersen thermostat.
    call random_seed(size = seed_size)
    allocate(seed(seed_size))

    if (test_mode == "ON") then
        seed = 123456789    ! Set arbitrary seed to all elements.
        call random_seed(put = seed)
    endif

    ! Force all threads to start thermalization at the same time.
    call mpi_barrier(MPI_COMM_WORLD, ierr)

    do step = 1, equilibration_step_num
        if (mod(step, 10) == 0) then
            call compute_verlet_list(positions, verlet_list, n_neighbors)
        end if

        call velocity_verlet(positions, velocities, lj_potential, verlet_list, n_neighbors, counts, displs)
        call andersen_thermostat(velocities, counts, displs)
    end do

    deallocate(seed)

    if (rank == 0) then
        call cpu_time(tcpuend)                  ! Stop counting cpu time.
        call system_clock(tclockend)            ! Stop tick counting.

        print *
        print *, 'Generation of the initial system configuration complete.'
        print *, 'Cputime: ', tcpuend - tcpustart, 's'
        print *, 'Wallclock time: ', real(tclockend - tclockstart) / clock_rate, 's'
    end if

    !
    ! System evolution.
    !

    if (rank == 0) then
        print *
        print *, 'Computing evolution of the particle system...'

        call cpu_time(tcpustart)
        call system_clock(tclockstart)

        ! Create a new positions_file or replace the existing one.
        ! This needs to be done previously, as the write_positions_xyz subroutine
        ! appends the positions to an existing file.
        positions_file = 'positions.xyz'
        open(4, file = positions_file, status = 'replace')

        thermodynamics_file = 'thermodynamics.dat'
        open(12, file = thermodynamics_file, status = 'replace')
        write(12, *) '# time, lj_potential, kinetic_energy, total_energy, temperature_inst'
    end if

    time = 0
    do step = 1, step_num
        time = time + timestep

        if (mod(step, 10) == 0) then
            call compute_verlet_list(positions, verlet_list, n_neighbors)
        end if

        call velocity_verlet(positions, velocities, lj_potential, verlet_list, n_neighbors, counts, displs)
        call andersen_thermostat(velocities, counts, displs)

        call compute_total_kinetic_energy(velocities, kinetic_energy)

        total_energy = lj_potential + kinetic_energy
        temperature_inst = instantaneous_temperature(kinetic_energy)

        if (rank == 0) then
            write(12, '(f12.6, 2x, f12.6, 2x, f12.6, 2x, f12.6, 2x, f12.6)') &
                time, lj_potential, kinetic_energy, total_energy, temperature_inst

            call write_positions_xyz(time, positions, positions_file)
        end if
    end do

    if (rank == 0) then
        close(4)
        close(12)
    end if

    deallocate(positions, velocities, forces)

    if (infinite_distance) then
        print *, 'Particles escaping really far from the system box have been detected during the simulation.'
        print *, 'Your system might present unstabilities.'
        print *, 'Check your simulation results and consider using a smaller timestep.'
    end if

    if (rank == 0) then
        call cpu_time(tcpuend)
        call system_clock(tclockend)

        print *
        print *, 'Study of the system evolution (production loops) complete.'
        print *, 'Particle trajectories saved to ', positions_file
        print *, 'Measures of thermodynamical variables saved to ', thermodynamics_file
        print *, 'Cputime: ', tcpuend - tcpustart, 's'
        print *, 'Wallclock time: ', real(tclockend - tclockstart) / clock_rate, 's'
    end if

    !
    ! Post-trajectory analysis (RDF and RMSD computation).
    !

    if (rank == 0) then
        call cpu_time(tcpustart)
        call system_clock(tclockstart)

        print *
        print *, 'Performing post-trajectory analysis...'
    end if

    allocate( &
        x(part_num, step_num), y(part_num, step_num), &
        z(part_num, step_num), time_points(step_num) &
    )

    call read_trajectory(positions_file, x, y, z, time_points)

    if (rank == 0) then
	call cpu_time(tcpuend)
        call system_clock(tclockend)

	print *, 'Cputime: ', tcpuend-tcpustart, 's'
        print *, 'Wallclock time: ', real(tclockend - tclockstart) / clock_rate, 's'
	print *
    endif

    rdf_file = 'rdf.dat'
    rmsd_file = 'rmsd.dat'

    if (rank == 0) then
        call cpu_time(tcpustart)
        call system_clock(tclockstart)
    endif

    call compute_rdf(x, y, z, rdf_file)
    
    if (rank == 0) then
	call cpu_time(tcpuend)
        call system_clock(tclockend)

	print *, 'Cputime: ', tcpuend-tcpustart, 's'
        print *, 'Wallclock time: ', real(tclockend - tclockstart) / clock_rate, 's'
	print *
    endif

    
    if (rank == 0) then
        call cpu_time(tcpustart)
        call system_clock(tclockstart)
    endif
	
    call compute_rmsd(x, y, z, time_points, rmsd_file)

    deallocate(x, y, z, time_points)

    if (rank == 0) then
        call cpu_time(tcpuend)
        call system_clock(tclockend)

        print *, 'Cputime: ', tcpuend-tcpustart, 's'
        print *, 'Wallclock time: ', real(tclockend - tclockstart) / clock_rate, 's'
	print *
    end if

    call mpi_barrier(MPI_COMM_WORLD, ierr)
    call mpi_finalize(ierr)
end program vdw_gas
