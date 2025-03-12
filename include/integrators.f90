!
! integrators.f90
! Molecular Dynamics simulation of 3D Van der Waals gas
! Oriol Miro
!
! Integrators for Molecular Dynamics.
! A 3-dimensional cubic system is assumed.
!

! Verlet integrator, updates positions, positions_old and velocities.

subroutine verlet(part_num, dt, system_size, cutoff, positions, positions_old, velocities)
    implicit none

    integer, intent(in) :: part_num
    real, intent(in) :: dt, system_size, cutoff
    real, allocatable, intent(inout) :: positions(:, :), positions_old(:, :)
    real, allocatable, intent(out) :: velocities(:, :)
    
    real, allocatable :: positions_aux(:, :), forces(:, :)

    allocate( &
        positions(part_num, 3), &
        positions_old(part_num, 3), &
        positions_aux(part_num, 3), &
        velocities(part_num,3) &
    )
    
    call compute_forces(part_num, positions, forces, system_size, cutoff)
    positions_aux = positions
    positions = 2*positions - positions_old + forces*dt*dt
    positions_old = positions_aux

    call apply_pbc(positions, system_size)
end subroutine verlet

! Velocity Verlet integrator, updates positions and velocities.
subroutine velocity_verlet(dt, part_num, system_size, cutoff, positions, velocities)
    implicit none

    integer, intent(in) :: part_num
    real, intent(in) :: dt, system_size, cutoff
    real, allocatable, intent(inout) :: positions(:, :), velocities(:, :)
    
    real, allocatable :: forces(:, :)

    allocate( &
        positions(part_num, 3), &
        velocities(part_num, 3) &
    )

    call compute_forces(part_num, positions, forces, system_size, cutoff)
    positions = positions + velocities*dt + 0.5 * forces*dt*dt
    
    call apply_pbc(positions, system_size)
    
    call compute_forces(part_num, positions, forces, system_size, cutoff)
    velocities = velocities + 0.5 * forces*dt
end subroutine velocity_verlet

! Euler integrator, updates positions and velocities.
subroutine euler(dt, part_num, system_size, cutoff, positions, velocities)
    implicit none

    integer,intent(in) :: part_num
    real, intent(in) :: dt, system_size, cutoff
    real, allocatable, intent(inout) :: positions(:, :), velocities(:, :)
    
    real, allocatable :: forces(:, :)
    
    allocate( &
        positions(part_num, 3), &
        velocities(part_num, 3) &
    )

    call compute_forces(part_num, positions, forces, system_size, cutoff)
    positions = positions + velocities * dt + 0.5 * forces*dt*dt
    velocities = velocities + forces*dt
    
    call apply_pbc(positions, system_size)
    velocities = velocities + forces * dt
end subroutine euler


