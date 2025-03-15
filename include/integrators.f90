!
! integrators.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Oriol Miro
!
! Integrators for Molecular Dynamics.
! A 3-dimensional cubic system is assumed.
!

subroutine verlet(part_num, dt, system_size, cutoff, positions, positions_old, velocities, lj_potential)
    implicit none

    integer, intent(in) :: part_num
    real, intent(in) :: dt, system_size, cutoff
    real, allocatable, intent(inout) :: positions(:, :), positions_old(:, :), velocities(:, :)

    real,intent(out) :: lj_potential
    real, allocatable :: positions_aux(:, :), forces(:, :)

    call compute_forces(part_num, positions, forces, lj_potential, system_size, cutoff)
    positions_aux = positions
    positions = 2*positions - positions_old + forces*dt*dt
    positions_old = positions_aux

    call apply_pbc(positions, system_size)
end subroutine verlet

subroutine velocity_verlet(dt, part_num, system_size, cutoff, positions, velocities, lj_potential)
    implicit none

    integer, intent(in) :: part_num
    real, intent(in) :: dt, system_size, cutoff
    real, allocatable, intent(inout) :: positions(:, :), velocities(:, :)

    real, intent(out) :: lj_potential
    real, allocatable :: forces(:, :)

    call compute_forces(part_num, positions, forces, lj_potential, system_size, cutoff)
    positions = positions + velocities*dt + 0.5 * forces*dt*dt

    call apply_pbc(positions, system_size)

    call compute_forces(part_num, positions, forces, lj_potential, system_size, cutoff)
    velocities = velocities + 0.5 * forces*dt
end subroutine velocity_verlet

subroutine euler(dt, part_num, system_size, cutoff, positions, velocities, lj_potential)
    implicit none

    integer,intent(in) :: part_num
    real, intent(in) :: dt, system_size, cutoff
    real, allocatable, intent(inout) :: positions(:, :), velocities(:, :)

    real, intent(out) :: lj_potential
    real, allocatable :: forces(:, :)

    call compute_forces(part_num, positions, forces, lj_potential, system_size, cutoff)
    positions = positions + velocities * dt + 0.5 * forces*dt*dt
    velocities = velocities + forces*dt

    call apply_pbc(positions, system_size)
end subroutine euler
