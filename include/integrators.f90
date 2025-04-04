!
! integrators.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Oriol Miro
!
! Integrators for Molecular Dynamics.
! A 3-dimensional cubic system is assumed.
!

module integrators
    use global_vars
    use lj_forces
    use geometry

    implicit none

    contains
        subroutine verlet(positions, positions_old, velocities, lj_potential)
            implicit none

            real, allocatable, intent(inout) :: positions(:, :), positions_old(:, :), velocities(:, :)
            real, intent(out) :: lj_potential

            real, allocatable :: positions_aux(:, :), forces(:, :)

            allocate(positions_aux(part_num, 3))

            call compute_forces(positions, forces, lj_potential)
            positions_aux = positions
            positions = 2*positions - positions_old + forces*timestep*timestep
            positions_old = positions_aux
            velocities = (positions - positions_old)/timestep

            call apply_pbc(positions)

            deallocate(positions_aux)
        end subroutine verlet

        subroutine velocity_verlet(positions, velocities, lj_potential)
            implicit none

            real, allocatable, intent(inout) :: positions(:, :), velocities(:, :)
            real, intent(out) :: lj_potential

            real, allocatable :: forces(:, :)

            call compute_forces(positions, forces, lj_potential)
            positions = positions + velocities*timestep + 0.5 * forces*timestep*timestep
            velocities = velocities + 0.5 * forces*timestep

            call apply_pbc(positions)

            call compute_forces(positions, forces, lj_potential)
            velocities = velocities + 0.5 * forces*timestep
        end subroutine velocity_verlet

        subroutine euler(positions, velocities, lj_potential)
            implicit none

            real, allocatable, intent(inout) :: positions(:, :), velocities(:, :)
            real, intent(out) :: lj_potential

            real, allocatable :: forces(:, :)

            call compute_forces(positions, forces, lj_potential)
            positions = positions + velocities*timestep + 0.5 * forces*timestep*timestep
            velocities = velocities + forces*timestep

            call apply_pbc(positions)
        end subroutine euler
end module integrators
