!
! integrators.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Oriol Miro
!
! Integrators for Molecular Dynamics.
! A 3-dimensional cubic system is assumed.
!

module integrators
    use lj_forces
    use geometry
    use global_vars

    implicit none

    contains
        subroutine verlet(dt, positions, positions_old, velocities, forces)
            implicit none

            real, intent(in) :: dt
            real, allocatable, intent(inout) :: positions(:), positions_old(:), velocities(:)

            real, allocatable :: forces(:), positions_aux(:)

            allocate(positions_aux(part_num))

            positions_aux = positions
            positions = 2*positions - positions_old + forces*dt*dt
            positions_old = positions_aux
            velocities = (positions - positions_old)/dt

            call apply_pbc(positions)
        end subroutine verlet

        subroutine velocity_verlet(dt, positions, velocities, forces)
            implicit none

            real, intent(in) :: dt
            real, allocatable, intent(inout) :: positions(:), velocities(:)

            real, allocatable :: forces(:)

            positions = positions + velocities*dt + 0.5 * forces*dt*dt

            call apply_pbc(positions)

            velocities = velocities + 0.5 * forces*dt
        end subroutine velocity_verlet

        subroutine euler(dt, positions, velocities, forces)
            implicit none

            real, intent(in) :: dt
            real, allocatable, intent(inout) :: positions(:), velocities(:)

            real, allocatable :: forces(:)

            positions = positions + velocities * dt + 0.5 * forces*dt*dt
            velocities = velocities + forces*dt

            call apply_pbc(positions)
        end subroutine euler
end module integrators
