!
! integrators.f90
! Equilibrium Monte Carlo simulation of the 2D Ising model
! Oriol Miro
!
! Integrators for Molecular Dynamics.
! A 3-dimensional cubic system is assumed.
!

subroutine verlet(part_num, dt, system_size, cutoff, positions, positions_old, forces)
    implicit none

    integer, intent(in) :: part_num
    real, intent(in) :: dt, system_size, cutoff
    real, allocatable, intent(inout) :: positions(:,:), positions_old(:,:)
    real, allocatable, intent(out) :: forces(:,:)
    real, allocatable :: positions_aux(:,:)


    allocate(positions(part_num, 3), positions_old(part_num, 3), positions_aux(part_num, 3))
    call compute_forces(part_num, positions, forces, system_size, cutoff)
    positions_aux = positions
    positions = r*positions - positions_old + forces*dt*dt
    positions_old = positions_aux

end subroutine verlet