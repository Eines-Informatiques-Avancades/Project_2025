!
! integrators.f90
! Equilibrium Monte Carlo simulation of the 2D Ising model
! Oriol Miro
!
! Integrators for Molecular Dynamics.
! A 3-dimensional cubic system is assumed.
!

subroutine verlet(positions, positions_old, dt, system_size,cutoff, forces, pot)
    implicit none

    real, intent(in) :: dt, system_size, cutoff
    real, intent(inout) :: positions(:,:), positions_old(:,:)
    real, intent(out) :: forces(:,:), pot