!
! subroutines.f90
! Equilibrium Monte Carlo simulation of the 2D Ising model
! Ricard Rodriguez i Dargallo
!
! This is part of the "Molecular Dynamics Simulation of a Van der Waals Gas"
! program.
!

module subroutines
    implicit none

    contains
        include 'include/initial_conf.f90'
        include 'include/geometry.f90'
        include 'include/io.f90'
        include 'include/forces.f90'
        include 'include/integrators.f90'
        include 'include/thermostat.f90'
        include 'include/thermodynamics.f90'
        include 'include/post_trajectory_analysis.f90'

end module subroutines
