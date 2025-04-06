!
! thermodynamics.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Ricard Rodriguez
!
! Compute different thermodynamical variables of the system.
!

module thermodynamics
    use global_vars

    implicit none

    contains
        ! Compute the system's kinetic energy in an specific instant, given a matrix of
        ! velocities.
        subroutine compute_total_kinetic_energy(velocities, kinetic_energy)
            implicit none

            real(8), allocatable, intent(in) :: velocities(:, :)
            real(8), intent(out) :: kinetic_energy

            integer :: i
            real(8) :: velocity_norm_sq, kinetic_energy_part

            kinetic_energy = 0
            do i = 1, part_num
                velocity_norm_sq = velocities(i, 1)**2 + velocities(i, 2)**2 + velocities(i, 3)**2
                kinetic_energy_part = 0.5*velocity_norm_sq
                kinetic_energy = kinetic_energy + kinetic_energy_part
            end do
        end subroutine compute_total_kinetic_energy

        ! Compute the system's instantaneous temperature, given the system's particle
        ! number and the system's total kinetic energy.
        function instantaneous_temperature(kinetic_energy) result(temperature_inst)
            implicit none

            integer :: n_f
            real(8) :: kinetic_energy, temperature_inst

            n_f = 3*part_num - 3
            temperature_inst = (2*kinetic_energy)/n_f
        end function
end module thermodynamics
