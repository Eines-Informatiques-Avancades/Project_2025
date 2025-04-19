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
        subroutine compute_total_kinetic_energy(velocities, total_kinetic_energy)
            use mpi

            implicit none

            real(8), allocatable, intent(in) :: velocities(:, :)
            real(8), intent(out) :: total_kinetic_energy

            integer :: i
            real(8) :: velocity_norm_sq, kinetic_energy


            kinetic_energy = 0
            total_kinetic_energy = 0
            do i = start_part, end_part
                velocity_norm_sq = velocities(i, 1)**2 + velocities(i, 2)**2 + velocities(i, 3)**2
                kinetic_energy = kinetic_energy + 0.5*velocity_norm_sq
            end do

            call mpi_allreduce( &
                kinetic_energy, total_kinetic_energy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr &
            )
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
