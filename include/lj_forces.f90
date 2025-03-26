!
! lj_forces.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Alejandro Díaz, Ricard Rodríguez
!
! Compute forces derived from a Lennard-Jones potential.
!

module lj_forces
    use geometry

    implicit none

    contains
        subroutine compute_forces(positions, forces)
            implicit none

            real, allocatable, intent(in) :: positions(:)
            real, allocatable, intent(out) :: forces(:)

            integer :: i, j
            real :: r, f, r_diff

            allocate(forces(part_num))

            ! In reduced units: epsilon = sigma = 1.
            forces = 0.0

            do i = 1, part_num - 1
                do j = i + 1, part_num
                    ! Compute distance between particles i and j.
                    r_diff = positions(i) - positions(j)

                    ! Apply PBC
                    r_diff = pbc(r_diff, system_size)

                    r = abs(r_diff)

                    if (r < cutoff .and. r > 0.0) then

                        ! Force related to a Lennard-Jones potential.
                        f = 48.0 / (r**14) - 24.0 / (r**8)

                        ! Update forces.
                        forces(i) = forces(i) + f * r_diff
                        forces(j) = forces(j) - f * r_diff
                    end if
                end do
            end do
        end subroutine compute_forces
end module lj_forces
