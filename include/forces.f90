!
! forces.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Alejandro Díaz, Ricard Rodríguez
!
! Compute forces derived from a Lennard-Jones potential.
!

subroutine compute_forces(part_num, positions, forces, system_size, cutoff)
    implicit none

    integer, intent(in) :: part_num
    real, intent(in) :: system_size, cutoff
    real, allocatable, intent(in) :: positions(:, :)
    real, allocatable, intent(out) :: forces(:, :)

    integer :: i, j, k
    real :: r, f
    real :: r_vec(3)

    allocate(forces(part_num, 3))

    ! In reduced units: epsilon = sigma = 1.
    forces = 0.0

    do i = 1, part_num - 1
        do j = i + 1, part_num
            ! Compute distance between particles i and j.
            r_vec(:) = positions(i, :) - positions(j, :)

            ! Apply PBC
            do k = 1, 3
                r_vec(k) = pbc(r_vec(k), system_size)
            end do

            r = sqrt(dot_product(r_vec, r_vec))

            if (r < cutoff .and. r > 0.0) then
                ! Force related to a Lennard-Jones potential.
                f = 48.0 / (r**14) - 24.0 / (r**8)

                ! Update forces.
                forces(i, :) = forces(i, :) + f * r_vec(:)
                forces(j, :) = forces(j, :) - f * r_vec(:)
            end if
        end do
    end do
end subroutine compute_forces
