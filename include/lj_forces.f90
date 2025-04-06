!
! lj_forces.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Alejandro Díaz, Ricard Rodríguez
!
! Compute forces derived from a Lennard-Jones potential.
!

module lj_forces
    use global_vars
    use geometry

    implicit none

    contains
        subroutine compute_forces(positions, forces, lj_potential)
            implicit none

            real(8), allocatable, intent(in) :: positions(:, :)
            real(8), allocatable, intent(out) :: forces(:, :)
            real(8), intent(out) :: lj_potential

            integer :: i, j, k
            real(8) :: r, f
            real(8) :: r_vec(3)

            allocate(forces(part_num, 3))

            ! In reduced units: epsilon = sigma = 1.
            forces = 0.0
            lj_potential = 0.0

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
                        ! Lennard-Jones potential
                        lj_potential = lj_potential + 4.0*(1.0/(r**12) - 1.0/(r**6))

                        ! Force related to a Lennard-Jones potential.
                        f = 48.0/(r**14) - 24.0/(r**8)

                        ! Update forces.
                        forces(i, :) = forces(i, :) + f * r_vec(:)
                        forces(j, :) = forces(j, :) - f * r_vec(:)
                    end if
                end do
            end do
        end subroutine compute_forces
end module lj_forces
