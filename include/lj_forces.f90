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
        subroutine compute_forces(positions_x, positions_y, positions_z, forces_x, forces_y, forces_z)
            implicit none

            real, allocatable, intent(in) :: positions_x(:), positions_y(:), positions_z(:)
            real, allocatable :: positions(:, :)
            real, allocatable, intent(out) :: forces_x(:), forces_y(:), forces_z(:) 
            real, allocatable :: forces(:, :)

            integer :: i, j, k, l, m
            real :: r, f
            real :: r_vec(3)

            allocate(forces(part_num, 3))
            allocate(forces_x(part_num), forces_y(part_num), forces_z(part_num))
            allocate(positions(part_num, 3))

            do l = 1, part_num
                positions(l,1) = positions_x(l) 
                positions(l,2) = positions_y(l) 
                positions(l,3) = positions_z(l) 
            end do

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

            do m = 1, part_num
                forces_x(m) = forces(m,1)
                forces_y(m) = forces(m,2)
                forces_z(m) = forces(m,3)
            end do
        end subroutine compute_forces
end module lj_forces
