!
! lj_potentials.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Alejandro Díaz, Ricard Rodríguez
!
! Compute Lennard-Jones potential.
!

module lj_potentials
    use geometry

    implicit none

    contains
        subroutine compute_lj(positions_x, positions_y, positions_z, lj_potential)
            implicit none

            real, allocatable, intent(in) :: positions_x(:), positions_y(:), positions_z(:)
            real, allocatable :: positions(:, :)
            real, intent(out) :: lj_potential

            integer :: i, j, k, l
            real :: r
            real :: r_vec(3)

            allocate(positions(part_num, 3))

            do l = 1, part_num
                positions(l,1) = positions_x(l) 
                positions(l,2) = positions_y(l) 
                positions(l,3) = positions_z(l) 
            end do

            ! In reduced units: epsilon = sigma = 1.
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
                        lj_potential = lj_potential + 4.0 * (1.0 / (r**12) - 1.0 / (r**6))
                    end if
                end do
            end do
        end subroutine compute_lj
end module lj_potentials
