!
! initial_conf.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Ricard Rodriguez
!
! Subroutines related to generating the system's initial configuration.
! A 3-dimensional cubic system is assumed.
!

module initial_conf
    use global_vars

    implicit none

    contains
        ! Generate an initial lattice for the system (either SC or BCC).
        subroutine gen_initial_conf(part_density, positions)
            implicit none

            real(8), intent(out) :: part_density
            real(8), allocatable, intent(out) :: positions(:, :)

            integer :: i, j, k, uc_onedim, part_index
            real(8) :: lattice_spacing

            print *, 'Initializing configuration for ', trim(lattice_type), ' lattice...'

            allocate(positions(part_num, 3))

            ! Generate the specified lattice type for the specified number of atoms with
            ! the specified particle density.
            select case (trim(lattice_type))
                ! Simple cubic lattice.
                case("SC")

                    ! Lattice parameters.
                    uc_onedim = int(part_num**(1./3.))
                    part_density = (part_num/system_size**(3.))
                    lattice_spacing = system_size/uc_onedim

                    part_index = 1
                    do i = 1, uc_onedim
                        do j = 1, uc_onedim
                            do k = 1, uc_onedim
                                positions(part_index, 1) = i - 1
                                positions(part_index, 2) = j - 1
                                positions(part_index, 3) = k - 1

                                part_index = part_index + 1
                            end do
                        end do
                    end do

                ! Body-centered cubic lattice.
                case("BCC")

                    ! Lattice parameters.
                    uc_onedim = int((part_num/2)**(1./3.))
                    part_density = (part_num/system_size**(3.))
                    lattice_spacing = system_size/uc_onedim

                    part_index = 1
                    do i = 1, uc_onedim
                        do j = 1, uc_onedim
                            do k = 1, uc_onedim
                                positions(part_index, 1) = i - 1
                                positions(part_index, 2) = j - 1
                                positions(part_index, 3) = k - 1

                                part_index = part_index + 1

                                ! Add central atom
                                positions(part_index, 1) = i - 1 + 0.5
                                positions(part_index, 2) = j - 1 + 0.5
                                positions(part_index, 3) = k - 1 + 0.5

                                part_index = part_index + 1
                            end do
                        end do
                    end do

                case default
                    print *, 'Invalid value for lattice_type: ', lattice_type
                    stop
            end select

            positions(:, :) = lattice_spacing * positions(:, :)
        end subroutine gen_initial_conf

        ! Generate a bimodal distribution of velocities, where particles either have a
        ! velocity of -v_norm or +v_norm, which is assigned randomly.
        ! v_norm absolute value is given by the temperature.
        subroutine gen_velocities_bimodal_distr(velocities)
            implicit none

            real(8), intent(out), allocatable :: velocities(:, :)

            integer :: i, j
            real(8) :: random_num, v_norm

            allocate(velocities(part_num, 3))

            v_norm = sqrt(temperature)

            ! Assign velocities to each coordinate of each particle based on generating
            ! a random number between 0 and 1.
            do i = 1, part_num
                do j = 1, 3
                    call random_number(random_num)

                    if (random_num > 0.5) then
                        velocities(i, j) = v_norm
                    else
                        velocities(i, j) = -v_norm
                    end if
                end do
            end do
        end subroutine
end module initial_conf
