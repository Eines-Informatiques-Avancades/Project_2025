subroutine gen_initial_conf(lattice_type, system_size, part_num, part_density, positions)
    implicit none

    integer, intent(in) ::  part_num
    real, intent(in) :: part_density
    character(5), intent(in) :: lattice_type
    real, intent(out) :: system_size
    real, allocatable, intent(out) :: positions(:, :)

end subroutine gen_initial_conf
