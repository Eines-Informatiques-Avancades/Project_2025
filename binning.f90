!
! binning.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Ricard Rodriguez, Joan Serrano
!
! Binning program corresponding to the final project for the Advanced Computer
! Tools subject.
!

program binning
    implicit none

    character(50) :: output_file, input_datafile
    integer(8) :: i, j, k, num_values, step, max_num_bins, num_bins, bin_size, bin_start, bin_end
    real(8) :: sum_bin_energy, var_energy, sum_var_energy, mean_energy, step_value
    real(8), allocatable :: energy(:), avg_bin_energy(:)

    print *, 'Binning data analysis for LJ potential data.'

    ! Define the input file name
    input_datafile = 'lj_potential.dat'
    print *, 'Reading data from ', trim(adjustl(input_datafile))

    open(4, file = input_datafile, status = 'old')

    ! Count the number of values in the file
    num_values = 0
    do
        read(4, *, end=10)
        num_values = num_values + 1
    end do
10  continue
    rewind(4)

    allocate(energy(num_values))

    ! Read data from file
    do i = 1, num_values
        read(4, *) step_value, energy(i)  ! step_value is now real(8)
    end do
    close(4)

    output_file = 'binning_lj.out'
    open(6, file = output_file)
    write(6, *) '# Bin size, Mean energy, Energy variance, Energy stddev'

    ! Determine the maximum amount of bins that can be used
    max_num_bins = 1
    do while (num_values/(2**max_num_bins) > 100)
        max_num_bins = max_num_bins + 1
    end do

    ! Iterate over the possible bin sizes
    do i = 0, max_num_bins - 1
        bin_size = 2**i
        num_bins = num_values/bin_size

        allocate(avg_bin_energy(num_bins))

        ! Iterate over all bins (for current bin size)
        do j = 1, num_bins
            bin_start = (j - 1)*bin_size + 1
            bin_end = j*bin_size

            sum_bin_energy = 0

            do k = bin_start, bin_end
                sum_bin_energy = sum_bin_energy + energy(k)
            end do

            avg_bin_energy(j) = sum_bin_energy/bin_size
        end do

        mean_energy = sum(avg_bin_energy)/num_bins

        ! Compute variance
        sum_var_energy = 0
        do j = 1, num_bins
            sum_var_energy = sum_var_energy + (avg_bin_energy(j) - mean_energy)**2
        end do

        var_energy = sum_var_energy/real(((int(num_bins, 16) - 1)*int(num_bins, 16)), 8)

        write(6, *) bin_size, mean_energy, var_energy, sqrt(var_energy)

        deallocate(avg_bin_energy)
    end do

    deallocate(energy)
    close(6)

end program binning

