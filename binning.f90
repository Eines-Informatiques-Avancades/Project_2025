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
    integer(8) :: i, j, k, num_values, max_num_bins, num_bins, bin_size, bin_start, bin_end
    real(8) :: sum_bin_variable, var_variable, sum_var_variable, mean_variable, step_value
    real(8), allocatable :: variable_name(:), avg_bin_variable_name(:)

    print *, 'Binning data analysis for LJ potential data.'

    ! Define the input file name
    input_datafile = 'output/lj_potential.dat'
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

    allocate(variable_name(num_values))

    ! Read data from file
    do i = 1, num_values
        read(4, *) step_value, variable_name(i)  ! step_value is now real(8)
    end do
    close(4)

    output_file = 'output/binning_lj.out'
    open(6, file = output_file)
    write(6, *) '# Bin size, Mean, Variance, stddev'

    ! Determine the maximum amount of bins that can be used
    max_num_bins = 1
    do while (num_values/(2**max_num_bins) > 100)
        max_num_bins = max_num_bins + 1
    end do

    ! Iterate over the possible bin sizes
    do i = 0, max_num_bins - 1
        bin_size = 2**i
        num_bins = num_values/bin_size

        allocate(avg_bin_variable_name(num_bins))

        ! Iterate over all bins (for current bin size)
        do j = 1, num_bins
            bin_start = (j - 1)*bin_size + 1
            bin_end = j*bin_size

            sum_bin_variable = 0

            do k = bin_start, bin_end
                sum_bin_variable = sum_bin_variable + variable_name(k)
            end do

            avg_bin_variable_name(j) = sum_bin_variable/bin_size
        end do

        mean_variable = sum(avg_bin_variable_name)/num_bins

        ! Compute variance
        sum_var_variable = 0
        do j = 1, num_bins
            sum_var_variable = sum_var_variable + (avg_bin_variable_name(j) - mean_variable)**2
        end do

        var_variable = sum_var_variable/real(((int(num_bins, 16) - 1)*int(num_bins, 16)), 8)

        write(6, *) bin_size, mean_variable, var_variable, sqrt(var_variable)

        deallocate(avg_bin_variable_name)
    end do

    deallocate(variable_name)
    close(6)

end program binning

