!=======================================================================
! mpi_binning_parallel.f90
!
! MPI program that:
!   1) Reads in parallel (MPI I/O) a text file where:
!      - The first line is the header (column names).
!      - From the second line onward, each row has a fixed length.
!   2) Reconstructs the global matrix (stored as a 1D vector in
!      row-major order) using MPI_Allgatherv.
!   3) Extracts each column and performs parallel binning,
!      writing results only on rank 0.
!
! Compile with: mpif90 binning_parallel.f90 -o binning_parallel
! Run with: mpirun -np 2 ./binning_parallel
!=======================================================================
program mpi_binning
   use mpi
   implicit none

   !--------------------------------------------------------------------
   ! 1. Global declarations and parameters
   !--------------------------------------------------------------------
   integer, parameter :: MAX_COLS = 200
   character(len=64), dimension(MAX_COLS) :: col_names
   integer :: n_cols

   integer :: ierr, my_rank, nprocs
   integer(MPI_OFFSET_KIND) :: file_size_offset, offset
   integer :: fh                 ! MPI file handle
   integer :: status(MPI_STATUS_SIZE)

   ! Timing variables
   real(8) :: t_total_start, t_total_end
   real(8) :: t_read_start, t_read_end
   real(8) :: t_binning_start, t_binning_end

   ! Files
   character(len=256) :: input_file, output_file
   integer :: unit_out

   ! For reading header and second line (fixed buffers of 1024)
   character(len=1), dimension(1024) :: header_buf, row_buf
   character(len=:), allocatable :: header_line, row_line
   integer :: header_length, row_length

   ! Number of rows and distribution parameters
   integer :: n_rows, local_rows, start_row
   integer :: base, rem

   ! Local buffer (1D character vector)
   integer :: buf_len
   character(len=1), allocatable :: local_buffer(:)

   ! Data: stored in a 1D vector in row-major order.
   real(8), allocatable :: local_vec(:), global_vec(:)

   ! For MPI_Allgatherv
   integer, allocatable :: recv_counts(:), displs(:)
   integer :: local_count

   ! Auxiliary indices
   integer :: i, j, pos, len_last

   ! Auxiliary variable for extracting each column
   real(8), allocatable :: col_vec(:)

   !--------------------------------------------------------------------
   ! 2. Program start
   !--------------------------------------------------------------------
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
   t_total_start = MPI_Wtime()

   input_file  = 'thermodynamics.dat'
   output_file = 'binning_' // trim(input_file)

   !--------------------------------------------------------------------
   ! 3. Open file in parallel mode (MPI I/O)
   !--------------------------------------------------------------------
   call MPI_File_open(MPI_COMM_WORLD, trim(input_file), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
   if (ierr /= 0) then
      if (my_rank == 0) print *, "ERROR: Could not open the file:", trim(input_file)
      call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
   end if

   !--------------------------------------------------------------------
   ! 4. Read the header (rank 0 only)
   !--------------------------------------------------------------------
   if (my_rank == 0) then
      offset = 0_MPI_OFFSET_KIND
      call MPI_File_seek(fh, offset, MPI_SEEK_SET, ierr)
      call MPI_File_read(fh, header_buf, 1024, MPI_CHARACTER, status, ierr)
      header_line = array_to_string(header_buf)
      header_length = index(header_line, char(10))
      if (header_length == 0) header_length = len_trim(header_line)
      print *, "[Rank 0 Debug] header_length =", header_length
      print *, "[Rank 0 Debug] Header = |", trim(header_line), "|"
   end if
   if (my_rank /= 0) then
      allocate(character(len=1024) :: header_line)
   end if
   call MPI_Bcast(header_line, 1024, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(header_length, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call parse_header(header_line(1:header_length), col_names, n_cols)

if (n_cols > 0) then
   ! Obtain the current length without trailing spaces
   len_last = len_trim(col_names(n_cols))
   if (len_last > 0) then
      col_names(n_cols) = col_names(n_cols)(1:len_last-1)
   end if
end if
   !--------------------------------------------------------------------
   ! 5. Read the second line (rank 0 only) to determine row_length
   !--------------------------------------------------------------------
   if (my_rank == 0) then
      offset = header_length
      call MPI_File_seek(fh, offset, MPI_SEEK_SET, ierr)
      call MPI_File_read(fh, row_buf, 1024, MPI_CHARACTER, status, ierr)
      row_line = array_to_string(row_buf)
      row_length = index(row_line, char(10))
      if (row_length == 0) then
         row_length = len_trim(row_line)
      end if
      print *, "[Rank 0 Debug] row_line = |", trim(row_line), "|"
      print *, "[Rank 0 Debug] row_length =", row_length
   end if
   if (my_rank /= 0) then
      allocate(character(len=1024) :: row_line)
   end if
   call MPI_Bcast(row_length, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   !--------------------------------------------------------------------
   ! 6. Get the total file size and calculate n_rows (data rows)
   !--------------------------------------------------------------------
   call MPI_File_get_size(fh, file_size_offset, ierr)
   if (row_length <= 0) then
      if (my_rank == 0) print *, "ERROR: row_length <= 0. Check the file."
      call MPI_Abort(MPI_COMM_WORLD, 111, ierr)
   end if
   n_rows = int((file_size_offset - header_length) / row_length)
   if (n_rows <= 0) then
      if (my_rank == 0) print *, "ERROR: No data rows were found."
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
   end if
   if (my_rank == 0) then
      print *, "[Rank 0 Debug] file_size_offset =", file_size_offset, &
               " => n_rows =", n_rows, " n_cols =", n_cols
   end if

   !--------------------------------------------------------------------
   ! 7. Distribute rows among processes (base+rem scheme)
   !--------------------------------------------------------------------
   base = n_rows / nprocs
   rem  = mod(n_rows, nprocs)
   if (my_rank < rem) then
      local_rows = base + 1
      start_row  = my_rank * (base + 1)
   else
      local_rows = base
      start_row  = rem*(base+1) + (my_rank - rem)*base
   end if
   offset  = header_length + start_row * row_length
   buf_len = local_rows * row_length
   if (offset + buf_len > file_size_offset) then
      buf_len    = file_size_offset - offset
      local_rows = buf_len / row_length
   end if
   print *, "[Rank", my_rank, " Debug] offset =", offset, " buf_len =", buf_len, &
            " local_rows =", local_rows, " row_length =", row_length

   !--------------------------------------------------------------------
   ! 8. Each process reads its assigned block (stored in a 1D vector)
   !--------------------------------------------------------------------
   if (local_rows > 0) then
      allocate(local_buffer(buf_len))
      t_read_start = MPI_Wtime()
      call MPI_File_seek(fh, offset, MPI_SEEK_SET, ierr)
      call MPI_File_read(fh, local_buffer, buf_len, MPI_CHARACTER, status, ierr)
      t_read_end = MPI_Wtime()
      call MPI_File_close(fh, ierr)
   else
      allocate(local_buffer(0))
   end if

   ! Convert the read data into a 1D real vector in row-major order.
   allocate(local_vec(local_rows * n_cols))
   do i = 1, local_rows
      pos = (i - 1) * row_length + 1
      call parse_fixed_line(local_buffer(pos:pos+row_length-1), &
                             local_vec((i-1)*n_cols+1:i*n_cols), n_cols)
   end do
   deallocate(local_buffer)

   !--------------------------------------------------------------------
   ! 9. Reconstruct the global vector using MPI_Allgatherv
   !--------------------------------------------------------------------
   local_count = local_rows * n_cols
   allocate(recv_counts(nprocs), displs(nprocs))
   call MPI_Allgather(local_count, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
   displs(1) = 0
   do i = 2, nprocs
      displs(i) = displs(i-1) + recv_counts(i-1)
   end do
   allocate(global_vec(n_rows * n_cols))
   call MPI_Allgatherv(local_vec, local_count, MPI_DOUBLE_PRECISION, &
                       global_vec, recv_counts, displs, MPI_DOUBLE_PRECISION, &
                       MPI_COMM_WORLD, ierr)
   deallocate(local_vec)

   !--------------------------------------------------------------------
   ! 10. Extract each column from the global vector and perform binning
   !--------------------------------------------------------------------
   if (my_rank == 0) then
      open(newunit=unit_out, file=output_file, status='replace', action='write')
      write(unit_out, '(A)') "# Binning results for multicolumn data"
      write(unit_out, '(A)') "# Format: ColName, BinSize, Mean, Variance, StdDev"
   end if

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   t_binning_start = MPI_Wtime()
   do j = 2, n_cols
      allocate(col_vec(n_rows))
      do i = 1, n_rows
         col_vec(i) = global_vec((i-1)*n_cols + j)
      end do
      call do_binning_for_column(col_vec, n_rows, col_names(j), unit_out, my_rank, nprocs)
      deallocate(col_vec)
   end do
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   t_binning_end = MPI_Wtime()

   if (my_rank == 0) then
      close(unit_out)
   end if

   !--------------------------------------------------------------------
   ! 11. Print timings and finalize
   !--------------------------------------------------------------------
   if (my_rank == 0) then
      t_total_end = MPI_Wtime()
      print *, "============================================="
      print *, "[Rank 0] Total time:   ", t_total_end - t_total_start, " sec"
      print *, "[Rank 0] Read time:    ", t_read_end - t_read_start,   " sec"
      print *, "[Rank 0] Binning time: ", t_binning_end - t_binning_start, " sec"
      print *, "[Rank 0] Results in:   ", trim(output_file)
      print *, "============================================="
   end if

   deallocate(global_vec, recv_counts, displs)
   call MPI_Finalize(ierr)

contains
   !---------------------------------------------------------------------
   ! array_to_string: Converts an array of single characters into a Fortran string
   !---------------------------------------------------------------------
   function array_to_string(arr) result(str)
      implicit none
      character(len=1), dimension(:), intent(in) :: arr
      character(len=size(arr)) :: str
      integer :: i
      do i = 1, size(arr)
         str(i:i) = arr(i)
      end do
   end function array_to_string

   !---------------------------------------------------------------------
   ! remove_cr: Removes CR characters (ASCII 13) from a string.
   !---------------------------------------------------------------------
   function remove_cr(str) result(out)
      implicit none
      character(len=*), intent(in) :: str
      character(len=len(str)) :: out
      integer :: i, j
      j = 0
      do i = 1, len(str)
         if (str(i:i) /= char(13)) then
            j = j + 1
            out(j:j) = str(i:i)
         end if
      end do
      out = out(:j)
   end function remove_cr

   !---------------------------------------------------------------------
   ! parse_header: Parses the header (first line) and extracts column names.
   !---------------------------------------------------------------------
   subroutine parse_header(line_in, col_names, n_cols)
      implicit none
      character(len=*), intent(in) :: line_in
      character(len=64), dimension(:), intent(out) :: col_names
      integer, intent(out) :: n_cols
      character(len=1024) :: tmp_line
      integer :: pos, len_line, token_len, i
      character(len=64) :: token

      tmp_line = adjustl(line_in)
      if (len_trim(tmp_line) > 0 .and. tmp_line(1:1) == '#') then
         tmp_line = adjustl(tmp_line(2:))
      end if
      ! Replace commas with spaces
      do i = 1, len(tmp_line)
         if (tmp_line(i:i) == ',') then
            tmp_line(i:i) = ' '
         end if
      end do
      len_line = len_trim(tmp_line)
      n_cols   = 0
      pos      = 1
      do while (pos <= len_line)
         token_len = index(tmp_line(pos:), ' ')
         if (token_len == 0) then
            token = tmp_line(pos:)
            pos = len_line + 1
         else
            token = tmp_line(pos: pos+token_len-2)
            pos = pos + token_len
         end if
         if (len_trim(token) > 0) then
            n_cols = n_cols + 1
            col_names(n_cols) = trim(token)
         end if
      end do
   end subroutine parse_header

   !---------------------------------------------------------------------
   ! parse_fixed_line: Parses a fixed line (array of single characters)
   ! and converts it into a vector of reals (double precision).
   !---------------------------------------------------------------------
   subroutine parse_fixed_line(line_arr, arr, n_cols)
      implicit none
      character(len=1), dimension(:), intent(in) :: line_arr
      real(8), intent(out) :: arr(n_cols)
      integer, intent(in) :: n_cols
      character(len=:), allocatable :: cleaned_line

      cleaned_line = remove_cr(trim(array_to_string(line_arr)))
      read(cleaned_line, *) arr
   end subroutine parse_fixed_line

   !---------------------------------------------------------------------
   ! do_binning_for_column: Performs binning for one column
   !---------------------------------------------------------------------
   subroutine do_binning_for_column(col_data, n_values, col_name, unit_out, my_rank, nprocs)
      use mpi
      implicit none
      real(8), intent(in) :: col_data(:)
      integer, intent(in) :: n_values
      character(len=*), intent(in) :: col_name
      integer, intent(in) :: unit_out, my_rank, nprocs

      integer(8) :: i, j, k
      integer(8) :: max_num_bins, num_bins, bin_size
      integer(8) :: bin_start, bin_end
      real(8) :: sum_bin, mean_bin
      real(8) :: mean_global, var_value, sum_diff
      real(8), allocatable :: avg_bin(:)
      integer(8) :: level

      ! Determine the maximum binning level (powers of 2) ensuring that
      ! each level has at least 100 data points
      max_num_bins = 1
      do while (n_values / (2**max_num_bins) > 100)
         max_num_bins = max_num_bins + 1
      end do

      do level = 0, max_num_bins - 1
         bin_size = 2**level
         num_bins = n_values / bin_size
         if (num_bins < 2) cycle

         allocate(avg_bin(num_bins))
         ! Calculate the mean in each bin
         do j = 1, num_bins
            bin_start = (j - 1)*bin_size + 1
            bin_end   = j*bin_size
            sum_bin = 0.0d0
            do k = bin_start, bin_end
               sum_bin = sum_bin + col_data(k)
            end do
            avg_bin(j) = sum_bin / bin_size
         end do

         ! Calculate the global mean of the binned values
         mean_global = sum(avg_bin) / num_bins

         ! Calculate the sum of squared differences
         sum_diff = 0.0d0
         do j = 1, num_bins
            sum_diff = sum_diff + (avg_bin(j) - mean_global)**2
         end do

         ! Reference formula: var = sum_diff / ((num_bins-1)*num_bins)
         var_value = sum_diff / real((num_bins - 1)*num_bins, 8)
         if (var_value < 0.d0) var_value = 0.d0

         if (my_rank == 0) then
            write(unit_out, '(A, 1X, I12, 1X, F15.8, 1X, F15.8, 1X, F15.8)') &
                trim(col_name), bin_size, mean_global, var_value, sqrt(var_value)
         end if

         deallocate(avg_bin)
      end do
   end subroutine do_binning_for_column

end program mpi_binning

