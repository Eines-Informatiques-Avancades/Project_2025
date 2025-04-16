program mpi_binning
   use mpi
   implicit none

   !--------------------------------------------------------------------
   ! 1. Global declarations and parameters
   !--------------------------------------------------------------------
   integer, parameter :: MAX_COLS = 200
   character(len=64), dimension(MAX_COLS) :: col_names
   integer :: n_cols

   ! MPI variables
   integer :: ierr, ierror
   integer :: my_rank, nprocs

   ! We use MPI_OFFSET_KIND for variables related to long offsets
   integer(MPI_OFFSET_KIND) :: file_size_offset, offset

   ! MPI file handler and status
   integer :: fh
   integer :: status(MPI_STATUS_SIZE)

   ! Times
   real(8) :: t_total_start, t_total_end
   real(8) :: t_read_start, t_read_end
   real(8) :: t_binning_start, t_binning_end

   ! File names
   character(len=256) :: input_file, output_file
   integer :: unit_out

   ! Buffers to read the header and the second line
   character(len=1), dimension(1024) :: header_buf, row_buf
   character(len=1024) :: header_line, row_line
   integer :: header_length, row_length

   ! Number of rows and data to distribute
   integer(MPI_OFFSET_KIND) :: n_rows, local_rows, start_row
   integer(MPI_OFFSET_KIND) :: base, rem, buf_len

   ! Local buffer read
   character(len=1), allocatable :: local_buffer(:)

   ! Vectors for data
   real(8), allocatable :: local_vec(:), global_vec(:)

   ! For MPI_Allgatherv
   integer, allocatable :: recv_counts(:), displs(:)
   integer :: local_count

   ! Auxiliary indices
   integer :: i, j, pos, len_last

   ! To extract columns
   real(8), allocatable :: col_vec(:)

   ! We initialize time variables to avoid warnings for uninitialized use
   t_read_start = 0.0d0
   t_read_end   = 0.0d0

   !--------------------------------------------------------------------
   ! 2. Start of the program
   !--------------------------------------------------------------------
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

   t_total_start = MPI_Wtime()

   input_file  = 'thermodynamics.dat'
   output_file = 'binning_' // trim(input_file)

   !--------------------------------------------------------------------
   ! 3. Open file in parallel (MPI I/O)
   !--------------------------------------------------------------------
   call MPI_File_open(MPI_COMM_WORLD, trim(input_file), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
   if (ierr /= 0) then
      if (my_rank == 0) then
         print *, "ERROR: Could not open the file:", trim(input_file)
      end if
      ! errorcode = ierr, the third argument is ierror (different)
      call MPI_Abort(MPI_COMM_WORLD, ierr, ierror)
   end if

   !--------------------------------------------------------------------
   ! 4. Read the header (only rank 0)
   !--------------------------------------------------------------------
   if (my_rank == 0) then
      offset = int(0, kind=MPI_OFFSET_KIND)
      call MPI_File_seek(fh, offset, MPI_SEEK_SET, ierr)
      call MPI_File_read(fh, header_buf, 1024, MPI_CHARACTER, status, ierr)

      header_line = array_to_string(header_buf)
      header_length = index(header_line, char(10))
      if (header_length == 0) header_length = len_trim(header_line)

      print *, "[Rank 0 Debug] header_length =", header_length
      print *, "[Rank 0 Debug] Header = |", trim(header_line), "|"
   end if

   ! Broadcast the header line and its length to everyone
   call MPI_Bcast(header_line, 1024, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(header_length, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   call parse_header(header_line(1:header_length), col_names, n_cols)

   ! Adjust the last column name, removing a character if needed
   if (n_cols > 0) then
      len_last = len_trim(col_names(n_cols))
      if (len_last > 1) then
         col_names(n_cols) = col_names(n_cols)(1:len_last-1)
      else
         col_names(n_cols) = ''
      end if
   end if

   !--------------------------------------------------------------------
   ! 5. Read the second line (only rank 0) to determine row_length
   !--------------------------------------------------------------------
   if (my_rank == 0) then
      offset = int(header_length, kind=MPI_OFFSET_KIND)
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
   call MPI_Bcast(row_length, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   !--------------------------------------------------------------------
   ! 6. Get the total size and calculate n_rows
   !--------------------------------------------------------------------
   call MPI_File_get_size(fh, file_size_offset, ierr)

   if (row_length <= 0) then
      if (my_rank == 0) then
         print *, "ERROR: row_length <= 0. Check the file."
      end if
      call MPI_Abort(MPI_COMM_WORLD, 111, ierror)
   end if

   ! n_rows = (file_size_offset - header_length) / row_length  (in 64 bits)
   n_rows = ( file_size_offset - int(header_length, MPI_OFFSET_KIND) ) / &
            int(row_length, MPI_OFFSET_KIND)

   if (n_rows <= 0) then
      if (my_rank == 0) then
         print *, "ERROR: No data rows were found."
      end if
      call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
   end if

   if (my_rank == 0) then
      print *, "[Rank 0 Debug] file_size_offset =", file_size_offset, &
               " => n_rows =", n_rows, " n_cols =", n_cols
   end if

   !--------------------------------------------------------------------
   ! 7. Distribute rows among processes (base+remainder method)
   !--------------------------------------------------------------------
   ! base = n_rows / nprocs  => nprocs must be in 64 bits
   base = n_rows / int(nprocs, MPI_OFFSET_KIND)
   rem  = mod(n_rows, int(nprocs, MPI_OFFSET_KIND))

   if (my_rank < rem) then
      local_rows = base + 1_MPI_OFFSET_KIND
      start_row  = int(my_rank, MPI_OFFSET_KIND)*(base + 1_MPI_OFFSET_KIND)
   else
      local_rows = base
      start_row  = rem*(base + 1_MPI_OFFSET_KIND) + &
                   (int(my_rank, MPI_OFFSET_KIND) - rem)*base
   end if

   offset  = int(header_length, MPI_OFFSET_KIND) + &
             start_row * int(row_length, MPI_OFFSET_KIND)

   buf_len = local_rows * int(row_length, MPI_OFFSET_KIND)

   if (offset + buf_len > file_size_offset) then
      buf_len = file_size_offset - offset
      local_rows = buf_len / int(row_length, MPI_OFFSET_KIND)
   end if

   print *, "[Rank", my_rank, " Debug] offset =", offset, &
            " buf_len =", buf_len, " local_rows =", local_rows, &
            " row_length =", row_length

   !--------------------------------------------------------------------
   ! 8. Read block assigned to each process
   !--------------------------------------------------------------------
   if (local_rows > 0) then
      allocate(local_buffer( int(buf_len) ))
      t_read_start = MPI_Wtime()

      call MPI_File_seek(fh, offset, MPI_SEEK_SET, ierr)
      call MPI_File_read(fh, local_buffer, int(buf_len), MPI_CHARACTER, status, ierr)

      t_read_end = MPI_Wtime()
   else
      allocate(local_buffer(0))
   end if

   ! Close the file (collective)
   call MPI_File_close(fh, ierr)

   ! Convert the read data to a real vector (row-major)
   allocate(local_vec( int(local_rows)*n_cols ))
   do i = 1, int(local_rows)
      pos = (i-1)*row_length + 1
      call parse_fixed_line(local_buffer(pos:pos+row_length-1), &
                            local_vec((i-1)*n_cols+1 : i*n_cols), n_cols)
   end do
   if (allocated(local_buffer)) deallocate(local_buffer)

   !--------------------------------------------------------------------
   ! 9. Reconstruct the global vector with MPI_Allgatherv
   !--------------------------------------------------------------------
   local_count = int(local_rows)*n_cols

   allocate(recv_counts(nprocs), displs(nprocs))
   call MPI_Allgather(local_count, 1, MPI_INTEGER, &
                      recv_counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

   displs(1) = 0
   do i = 2, nprocs
      displs(i) = displs(i-1) + recv_counts(i-1)
   end do

   allocate(global_vec( int(n_rows)*n_cols ))
   call MPI_Allgatherv(local_vec, local_count, MPI_DOUBLE_PRECISION, &
                       global_vec, recv_counts, displs, MPI_DOUBLE_PRECISION, &
                       MPI_COMM_WORLD, ierr)
   if (allocated(local_vec)) deallocate(local_vec)

   !--------------------------------------------------------------------
   ! 10. Extract columns and perform binning
   !--------------------------------------------------------------------
   if (my_rank == 0) then
      open(newunit=unit_out, file=trim(output_file), status='replace', action='write')
      write(unit_out, '(A)') "# Binning results for multicolumn data"
      write(unit_out, '(A)') "# Format: ColName, BinSize, Mean, Variance, StdDev"
   end if

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   t_binning_start = MPI_Wtime()

   do j = 2, n_cols
      allocate(col_vec( int(n_rows) ))
      do i = 1, int(n_rows)
         col_vec(i) = global_vec( (i-1)*n_cols + j )
      end do

      call do_binning_for_column(col_vec, int(n_rows), col_names(j), &
                                 unit_out, my_rank)
      deallocate(col_vec)
   end do

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   t_binning_end = MPI_Wtime()

   if (my_rank == 0) then
      close(unit_out)
   end if

   !--------------------------------------------------------------------
   ! 11. Print times and finish
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

   if (allocated(global_vec)) deallocate(global_vec)
   if (allocated(recv_counts)) deallocate(recv_counts)
   if (allocated(displs))     deallocate(displs)

   call MPI_Finalize(ierr)

contains

   !---------------------------------------------------------------------
   ! array_to_string: converts an array of chars (len=1) to a Fortran string
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
   ! remove_cr: removes CR (ASCII 13) from a string.
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
   ! parse_header: extracts column names from the first line.
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
   ! parse_fixed_line: parses a fixed-length line (array of chars) into a vector of reals
   !---------------------------------------------------------------------
   subroutine parse_fixed_line(line_arr, arr, n_cols)
      implicit none
      character(len=1), dimension(:), intent(in) :: line_arr
      real(8), intent(out) :: arr(n_cols)
      integer, intent(in)  :: n_cols
      character(len=:), allocatable :: cleaned_line

      cleaned_line = remove_cr(trim(array_to_string(line_arr)))
      read(cleaned_line, *) arr
   end subroutine parse_fixed_line

   !---------------------------------------------------------------------
   ! do_binning_for_column: Performs binning for 1 column
   !---------------------------------------------------------------------
   subroutine do_binning_for_column(col_data, n_values, col_name, unit_out, my_rank)
      use mpi
      implicit none
      real(8), intent(in) :: col_data(:)
      integer, intent(in) :: n_values
      character(len=*), intent(in) :: col_name
      integer, intent(in) :: unit_out, my_rank

      integer(8) :: j, k
      integer(8) :: max_num_bins, num_bins, bin_size
      integer(8) :: bin_start, bin_end
      real(8) :: sum_bin, mean_global, var_value, sum_diff
      real(8), allocatable :: avg_bin(:)
      integer(8) :: level

      ! Determine binning levels (powers of 2) with at least ~100 points
      max_num_bins = 1
      do while (n_values / (2**max_num_bins) > 100)
         max_num_bins = max_num_bins + 1
      end do

      do level = 0, max_num_bins - 1
         bin_size = 2**level
         num_bins = n_values / bin_size
         if (num_bins < 2) cycle

         allocate(avg_bin(num_bins))

         ! Mean in each bin
         do j = 1, num_bins
            bin_start = (j - 1)*bin_size + 1
            bin_end   = j*bin_size
            sum_bin = 0.0d0
            do k = bin_start, bin_end
               sum_bin = sum_bin + col_data(k)
            end do
            avg_bin(j) = sum_bin / bin_size
         end do

         ! Calculate global mean
         mean_global = sum(avg_bin) / num_bins

         ! Calculate variance (sum of differences^2)
         sum_diff = 0.0d0
         do j = 1, num_bins
            sum_diff = sum_diff + (avg_bin(j) - mean_global)**2
         end do

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

