PROGRAM visualize
  IMPLICIT NONE
  INTEGER :: frame, n, i, ios
  REAL(4), ALLOCATABLE :: x(:), y(:), z(:)
  CHARACTER(len=10) :: t  ! Atom type (limited length)

  ! Open the trajectory file
  OPEN(1, FILE="traj.xyz", STATUS="OLD", ACTION="READ")

  ! Read frames until end of file
  DO 
    ! Read frame number
    READ(1, *, IOSTAT=ios) frame
    IF (ios /= 0) EXIT  ! Stop reading if EOF

    ! Read number of atoms
    READ(1, *, IOSTAT=ios) n
    IF (ios /= 0) EXIT

    ! Allocate arrays
    IF (.NOT. ALLOCATED(x)) ALLOCATE(x(n), y(n), z(n))

    PRINT *, "Frame:", frame
    PRINT *, "Number of atoms:", n

    ! Read atomic data
    DO i = 1, n
      READ(1, *, IOSTAT=ios) t, x(i), y(i), z(i)
      IF (ios /= 0) EXIT
      PRINT *, "Atom + XYZ:", TRIM(t), x(i), y(i), z(i)
    END DO

  END DO

  ! Close file
  CLOSE(1)
  PRINT *, "Finished reading trajectory."

END PROGRAM visualize

