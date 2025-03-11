! this module is used for the storage of the n (number of atoms/particles), t(simulation time), frame  and positions information.
! also it unsures that these variables can be accessed by the subroutine  
MODULE VAR
  IMPLICIT NONE
  INTEGER :: n, t
  INTEGER, ALLOCATABLE :: frame(:)  
  REAL, ALLOCATABLE :: x(:,:), y(:,:), z(:,:)  
  REAL :: d, a, v
END MODULE VAR

! test main program
PROGRAM MAIN
  USE VAR
  IMPLICIT NONE

  ! n: number of atoms/particles = d * v
  ! t: simulation time
  ! d: density = n/v
  ! a: lenght of cube
  ! v: volume
  t = 5  
  d = 0.8
  a = 1.55
  v = a**3
  ! NINT is necessary to truncate the nearest integer number
  n = NINT(d * v)
  ! allocate memory for x y z variable that is labeled following the n(atom1,atom2,...) and the time frame (time1,time2,...)
  ALLOCATE(x(n, t), y(n, t), z(n, t), frame(t))  

  ! reading trajectory file (xyz format)
  CALL READING()

  ! test to proff that the program have access to all the stored xyz information
  !CALL TEST1()

  ! second test to access certain frame xyz information
  !CALL TEST2()
  print*,n
  CALL CALCULATE_RDF()

  CALL CALCULATE_RMSD()

  ! at the end release the memory for these variables storage.
  DEALLOCATE(x, y, z, frame)
END PROGRAM MAIN

SUBROUTINE READING()
  USE VAR
  IMPLICIT NONE
  INTEGER :: i, j, ios, f

  ! open trajectory file
  OPEN(1, FILE="traj.xyz", STATUS="OLD", ACTION="READ")

  ! read each lines information 
  ! i.e. at time j=1 or frame=1, read n atoms xyz information, then for j=2,3,... the same
  DO j = 1, t
    DO i = 1, n
      ! this 'xyz' format is： frame x y z
      ! these positions are vectorized following the the atom number (i=atom1,atom2,...) and the frame that is registered (j=1,2,...)
      READ(1, *, IOSTAT=ios) f, x(i, j), y(i, j), z(i, j)

      ! store the frame information (it's enough with register first atom frame)
      IF (i == 1) THEN
        frame(j) = f
      END IF
    END DO
  END DO

  ! close the file
  CLOSE(1)
  PRINT *, "Finished reading trajectory."

END SUBROUTINE READING

! test to proff that the program have access to all the stored xyz information
SUBROUTINE TEST1()
  USE VAR
  IMPLICIT NONE
  INTEGER :: i, j

  PRINT *, "Processing stored coordinates:"
  DO j = 1, t
    PRINT *, "Frame:", frame(j)  
    DO i = 1, n
      PRINT *,": (", x(i, j), y(i, j), z(i, j), ")"
    END DO
  END DO
END SUBROUTINE TEST1

! second test to access certain frame xyz information
SUBROUTINE TEST2()
  USE VAR
  IMPLICIT NONE
  INTEGER :: i, j

  PRINT *, "Testing specific frame data:"
  DO j = 1, t
    IF (frame(j) == 3) THEN  ! let's test for frame=3
      PRINT *, "Frame:", frame(j)
      DO i = 1, n
        PRINT*, 'x', x(i, j), 'y', y(i, j), 'z', z(i, j)
      END DO
    END IF
  END DO
END SUBROUTINE TEST2

! with stored datas, compute rdf and rmsd

! rdf
SUBROUTINE CALCULATE_RDF()
  USE VAR
  IMPLICIT NONE
  INTEGER :: i, j, k, fi
  REAL(4), PARAMETER :: rmax = 1.0  ! maximum radius
  REAL(4), PARAMETER :: dr = 0.1      ! rangesteps
  INTEGER :: bins                     ! bin number
  REAL(4), ALLOCATABLE :: rdf(:), r_values(:) 
  REAL(4) :: r, dx, dy, dz, dv, density  ! dx,dy,dz positions variations
  INTEGER :: bi                          ! bi bv are respectively bin index and spherical volume
  REAL(4) :: volume

  bins = INT(rmax / dr)
  ALLOCATE(rdf(bins), r_values(bins))
  rdf = 0.0  

 
  DO fi = 1, t       ! loop for all frames
    DO i = 1, n-1    ! loop for all atoms  i.e. i=1 j=2 dx=x1-x2 ...
      DO j = i+1, n  
        ! x y z variation between frames to calculate r
        dx = x(j, fi) - x(i, fi)
        dy = y(j, fi) - y(i, fi)
        dz = z(j, fi) - z(i, fi)
        r = SQRT(dx**2 + dy**2 + dz**2)     ! r = sqrt(dx^2 + dy^2 + dz^2)
        PRINT *, 'r:',r                     ! check the calculated r
        ! consider  sphere to calculate rdf
        IF (r < rmax) THEN
          bi = INT(r / dr) + 1  ! there a different portions/zones of spheres, each zone/portions labeled as bi=1,2,3,...
          rdf(bi) = rdf(bi) + 1 ! bi correspond to the zone of sphere that this r belongs and we will increase the size to accumulate the rdf
        END IF
      END DO
    END DO
  END DO

  volume = (4.0/3.0) * 3.1415926 * rmax**3   ! this volume is not the entire cubic volume, it's refering the rdf related spherical volume
  density = n  / volume             ! therefore, it's necessary to set a maximum radius 

  DO k = 1, bins
    r_values(k) = k * dr
    dv = 4.0 * 3.14159 * r_values(k)**2 * dr  ! volume between r to r + dr to normalize the rdf
    rdf(k) = rdf(k) / (density * n * dv )  ! normalize the rdf
  END DO

  ! store the rdf information
  OPEN(2, FILE="rdf_data.txt", STATUS="REPLACE")
  DO k = 1, bins
    WRITE(2, *) r_values(k), rdf(k)
  END DO
  CLOSE(2)
  PRINT *, "RDF calculation completed and saved to rdf_data.txt"
  DEALLOCATE(rdf, r_values)

END SUBROUTINE CALCULATE_RDF

! rmsd
SUBROUTINE CALCULATE_RMSD()
  USE VAR
  IMPLICIT NONE
  INTEGER :: i, j
  REAL(4), ALLOCATABLE :: rmsd(:)
  REAL(4) :: dx, dy, dz, sum_sq

 
  ALLOCATE(rmsd(t))

  
  DO j = 1, t   ! for all time frame
    sum_sq = 0.0  ! summation using the acumulation 
    DO i = 1, n   
      dx = x(i, j) - x(i, 1)  ! here difference of positions is between the xj and x1 (as reference and changeble)
      dy = y(i, j) - y(i, 1)
      dz = z(i, j) - z(i, 1)
      sum_sq = sum_sq + (dx**2 + dy**2 + dz**2)
    END DO
    rmsd(j) = SQRT(sum_sq / n)    !   rmsd=sqrt(summation(r-r')^2 / N)
  END DO

  ! store the obtained rmsd
  OPEN(2, FILE="rmsd_data.txt", STATUS="REPLACE")
  DO j = 1, t
    WRITE(2, *) frame(j), rmsd(j)
  END DO
  CLOSE(2)

  PRINT *, "RMSD calculation completed and saved to rmsd_data.txt"

  DEALLOCATE(rmsd)

END SUBROUTINE CALCULATE_RMSD




