!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MINVAL unit test.
!
! Written By: Jason Turner
! Last Updated: January 15, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM MINVAL_UNIT_TEST

IMPLICIT NONE

! Defines standard integer-, real-precision types.
INCLUDE 'integer_types.h'
INCLUDE 'real_types.h'

REAL(dp) :: start_time, &
& end_time

CALL CPU_TIME(start_time)
CALL MAIN
CALL CPU_TIME(end_time)

WRITE(*,'(A,F10.5,A)') 'Execution time: ', end_time - start_time, '.'
WRITE(*,*) 'MINVAL unit test complete. Normal termination.'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Main program.
!
! STRUCTURE: 1) Reads the NAMELIST. Initializes MPI, decomposes grid,
! calculates the reference point for each sub-grid, and fills in each sub-grid.
! 2) Obtains the minimum value across the grid, and each processor prints out
! out this minimum value.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE MAIN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER(qb) :: dim1_len, &
& dim2_len, &
& overlap, &
& order, &
& i, &
& j, &
& proc_id, &
& proc_count, &
& ierror
REAL(dp) :: dim1_ref, &
& dim2_ref, &
& dim1_spc, &
& dim2_spc
REAL(dp) :: min_val
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: test_grid
REAL(dp), PARAMETER :: pi_dp = 4.0_dp * ATAN(1.0_dp)
CHARACTER(len=15) :: output_format

NAMELIST /test_params/ dim1_len, dim2_len, overlap, dim1_ref, dim2_ref, &
& dim1_spc, dim2_spc
OPEN(1000, file = 'NAMELIST')
READ(1000, nml = test_params)
CLOSE(1000)

! Create format string for output grid.
WRITE(output_format,'(A,I0.8,A)') '(', dim1_len, 'F16.8)'

CALL INIT_MPI_EZP
CALL GET_ID_EZP(proc_id)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, proc_count, ierror)

! Output dimensions of grid using processor 0.
IF (proc_id .EQ. 0_qb) THEN
  PRINT *, 'dim1_len_total: ', dim1_len, ' dim2_len_total: ', dim2_len, &
  & ' overlap: ', overlap
  PRINT *, 'dim1_ref: ', dim1_ref, ' dim2_ref: ', dim2_ref, ' dim1_spc: ', &
  & dim1_spc, ' dim2_spc: ', dim2_spc
  PRINT *, 'order: ', order
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

! Fill in the test grid.
CALL DECOMP_GRID_EZP(dim2_len, overlap)
! Each processor prints out its dim2_len.
DO i = 0, proc_count-1
  IF (proc_id .EQ. i) THEN
    PRINT *, 'proc_id: ', proc_id, ' dim2_len: ', dim2_len
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  ELSE
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  END IF
END DO
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

CALL IDENTIFY_REF_POINT_EZP(dim2_len, dim2_ref, dim2_spc, overlap)
! Each processor prints out its reference point.
DO i = 0, proc_count-1
  IF (proc_id .EQ. i) THEN
    PRINT *, 'proc_id: ', proc_id, ' dim1_ref: ', dim1_ref, ' dim2_ref: ', &
    & dim2_ref
  ELSE
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  END IF
END DO
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

ALLOCATE(test_grid(dim1_len, dim2_len))
test_grid = (0.0_dp, 0.0_dp)
DO j = 1, dim2_len
  DO i = 1, dim1_len
    test_grid(i,j) = SIN((dim1_ref + REAL(i-1_qb, dp) * dim1_spc) * pi_dp)
  END DO
END DO
! Each processor prints out its test_grid.
DO i = 0, proc_count-1
  IF (proc_id .EQ. i) THEN
    PRINT *, 'proc_id: ', proc_id, ' test_grid: '
    DO j = 1, dim2_len
      WRITE(*, output_format) REAL(test_grid(:,j), dp)
    END DO
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  ELSE
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  END IF
END DO
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

CALL MINVAL_EZP(dim1_len, dim2_len, test_grid, min_val)
! Each processor prints out the maximum value of the grid..
DO i = 0, proc_count-1
  IF (proc_id .EQ. i) THEN
    PRINT *, 'proc_id: ', proc_id, ' min_val: ', min_val
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  ELSE
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  END IF
END DO
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

CALL FIN_MPI_EZP

RETURN

END SUBROUTINE MAIN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END PROGRAM
