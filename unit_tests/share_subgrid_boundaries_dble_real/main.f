!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SHARE_SUBGRID_BOUNDARIES_DBLE_REAL unit test.
!
! Written By: Jason Turner
! Last Updated: January 15, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM SHARE_SUBGRID_BOUNDARIES_DBLE_REAL_UNIT_TEST

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
WRITE(*,*) 'SHARE_SUBGRID_BOUNDARIES_DBLE_REAL unit test complete. ', &
& 'Normal termination.'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Main program.
!
! STRUCTURE: 1) Reads the NAMELIST. Initializes MPI, decomposes grid,
! calculates the reference point for each sub-grid, and fills in each sub-grid.
! 2) Change the interior of each sub-grid, and share them.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE MAIN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER(qb) :: dim1_len, &
& dim2_len, &
& overlap, &
& i, &
& j, &
& proc_id, &
& proc_count, &
& ierror
REAL(dp) :: dim1_ref, &
& dim2_ref, &
& dim1_spc, &
& dim2_spc
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: test_grid(:,:)
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
    test_grid(i,j) = (dim1_ref + REAL(i-1_qb, dp) * dim1_spc) &
    & + (dim2_ref + REAL(j-1_qb, dp) * dim2_spc)
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

! Adjust the interior of the sub-grid.
! Processor 0 has boundary at end of dim2.
IF (proc_id .EQ. 0) THEN
  DO j = 1, dim2_len-overlap
    DO i = 1, dim1_len
      test_grid(i,j) = 2 * (dim1_ref + REAL(i-1_qb, dp) * dim1_spc) &
      & + (dim2_ref + REAL(j-1_qb, dp) * dim2_spc)
    END DO
  END DO
! Last processor has boundary at start of dim2.
ELSE IF (proc_id .EQ. proc_count-1) THEN
  DO j = 1+overlap, dim2_len
    DO i = 1, dim1_len
      test_grid(i,j) = 2 * (dim1_ref + REAL(i-1_qb, dp) * dim1_spc) &
      & + (dim2_ref + REAL(j-1_qb, dp) * dim2_spc)
    END DO
  END DO
! All interior processors have boundary at start and end of dim2.
ELSE
  DO j = 1+overlap, dim2_len-overlap
    DO i = 1, dim1_len
      test_grid(i,j) = 2 * (dim1_ref + REAL(i-1_qb, dp) * dim1_spc) &
      & + (dim2_ref + REAL(j-1_qb, dp) * dim2_spc)
    END DO
  END DO
END IF

! Share sub-grid boundaries.
CALL SHARE_SUBGRID_BOUNDARIES_DBLE_REAL_EZP(dim1_len, dim2_len, overlap, &
& test_grid)

! Each processor prints out its test_grid.
DO i = 0, proc_count-1
  IF (proc_id .EQ. i) THEN
    PRINT *, 'proc_id: ', proc_id, ' output test_grid: '
    DO j = 1, dim2_len
      WRITE(*, output_format) REAL(test_grid(:,j), dp)
    END DO
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
