!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ZERO_PADDING_INV_DBLE_CMPLX unit test.
!
! Written By: Jason Turner
! Last Updated: January 28, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM ZERO_PADDING_INV_DBLE_CMPLX_TEST

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
WRITE(*,*) 'ZERO_PADDING_INV_DBLE_CMPLX unit test complete. Normal termination.'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Main program.
!
! STRUCTURE: 1) Reads the NAMELIST. Initializes MPI, decomposes grid,
! calculates the reference point for each sub-grid, and fills in each sub-grid.
! 2) FFT the test grid and print it.
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
& ierror, &
& scl_dim1_len, &
& scl_dim2_len
REAL(dp) :: dim1_ref, &
& dim2_ref, &
& dim1_spc, &
& dim2_spc
COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: test_grid, &
& scaled_test_grid
CHARACTER(len=14) :: output_format_1, &
output_format_2

NAMELIST /test_params/ dim1_len, dim2_len, overlap, dim1_ref, dim2_ref, &
& dim1_spc, dim2_spc
OPEN(1000, file = 'NAMELIST')
READ(1000, nml = test_params)
CLOSE(1000)

! Create format string for output grid.
WRITE(output_format_1,'(A,I0.8,A)') '(', 2_qb*dim1_len, 'F8.4)'

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
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
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
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
DO i = 0, proc_count-1
  IF (proc_id .EQ. i) THEN
    PRINT *, 'proc_id: ', proc_id, ' dim1_ref: ', dim1_ref, ' dim2_ref: ', &
    & dim2_ref
  ELSE
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  END IF
END DO
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

! Fill the test grid..
ALLOCATE(test_grid(dim1_len, dim2_len))
test_grid = 0.0_dp
DO j = 1, dim2_len
  DO i = 1, dim1_len
    test_grid(i, j) = (dim1_ref + REAL(i-1_qb, dp) * dim1_spc) &
    & + (dim2_ref + REAL(j-1_qb, dp) * dim2_spc)
  END DO
END DO

CALL SHARE_SUBGRID_BOUNDARIES_DBLE_CMPLX_EZP(dim1_len, dim2_len, overlap, &
& test_grid)

! Each processor prints out its test_grid.
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
DO i = 0, proc_count-1
  IF (proc_id .EQ. i) THEN
    PRINT *, 'proc_id: ', proc_id, ' test_grid: '
    DO j = 1, dim2_len
      WRITE(*, output_format_1) test_grid(:,j)
    END DO
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  ELSE
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  END IF
END DO
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

! Scale the test grid.
CALL ZERO_PADDING_GET_SHAPE_EZP(dim1_len, dim2_len, overlap, scl_dim1_len, &
  & scl_dim2_len)
! Each processor prints out its scl_dim_len.
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
DO i = 0, proc_count-1
  IF (proc_id .EQ. i) THEN
    PRINT *, 'proc_id: ', proc_id, ' scl_dim_len: ', scl_dim1_len, scl_dim2_len
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  ELSE
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  END IF
END DO
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

! Pad the grid with zeros.
ALLOCATE(scaled_test_grid(scl_dim1_len, scl_dim2_len))
WRITE(output_format_2,'(A,I0.8,A)') '(', 2_qb*scl_dim1_len, 'F8.4)'
scaled_test_grid = 0.0_dp
CALL ZERO_PADDING_DBLE_CMPLX_EZP(dim1_len, dim2_len, test_grid, &
  & scl_dim1_len, scl_dim2_len, scaled_test_grid, overlap)

! Each processor prints out its scaled_test_grid.
DO i = 0, proc_count-1
  IF (proc_id .EQ. i) THEN
    PRINT *, 'proc_id: ', proc_id, ' scaled_test_grid: '
    DO j = 1_qb, scl_dim2_len
      WRITE(*, output_format_2) scaled_test_grid(:,j)
    END DO
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  ELSE
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  END IF
END DO
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

! Remove the padded zeros.
CALL ZERO_PADDING_INV_DBLE_CMPLX_EZP(dim1_len, dim2_len, test_grid, &
  & scl_dim1_len, scl_dim2_len, scaled_test_grid, overlap)

! Each processor prints out its test_grid.
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
DO i = 0, proc_count-1
  IF (proc_id .EQ. i) THEN
    PRINT *, 'proc_id: ', proc_id, ' test_grid: '
    DO j = 1, dim2_len
      WRITE(*, output_format_1) test_grid(:,j)
    END DO
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  ELSE
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  END IF
END DO
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

CALL FIN_MPI_EZP

END SUBROUTINE MAIN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END PROGRAM
