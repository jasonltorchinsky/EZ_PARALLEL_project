!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ZERO_PADDING_GET_SHAPE unit test.
!
! Written By: Jason Turner
! Last Updated: January 24, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM ZERO_PADDING_GET_SHAPE_TEST

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
WRITE(*,*) 'ZERO_PADDING_GET_SHAPE unit test complete. Normal termination.'

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

NAMELIST /test_params/ dim1_len, dim2_len, overlap
OPEN(1000, file = 'NAMELIST')
READ(1000, nml = test_params)
CLOSE(1000)

CALL INIT_MPI_EZP
CALL GET_ID_EZP(proc_id)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, proc_count, ierror)

! Output dimensions of grid using processor 0.
IF (proc_id .EQ. 0_qb) THEN
  PRINT *, 'dim1_len_total: ', dim1_len, ' dim2_len_total: ', dim2_len, &
  & ' overlap: ', overlap
  PRINT *, 'scl_dim1_len_total: ', (3*dim1_len)/2, ' scl_dim2_len_total: ', &
  & (3*dim2_len)/2
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

! Decompose the grid.
CALL DECOMP_GRID_EZP(dim2_len, overlap)
! Each processor prints out its sub-grid size.
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
DO i = 0, proc_count-1
  IF (proc_id .EQ. i) THEN
    PRINT *, 'proc_id: ', proc_id, ' dim_len: ', dim1_len, dim2_len
  ELSE
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  END IF
END DO
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

! Get dimensions of the scaled grid.
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

CALL FIN_MPI_EZP

END SUBROUTINE MAIN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END PROGRAM
