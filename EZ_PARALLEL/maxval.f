!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAXVAL. Returns the largest value in a DOUBLE PRECISION matrix.
!
! ARGUMENTS: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - matrix: The array from which we want the largest value (DOUBLE PRECISION,
! DIMENSION(dim1_len, dim2_len)).
! - max_val: The DOUBLE PRECISION variable in which the maximum value will be
! stored (DOUBLE PRECISION).
!
! Written By: Jason Turner
! Last Updated: January 14, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE MAXVAL_EZP(dim1_len, dim2_len, matrix, max_val)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& proc_id, &
& proc_count, &
& ierror
DOUBLE PRECISION :: max_val, &
& proc_max_val
DOUBLE PRECISION :: matrix(dim1_len * dim2_len)

CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, proc_count, ierror)

! Check for errors in user input.
CALL ERROR_HANDLING(proc_id)

! Run in serial if there is only one processor.
IF (proc_count .EQ. 1) THEN
  CALL SERIAL(dim1_len, dim2_len, matrix, max_val)
ELSE
  CALL PARALLEL(dim1_len, dim2_len, matrix, max_val)
END IF

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SERIAL. Handles the execution if there is only one processor.
!
! STRUCTURE: 1) Calculate the maximum value of the grid using the native
! MAXVAL function.
!
! VARIABLES: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - matrix: The array frim which we want the largest value (DOUBLE PRECISION,
! DIMENSION(dim1_len, dim2_len)).
! - max_val: The DOUBLE PRECISION variable in which the maximum value will be
! stored (DOUBLE PRECISION).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE SERIAL(dim1_len, dim2_len, matrix, max_val)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len
DOUBLE PRECISION :: max_val
DOUBLE PRECISION :: matrix(dim1_len * dim2_len)

max_val = MAXVAL(matrix)

END SUBROUTINE SERIAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PARALLEL. Handles the execution if there is more than one processor.
!
! STRUCTURE: 1) Calculate the maximum value of the sub-grid using the native
! MAXVAL function.
! 2) Take the maximum value of these sub-grid maximum values using
! MPI_ALLREDUCE.
!
! VARIABLES: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - matrix: The array frim which we want the largest value (DOUBLE PRECISION,
! DIMENSION(dim1_len, dim2_len)).
! - max_val: The DOUBLE PRECISION variable in which the maximum value will be
! stored (DOUBLE PRECISION).
! - ierror: Integer for holding MPI error flag IDs (INTEGER).
! - proc_max_val: Maximum value of the sub-grid. (DOUBLE PRECISION).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE PARALLEL(dim1_len, dim2_len, matrix, max_val)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& ierror
DOUBLE PRECISION :: max_val, &
& proc_max_val
DOUBLE PRECISION :: matrix(dim1_len * dim2_len)

proc_max_val = MAXVAL(matrix)

CALL MPI_ALLREDUCE(proc_max_val, max_val, 1, MPI_DOUBLE_PRECISION, &
& MPI_MAX, MPI_COMM_WORLD, ierror)

END SUBROUTINE PARALLEL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ERROR_HANDLING. Checks for input error and other possible complications.
!
! ERRORS TO CHECK: N/A. (Assume that users will not enter negative matrix
! dimension lengths.)
!
! VARIABLES: - proc_id: Processor ID (INTEGER).
! - proc_err_flag - Flag is triggered if processor encounters an error
! (LOGICAL).
! - global_err_flag - Flag is triggered if any processor encouters an error
! (LOGICAL).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE ERROR_HANDLING(proc_id)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: proc_id
LOGICAL :: proc_err_flag, &
& global_err_flag

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

proc_err_flag = .FALSE.
global_err_flag = .FALSE.

! STOP all processors if any encounter an error.
CALL MPI_ALLREDUCE(proc_err_flag, global_err_flag, 1, MPI_LOGICAL, MPI_LOR, &
& MPI_COMM_WORLD, ierror)

IF (global_err_flag) THEN
  CALL MPI_FINALIZE(ierror)
  STOP
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

END SUBROUTINE ERROR_HANDLING
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE MAXVAL_EZP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
