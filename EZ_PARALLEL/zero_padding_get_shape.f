!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ZERO_PADDING_GET_SHAPE. Obtains the shape of the matrix needed for
! zero-padding (the 3/2-rule) for dealiasing Foutier transforms.
!
! ARGUMENTS: - dim1_len, dim2_len: The shape of the matrix (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - scl_dim1_len, scl_dim2_len: The integers that will hold the dimensions
! of the scaled grid (INTEGER).
!
!
! Written By: Jason Turner
! Last Updated: January 24, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE ZERO_PADDING_GET_SHAPE_EZP(dim1_len, dim2_len, overlap, &
& scl_dim1_len, scl_dim2_len)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& overlap, &
& scl_dim1_len, &
& scl_dim2_len, &
& proc_id, &
& proc_count, &
& ierror

CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, proc_count, ierror)

! Check for errors in user input.
CALL ERROR_HANDLING(proc_id)

! Run in serial if there is only one processor.
IF (proc_count .EQ. 1) THEN
  CALL SERIAL(dim1_len, dim2_len, scl_dim1_len, scl_dim2_len)
ELSE
  CALL PARALLEL(dim1_len, dim2_len, overlap, scl_dim1_len, scl_dim2_len, &
  & proc_id, proc_count)
END IF


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SERIAL. Handles the execution if there is only one processor.
!
! STRUCTURE: Scaled the dim1_len and dim2_len by 3/2.
!
! VARIABLES: - dim1_len, dim2_len: The shape of the matrix (INTEGER).
! - scl_dim1_len, scl_dim2_len: The integers that will hold the dimensions
! of the scaled grid (INTEGER).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE SERIAL(dim1_len, dim2_len, scl_dim1_len, scl_dim2_len)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len,  &
& scl_dim1_len, &
& scl_dim2_len

scl_dim1_len = (3 * dim1_len)/2
scl_dim2_len = (3 * dim2_len)/2

END SUBROUTINE SERIAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PARALLEL. Handles the execution if there is more than one processor.
!
! STRUCTURE: 1) Calculate what the sub-grid size would be if we were to divide
! the scaled grid in the way that we divided the original grid.
!
!
! VARIABLES: - dim1_len, dim2_len: The shape of the matrix (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - scl_dim1_len, scl_dim2_len: The integers that will hold the dimensions
! of the scaled grid (INTEGER).
! - proc_id: Processor ID (INTEGER).
! - proc_count: Number of processors (INTEGER).
! - dim2_len_interior: dim2_len of sub-grid not counting the overlapping region
! (INTEGER).
! - i: Counting index for DO loops (INTEGER).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE PARALLEL(dim1_len, dim2_len, overlap, scl_dim1_len, scl_dim2_len, &
& proc_id, proc_count)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& overlap, &
& scl_dim1_len, &
& scl_dim2_len, &
& dim2_len_total, &
& dim2_len_interior, &
& proc_id, &
& proc_count, &
& ierror, &
& i

scl_dim1_len = (3 * dim1_len)/2

! Calculate the size of the sub-grid without the overlap region.
IF (MOD(proc_id, proc_count-1) .EQ. 0) THEN
  dim2_len_interior = dim2_len - overlap
ELSE
  dim2_len_interior = dim2_len - 2 * overlap
END IF

CALL MPI_ALLREDUCE(dim2_len_interior, dim2_len_total, 1, MPI_INTEGER, MPI_SUM, &
  & MPI_COMM_WORLD, ierror)
dim2_len_total = (3 * dim2_len_total)/2

! Divide dim2_len of total grid as evenly as possible.
scl_dim2_len = dim2_len_total/proc_count
  ! Now check if dim1_slice of a given processor needs to take additional row,
  ! i.e., dim2_len % proc_count != 0.
IF (proc_id .LT. MODULO(dim2_len_total, proc_count)) THEN
  scl_dim2_len = scl_dim2_len + 1
END IF

! Now add overlap.
IF (MOD(proc_id, proc_count-1) .EQ. 0) THEN
  scl_dim2_len = scl_dim2_len + overlap
ELSE
  scl_dim2_len = scl_dim2_len + 2 * overlap
END IF

END SUBROUTINE PARALLEL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ERROR_HANDLING. Checks for input error and other possible complications.
!
! ERRORS TO CHECK: N/A. Assume dim_len are all positive.
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

END SUBROUTINE ZERO_PADDING_GET_SHAPE_EZP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
