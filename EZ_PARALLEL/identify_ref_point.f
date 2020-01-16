!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IDENTIFY_REF_POINT. Identifies the reference point of the subgrid, assuming
! a linear spacing of grid points. Since only the second dimension (dim2)
! changes upon decomposing the grid, this is all we need to recalculate.
!
! ARGUMENTS: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - spec_dim2_deriv: The array that will store the matrix used in calculating
! the derivative of the numerical solution along the second dimension of the
! grid (DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len)).
! - order: The order of the derviate desired (INTEGER).
!
! NOTES: - This assumes a linear spacing of grid points along dim2.
!
! Written By: Jason Turner
! Last Updated: January --, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE IDENTIFY_REF_POINT_EZP(dim2_len, dim2_ref, dim2_spc, overlap)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim2_len, &
& overlap, &
& proc_id, &
& proc_count, &
& ierror
DOUBLE PRECISION :: dim2_ref, &
& dim2_spc

CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, proc_count, ierror)

! Check for errors in user input.
CALL ERROR_HANDLING(proc_id)

! Run in serial if there is only one processor.
IF (proc_count .EQ. 1) THEN
  CALL SERIAL
ELSE
  CALL PARALLEL(dim2_len, dim2_ref, dim2_spc, overlap, proc_id, proc_count)
END IF


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SERIAL. Handles the execution if there is only one processor.
!
! STRUCTURE: N/A If there is only one processor, nothing needs to
! be done.
!
! VARIABLES: N/A
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE SERIAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

END SUBROUTINE SERIAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PARALLEL. Handles the execution if there is more than one processor.
!
! STRUCTURE:
!
! VARIABLES: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - proc_id: Processor ID (INTEGER).
! - proc_count: Total number of processors (INTEGER).
! - dim2_len_raw: Dim2 length of the interior of the sub-grid (INTEGER).
! - dim2_len_total: Dim2 length of the grid (INTEGER).
! - dim2_len_pre: Dim2 length of all processors with lesser ID, i.e., the global
! dim2 index for the non-overlap part of the sub-grid (INTEGER).
! - ierror: Integer for holding MPI error flag IDs (INTEGER).
! - i: Counting index used in DO loops (INTEGER).
! - dim2_len_list: Array of the dim2_len for all processors (INTEGER,
! DIMENSION(proc_count)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE PARALLEL(dim2_len, dim2_ref, dim2_spc, overlap, proc_id, proc_count)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim2_len, &
& overlap, &
& dim2_len_raw, &
& dim2_len_total, &
& dim2_len_pre, &
& proc_id, &
& proc_count, &
& ierror, &
& i
INTEGER, DIMENSION(proc_count) :: dim2_len_list
DOUBLE PRECISION :: dim2_ref, &
& dim2_spc

! Calculate dim2 of sub-grid to send along.
IF (MOD(proc_id, proc_count-1) .EQ. 0)  THEN
  ! First or last processor processor.
  dim2_len_raw = dim2_len - overlap
ELSE
  ! Middle processor.
  dim2_len_raw = dim2_len - 2 * overlap
END IF

CALL MPI_ALLREDUCE(dim2_len_raw, dim2_len_total, 1, MPI_INTEGER, MPI_SUM, &
& MPI_COMM_WORLD, ierror)

! Calculate dim2_len for each sub-grid, to avoid communication between processors.
  ! Divide dim2_len_total as evenly as possible.
dim2_len_raw = dim2_len_total/proc_count
  ! Now check if sub-grid of a given processor needs to take additional dim2_len,
  ! i.e., dim2_len_total % proc_count != 0.
DO i = 1, proc_count
  IF ((i-1) .LT. MODULO(dim2_len_total, proc_count)) THEN
    dim2_len_list(i) = dim2_len_raw + 1
    ! The indexing is weird here since processors start at index 0 but arrays
    ! start at index 1. Thus, entry ii of this corresponds to processor (i-1).
  ELSE
    dim2_len_list(i) = dim2_len_raw
  END IF
END DO

! Calculate the lesser dim2 grid indices of the sub-grid (the global coordinate
! of the sub-grid). All processors besides 0 need 1 more and to subtract the
! overlap to get the correct index.
dim2_len_pre = SUM(dim2_len_list(1:proc_id))
IF (proc_id .NE. 0) THEN
  dim2_len_pre = dim2_len_pre - overlap
END IF

dim2_ref = dim2_ref + dim2_len_pre * dim2_spc

END SUBROUTINE PARALLEL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ERROR_HANDLING. Checks for input error and other possible complications.
!
! ERRORS TO CHECK: N/A. All possible errors checked by subroutines before this.
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

END SUBROUTINE IDENTIFY_REF_POINT_EZP
