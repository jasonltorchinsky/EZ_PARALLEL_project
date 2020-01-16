!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! DECOMP_GRID. Defines a grid decomposition by dividing the grid as evenly as
! possible along the second dimension (dim2), with the desired amount of space
! needed to store the overlap regions (the regions that need to be communicated
! between sub-grids). For example, this will divide an (8,9) grid among 3
! processors with an overlap of 1 into three sub-grids of size (8,4), (8,5),
! (8,4).
!
! NOTES: - overlap is the number of rows a sub-grid needs to borrow from one of
! its neighbors. For example, a first-order center difference approximation
! of a second derivative needs the points i-1, i, i+1.
!
! ARGUMENTS: - dim2_len_total: The size of dim2 of the grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE DECOMP_GRID_EZP(dim2_len_total, overlap)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim2_len_total, &
& overlap, &
& proc_id, &
& proc_count, &
& ierror

CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, proc_count, ierror)

! Check for errors in user input.
CALL ERROR_HANDLING(dim2_len_total, overlap, proc_id, proc_count)

! Run in serial if there is only one processor.
IF (proc_count .EQ. 1) THEN
  CALL SERIAL
ELSE
  CALL PARALLEL(dim2_len_total, overlap, proc_id, proc_count)
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

CONTINUE

END SUBROUTINE SERIAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PARALLEL. Handles the execution if there is more than one processor.
!
! STRUCTURE: 1) Divide dim2_len as easily as possible, adding 1 where necessary,
! e.g., splitting dim2_len = 9 among two processors, processor 0 needs 5 and
! processor 1 needs 4.
!
! VARIABLES: - dim2_len_total: The size of dim2 of the grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - proc_id: Processor ID (INTEGER).
! - proc_count: Total number of processors (INTEGER).
! - dim2_len: Size of dim2 of the sub-grid (INTEGER).
! - ierror: Integer for holding MPI error flag IDs (INTEGER).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE PARALLEL(dim2_len_total, overlap, proc_id, proc_count)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim2_len_total, &
& overlap, &
& proc_id, &
& proc_count, &
& dim2_len, &
& ierror

! Divide dim2_len_total as evenly as possible, ignoring overlap.
dim2_len = dim2_len_total/proc_count

! Add 1 to dim2_len if the sub-grid would be big enough, see STRUCTURE for
! example.
IF (proc_id .LT. MODULO(dim2_len_total, proc_count)) THEN
  dim2_len = dim2_len + 1
END IF

! Increase dim2_len for overlap.
! Interior sub-grid.
IF (MOD(proc_id, proc_count-1) .NE. 0) THEN
  dim2_len = dim2_len + 2 * overlap
! Boundary sub-grid.
ELSE IF (MODULO(proc_id, proc_count-1) .EQ. 0) THEN
  dim2_len = dim2_len + overlap
END IF

dim2_len_total = dim2_len

END SUBROUTINE PARALLEL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ERROR_HANDLING. Checks for input error and other possible complications.
!
! ERRORS TO CHECK: - overlap too large (overlap >= dim2_len for any processor).
! - overlap < 0.
!
! VARIABLES: - dim2_len_total: The size of dim2 of the grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - proc_id: Processor ID (INTEGER).
! - proc_count: Total number of processors (INTEGER).
! - proc_err_flag - Flag is triggered if processor encounters an error
! (LOGICAL).
! - global_err_flag - Flag is triggered if any processor encouters an error
! (LOGICAL).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE ERROR_HANDLING(dim2_len_total, overlap, proc_id, proc_count)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim2_len_total, &
& overlap, &
& proc_id, &
& proc_count
LOGICAL :: proc_err_flag, &
& global_err_flag

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

proc_err_flag = .FALSE.
global_err_flag = .FALSE.

! overlap too large. Only processor 0 needs to check this.
IF ((proc_id .EQ. 0) .AND. (overlap .GE. dim2_len_total/proc_count)) THEN
  proc_err_flag = .TRUE.
  PRINT *, 'EZ_PARALLEL: Issue with subroutine call DECOMP_GRID. overlap: ', &
  & overlap, ' is greater than or equal to dim2_len ', &
  & dim2_len_total/proc_count, ' of some processor.'
  PRINT *, 'Please decrease the overlap, decrease the number of processors, ', &
  & 'or increase the size of the second dimension of the grid.'
END IF

! overlap < 1. Only processor 0 needs to check this.
IF ((proc_id .EQ. 0) .AND. (overlap .LT. 0)) THEN
  proc_err_flag = .TRUE.
  PRINT *, 'EZ_PARALLEL: Issue with subroutine call DECOMP_GRID. overlap: ', &
  & overlap, ' is negative.'
  PRINT *, 'Please increase the overlap.'
END IF

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

END SUBROUTINE DECOMP_GRID_EZP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
