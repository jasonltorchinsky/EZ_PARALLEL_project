!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SHARE_SUBGRID_BOUNDARIES_DBLE_REAL. Communicates the sub-grid boundary to
! neighboring sub-grids, for DOUBLE PRECISION sub-grids.
!
! ARGUMENTS: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - sub_grid: The sub-grid of the processor, or the grid if only one processor
! is in use (DOUBLE PRECISION, DIMENSION(dim1_len, dim2_len)).
!
! Written By: Jason Turner
! Last Updated: January 15, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DBLE_REAL_EZP(dim1_len, dim2_len, overlap, &
& sub_grid)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& overlap, &
& proc_id, &
& proc_count, &
& ierror
DOUBLE PRECISION, DIMENSION(dim1_len, dim2_len) :: sub_grid

CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, proc_count, ierror)

! Check for errors in user input.
CALL ERROR_HANDLING(proc_id)

! Run in serial if there is only one processor.
IF (proc_count .EQ. 1) THEN
  CALL SERIAL
ELSE
  CALL PARALLEL(dim1_len, dim2_len, overlap, sub_grid, proc_id, proc_count)
END IF

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SERIAL. Handles the execution if there is only one processor.
!
! STRUCTURE: N/A. One processor implies there is no sub-grid boundary to
! communicate.
!
! VARIABLES: - variable1: variable1 descripotion.
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
! STRUCTURE: 1) Set the MPI message tages and zero out the SEND/RECV boundary
! arrays.
! 2) SEND/RECV boundary based on the following cases: a) Interior subgrid,
! needs to SEND/RECV both of its boundary arrays, b) first processor only needs
! to SEND/RECV next boundary, and c) last processor only needs to SEND/RECV prev
! boundary. Even processors SEND then RECV, odd processors RECV then SEND to
! reduce buffering.
!
! VARIABLES: - dim1_len, dim2_len: The shape of the sub-sub_grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - sub_grid: The sub-grid of the processor, or the grid if only one processor
! is in use (DOUBLE PRECISION, DIMENSION(dim1_len, dim2_len)).
! - proc_id: Processor ID (INTEGER).
! - proc_count: Total number of processors (INTEGER).
! - prev_send_tag, prev_recv_tag: MPI message tags for SEND to/RECV from the
! previous processor ID (INTEGER).
! - next_send_tag, next_recv_tag: MPI message tags for SEND to/RECV from the
! next processor ID (INTEGER).
! - prev_send_boundary, prev_recv_boundary: The arrays for storing the values
! to SEND to/RECV from the previous processor ID (DOUBLE PRECISION,
! DIMENSION(dim1_len, overlap)).
! - next_send_boundary, next_recv_boundary: The arrays for storing the values
! to SEND to/RECV from the next processor ID (DOUBLE PRECISION,
! DIMENSION(dim1_len, overlap)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE PARALLEL(dim1_len, dim2_len, overlap, sub_grid, proc_id, proc_count)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& overlap, &
& proc_id, &
& proc_count, &
& prev_send_tag, &
& next_recv_tag, &
& next_send_tag, &
& prev_recv_tag, &
& ierror, &
& status(MPI_STATUS_SIZE)
DOUBLE PRECISION, DIMENSION(dim1_len, dim2_len) :: sub_grid
DOUBLE PRECISION, DIMENSION(dim1_len, overlap) :: next_send_boundary, &
& prev_send_boundary, &
& next_recv_boundary, &
& prev_recv_boundary

prev_send_tag = 12
next_recv_tag = prev_send_tag
next_send_tag = 93
prev_recv_tag = next_send_tag

next_send_boundary = 0.0
prev_send_boundary = 0.0
next_recv_boundary = 0.0
prev_recv_boundary = 0.0

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Communicate boundary points of my subgrid to neighbors. 3 Cases:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 1. Interior subgrid, needs to send both of its boundary vectors.
IF (MODULO(proc_id, proc_count-1) .NE. 0) THEN
  next_send_boundary = sub_grid(:,dim2_len-2*overlap+1:dim2_len-overlap)
  prev_send_boundary = sub_grid(:,overlap+1:2*overlap)

  ! For buffering, even processors SEND/RECV and odd ones RECV/SEND.
  IF (MODULO(proc_id, 2) .EQ. 0) THEN
    CALL MPI_SEND(next_send_boundary, dim1_len*overlap, MPI_DOUBLE_PRECISION, &
    &  proc_id+1, next_send_tag, MPI_COMM_WORLD, ierror)
    CALL MPI_SEND(prev_send_boundary, dim1_len*overlap, MPI_DOUBLE_PRECISION, &
    &  proc_id-1, prev_send_tag, MPI_COMM_WORLD, ierror)

    CALL MPI_RECV(next_recv_boundary, dim1_len*overlap, MPI_DOUBLE_PRECISION, &
    &  proc_id+1, next_recv_tag, MPI_COMM_WORLD, status, ierror)
    CALL MPI_RECV(prev_recv_boundary, dim1_len*overlap, MPI_DOUBLE_PRECISION, &
    &  proc_id-1, prev_recv_tag,  MPI_COMM_WORLD, status, ierror)

  ELSE IF (MODULO(proc_id, 2) .EQ. 1) THEN
    CALL MPI_RECV(next_recv_boundary, dim1_len*overlap, MPI_DOUBLE_PRECISION, &
    &  proc_id+1, next_recv_tag, MPI_COMM_WORLD, status, ierror)
    CALL MPI_RECV(prev_recv_boundary, dim1_len*overlap, MPI_DOUBLE_PRECISION, &
    &  proc_id-1, prev_recv_tag,  MPI_COMM_WORLD, status, ierror)

    CALL MPI_SEND(next_send_boundary, dim1_len*overlap, MPI_DOUBLE_PRECISION, &
    &  proc_id+1, next_send_tag, MPI_COMM_WORLD, ierror)
    CALL MPI_SEND(prev_send_boundary, dim1_len*overlap, MPI_DOUBLE_PRECISION, &
    &  proc_id-1, prev_send_tag, MPI_COMM_WORLD, ierror)
  END IF

  sub_grid(:,dim2_len-overlap+1:dim2_len) = next_recv_boundary
  sub_grid(:,1:overlap) = prev_recv_boundary

! 2. Processor 0 only needs to send one boundary, SEND/RECV.
ELSE IF (proc_id .EQ. 0) THEN
  next_send_boundary = sub_grid(:,dim2_len-2*overlap+1:dim2_len-overlap)

  CALL MPI_SEND(next_send_boundary, dim1_len*overlap, MPI_DOUBLE_PRECISION, &
  &  proc_id+1, next_send_tag, MPI_COMM_WORLD, ierror)

  CALL MPI_RECV(next_recv_boundary, dim1_len*overlap, MPI_DOUBLE_PRECISION, &
  &  proc_id+1, next_recv_tag, MPI_COMM_WORLD, status, ierror)

  sub_grid(:,dim2_len-overlap+1:dim2_len) = next_recv_boundary

! 3. Last processor only needs to send one boundary, but follows even SEND/RECV
! odd RECV/SEND.
ELSE IF (proc_id .EQ. proc_count-1) THEN
  prev_send_boundary = sub_grid(:,overlap+1:2*overlap)

  ! For buffering, even processors send/recv and odd ones recv/send.
  IF (MODULO(proc_id, 2) .EQ. 0) THEN
    CALL MPI_SEND(prev_send_boundary, dim1_len*overlap, MPI_DOUBLE_PRECISION, &
    & proc_id-1, prev_send_tag, MPI_COMM_WORLD, ierror)

    CALL MPI_RECV(prev_recv_boundary, dim1_len*overlap, MPI_DOUBLE_PRECISION, &
    & proc_id-1, prev_recv_tag,  MPI_COMM_WORLD, status, ierror)

  ELSE IF (MODULO(proc_id, 2) .EQ. 1) THEN
    CALL MPI_RECV(prev_recv_boundary, dim1_len*overlap, MPI_DOUBLE_PRECISION, &
    & proc_id-1, prev_recv_tag,  MPI_COMM_WORLD, status, ierror)

    CALL MPI_SEND(prev_send_boundary, dim1_len*overlap, MPI_DOUBLE_PRECISION, &
    & proc_id-1, prev_send_tag, MPI_COMM_WORLD, ierror)
  END IF

  sub_grid(:,1:overlap) = prev_recv_boundary

! 4. ERROR STOP.
ELSE
  PRINT *, 'EZ_PARALLEL: Issue with processor ', proc_id, ' in subroutine ', &
  & 'call SHARE_SUBGRID_BOUNDARIES_DBLE_REAL (which is also called in ', &
  & 'SHARE_SUBGRID_BOUNDARIES_DBLE_CMPLX). You should ', &
  & 'never actually see this error, it means that the aforementioned ', &
  & 'processor does not have an ID between 0 and ', proc_count, &
  & '. I am not sure what you did to get here nor how I messed up, but look ', &
  & 'in the share_subgrid_boundaries.f file in the EZ_PARALLEL ', &
  & 'directory, the IF statement there messed up. I am going to force an ', &
  & 'ERROR STOP now.'
  ERROR STOP 'EZ_PARALLEL: Good luck fixing this, and let me know how you do!'
END IF

END SUBROUTINE PARALLEL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ERROR_HANDLING. Checks for input error and other possible complications.
!
! ERRORS TO CHECK: - N/A. (Assume that users will not enter negative matrix
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

END SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DBLE_REAL_EZP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
