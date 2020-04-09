!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The EZ_PARALLEL sub-grid boundary communication subroutine for DOUBLE
!! PRECISION sub-grids.
!> Communicates the sub-grid boundary to neighboring sub-grids, for DOUBLE
!! PRECISION sub-grids.
!
!> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
!! number of "unique" grid decompositions.
!> \param[in] sub_grid The sub-grid belonging to the processor. This is stored
!! under the original variable for the grid, e.g., if the serial code uses
!! heat_grid as the grid, then the heat_grid should be passed to this subroutine
!! as the sub_grid argument.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DBLE(decomp_id, sub_grid)

  USE MPI
  USE EZ_PARALLEL_STRUCTS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: decomp_id
  DOUBLE PRECISION, INTENT(INOUT) :: sub_grid(:,:)

  ! Check for errors in user input.
  CALL ERROR_HANDLING_DBLE(decomp_id)

  ! Run in serial if there is only one processor.
  IF (proc_count .EQ. 1) THEN
     CALL SERIAL_DBLE
  ELSE
     CALL PARALLEL_DBLE(decomp_id, sub_grid)
  END IF

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The serial EZ_PARALLEL sub-grid boundary communication subroutine for
  !! DOUBLE PRECISION sub-grids. 
  !> When running in serial, no sub-grid boundary communication occurs. This
  !! subroutine is included for file consistency, and may be depreciated in the
  !! future.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE SERIAL_DBLE
    
    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

  END SUBROUTINE SERIAL_DBLE


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The parallel EZ_PARALLEL sub-grid boundary communication subroutine for
  !! DOUBLE PRECISION sub-grids.
  !> Communicates sub-grid boundary to neighboring sub-grids. Note, since
  !! we are not using periodic boundary conditions, the first and last
  !! sub-grids only have one neighbor.
  !
  !> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
  !! number of "unique" grid decompositions.
  !> \param[in] sub_grid The sub-grid belonging to the processor. This is stored
  !! under the original variable for the grid, e.g., if the serial code uses
  !! heat_grid as the grid, then the heat_grid should be passed to this
  !! subroutine as the sub_grid argument.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE PARALLEL_DBLE(decomp_id, sub_grid)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: decomp_id
    DOUBLE PRECISION, INTENT(INOUT) :: sub_grid(:,:)
    INTEGER :: ierror
    INTEGER :: status(MPI_STATUS_SIZE)

    ! First sub-grid only SEND/RECV to/from one neighbor.
    IF (proc_id .EQ. 0) THEN
       CALL MPI_SEND(sub_grid, 1, grid_decomps(decomp_id)%SEND_BOUNDARIES(2), &
            proc_id+1, proc_id+1, MPI_COMM_WORLD, ierror)
       
       CALL MPI_RECV(sub_grid, 1, grid_decomps(decomp_id)%RECV_BOUNDARIES(2), &
            proc_id+1, proc_id, MPI_COMM_WORLD, status, ierror)
       
    ! Last sub-grid only SEND/RECV to/from one neighbor.
    ELSE IF (proc_id .EQ. proc_count-1) THEN
       ! If proc_id even, SEND then RECV.
       IF (MOD(proc_id,2) .EQ. 0) THEN
          CALL MPI_SEND(sub_grid, 1, grid_decomps(decomp_id)%SEND_BOUNDARIES(1), &
               proc_id-1, proc_id-1, MPI_COMM_WORLD, ierror)
          
          CALL MPI_RECV(sub_grid, 1, grid_decomps(decomp_id)%RECV_BOUNDARIES(1), &
               proc_id-1, proc_id, MPI_COMM_WORLD, status, ierror)
          
       ! If proc_id odd, RECV then SEND.
       ELSE
          CALL MPI_RECV(sub_grid, 1, grid_decomps(decomp_id)%RECV_BOUNDARIES(1), &
               proc_id-1, proc_id, MPI_COMM_WORLD, status, ierror)
          
          CALL MPI_SEND(sub_grid, 1, grid_decomps(decomp_id)%SEND_BOUNDARIES(1), &
               proc_id-1, proc_id-1, MPI_COMM_WORLD, ierror)
       END IF
       
    ! Interior sub-grids SEND/RECV to/from two neighbors.
    ELSE
       ! If proc_id even, SEND then RECV.
       IF (MOD(proc_id,2) .EQ. 0) THEN
          CALL MPI_SEND(sub_grid, 1, grid_decomps(decomp_id)%SEND_BOUNDARIES(1), &
               proc_id-1, proc_id-1, MPI_COMM_WORLD, ierror)
          CALL MPI_SEND(sub_grid, 1, grid_decomps(decomp_id)%SEND_BOUNDARIES(2), &
               proc_id+1, proc_id+1, MPI_COMM_WORLD, ierror)
          
          CALL MPI_RECV(sub_grid, 1, grid_decomps(decomp_id)%RECV_BOUNDARIES(1), &
               proc_id-1, proc_id, MPI_COMM_WORLD, status, ierror)
          CALL MPI_RECV(sub_grid, 1, grid_decomps(decomp_id)%RECV_BOUNDARIES(2), &
               proc_id-1, proc_id, MPI_COMM_WORLD, status, ierror)
          
       ! If proc_id odd, RECV then SEND
       ELSE
          CALL MPI_RECV(sub_grid, 1, grid_decomps(decomp_id)%RECV_BOUNDARIES(1), &
               proc_id-1, proc_id, MPI_COMM_WORLD, status, ierror)
          CALL MPI_RECV(sub_grid, 1, grid_decomps(decomp_id)%RECV_BOUNDARIES(2), &
               proc_id-1, proc_id, MPI_COMM_WORLD, status, ierror)
          
          CALL MPI_SEND(sub_grid, 1, grid_decomps(decomp_id)%SEND_BOUNDARIES(1), &
               proc_id-1, proc_id-1, MPI_COMM_WORLD, ierror)
          CALL MPI_SEND(sub_grid, 1, grid_decomps(decomp_id)%SEND_BOUNDARIES(2), &
               proc_id+1, proc_id+1, MPI_COMM_WORLD, ierror)
       END IF
    END IF

  END SUBROUTINE PARALLEL_DBLE

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The EZ_PARALLEL sub-grid boundary communication error handling subroutine
  !! for DOUBLE PRECISION sub-grids.
  !> Aborts the program if the grid decomposition parameters trigger any of the
  !!following errors:
  !! <ul>
  !! <li> decomp_id...
  !!    <ol>
  !!    <li> is outside of (1, ..., decomp_count)
  !!    </ol>
  !! </ul>
  !
  !> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
  !! number of "unique" grid decompositions.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE ERROR_HANDLING_DBLE(decomp_id)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: decomp_id
    INTEGER :: error_code
    INTEGER :: ierror
    LOGICAL :: error_flag !< Error flag to stop all processors.

    error_flag = .FALSE.


    ! Check if decomp_id outside of (1, ..., decomp_count)
    IF ((decomp_id .LT. 1) .OR. (decomp_id .GT. decomp_count)) THEN
       IF (proc_id .EQ. 0) THEN
          PRINT *, "ERROR: Attempted to communicate sub-grid boundaries of ", &
               " sub-grids of grid with decomposition ID ", decomp_ID, &
               " which is outside of 1, ..., decomp_count. decomp_count: ", &
               decomp_count, "." 
       END IF
       error_flag = .TRUE.
    END IF

    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF

  END SUBROUTINE ERROR_HANDLING_DBLE

END SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DBLE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The EZ_PARALLEL sub-grid boundary communication subroutine for DOUBLE COMLPEX
!! sub-grids.
!> Communicates the sub-grid boundary to neighboring sub-grids, for DOUBLE
!! COMPLEX sub-grids.
!
!> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
!! number of "unique" grid decompositions.
!> \param[in] sub_grid The sub-grid belonging to the processor. This is stored
!! under the original variable for the grid, e.g., if the serial code uses
!! heat_grid as the grid, then the heat_grid should be passed to this subroutine
!! as the sub_grid argument.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DCMPLX(decomp_id, sub_grid)

  USE MPI
  USE EZ_PARALLEL_STRUCTS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: decomp_id
  DOUBLE COMPLEX, INTENT(INOUT) :: sub_grid(:,:)

  ! Check for errors in user input.
  CALL ERROR_HANDLING_DCMPLX(decomp_id)

  ! Run in serial if there is only one processor.
  IF (proc_count .EQ. 1) THEN
     CALL SERIAL_DCMPLX
  ELSE
     CALL PARALLEL_DCMPLX(decomp_id, sub_grid)
  END IF

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The serial EZ_PARALLEL sub-grid boundary communication subroutine for
  !! DOUBLE COMPLEX sub-grids. 
  !> When running in serial, no sub-grid boundary communication occurs. This
  !! subroutine is included for file consistency, and may be depreciated in the
  !! future.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE SERIAL_DCMPLX
    
    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

  END SUBROUTINE SERIAL_DCMPLX


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The parallel EZ_PARALLEL sub-grid boundary communication subroutine for
  !! DOUBLE COMPLEX sub-grids.
  !> Communicates sub-grid boundary to neighboring sub-grids. Note, since
  !! we are not using periodic boundary conditions, the first and last
  !! sub-grids only have one neighbor.
  !
  !> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
  !! number of "unique" grid decompositions.
  !> \param[in] sub_grid The sub-grid belonging to the processor. This is stored
  !! under the original variable for the grid, e.g., if the serial code uses
  !! heat_grid as the grid, then the heat_grid should be passed to this
  !! subroutine as the sub_grid argument.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE PARALLEL_DCMPLX(decomp_id, sub_grid)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: decomp_id
    DOUBLE COMPLEX, INTENT(INOUT) :: sub_grid(:,:)
    INTEGER :: ierror
    INTEGER :: status(MPI_STATUS_SIZE)

    ! First sub-grid only SEND/RECV to/from one neighbor.
    IF (proc_id .EQ. 0) THEN
       CALL MPI_SEND(sub_grid, 1, grid_decomps(decomp_id)%SEND_BOUNDARIES(2), &
            proc_id+1, proc_id+1, MPI_COMM_WORLD, ierror)
       CALL MPI_RECV(sub_grid, 1, grid_decomps(decomp_id)%RECV_BOUNDARIES(2), &
            proc_id+1, proc_id, MPI_COMM_WORLD, ierror)
       
    ! Last sub-grid only SEND/RECV to/from one neighbor.
    ELSE IF (proc_id .EQ. proc_count-1) THEN
       ! If proc_id even, SEND then RECV.
       IF (MOD(proc_id,2) .EQ. 0) THEN
          CALL MPI_SEND(sub_grid, 1, grid_decomps(decomp_id)%SEND_BOUNDARIES(1), &
               proc_id-1, proc_id-1, MPI_COMM_WORLD, ierror)
          
          CALL MPI_RECV(sub_grid, 1, grid_decomps(decomp_id)%RECV_BOUNDARIES(1), &
               proc_id-1, proc_id, MPI_COMM_WORLD, ierror)
          
       ! If proc_id odd, RECV then SEND.
       ELSE
          CALL MPI_RECV(sub_grid, 1, grid_decomps(decomp_id)%RECV_BOUNDARIES(1), &
               proc_id-1, proc_id, MPI_COMM_WORLD, ierror)
          
          CALL MPI_SEND(sub_grid, 1, grid_decomps(decomp_id)%SEND_BOUNDARIES(1), &
               proc_id-1, proc_id-1, MPI_COMM_WORLD, ierror)
       END IF
       
    ! Interior sub-grids SEND/RECV to/from two neighbors.
    ELSE
       ! If proc_id even, SEND then RECV.
       IF (MOD(proc_id,2) .EQ. 0) THEN
          CALL MPI_SEND(sub_grid, 1, grid_decomps(decomp_id)%SEND_BOUNDARIES(1), &
               proc_id-1, proc_id-1, MPI_COMM_WORLD, ierror)
          CALL MPI_SEND(sub_grid, 1, grid_decomps(decomp_id)%SEND_BOUNDARIES(2), &
               proc_id+1, proc_id+1, MPI_COMM_WORLD, ierror)
          
          CALL MPI_RECV(sub_grid, 1, grid_decomps(decomp_id)%RECV_BOUNDARIES(1), &
               proc_id-1, proc_id, MPI_COMM_WORLD, ierror)
          CALL MPI_RECV(sub_grid, 1, grid_decomps(decomp_id)%RECV_BOUNDARIES(2), &
               proc_id-1, proc_id, MPI_COMM_WORLD, ierror)
          
       ! If proc_id odd, RECV then SEND
       ELSE
          CALL MPI_RECV(sub_grid, 1, grid_decomps(decomp_id)%RECV_BOUNDARIES(1), &
               proc_id-1, proc_id, MPI_COMM_WORLD, ierror)
          CALL MPI_RECV(sub_grid, 1, grid_decomps(decomp_id)%RECV_BOUNDARIES(2), &
               proc_id-1, proc_id, MPI_COMM_WORLD, ierror)
          
          CALL MPI_SEND(sub_grid, 1, grid_decomps(decomp_id)%SEND_BOUNDARIES(1), &
               proc_id-1, proc_id-1, MPI_COMM_WORLD, ierror)
          CALL MPI_SEND(sub_grid, 1, grid_decomps(decomp_id)%SEND_BOUNDARIES(2), &
               proc_id+1, proc_id+1, MPI_COMM_WORLD, ierror)
       END IF
    END IF

  END SUBROUTINE PARALLEL_DCMPLX

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The EZ_PARALLEL sub-grid boundary communication error handling subroutine
  !! for DOUBLE COMPLEX sub-grids.
  !> Aborts the program if the grid decomposition parameters trigger any of the
  !!following errors:
  !! <ul>
  !! <li> decomp_id...
  !!    <ol>
  !!    <li> is outside of (1, ..., decomp_count)
  !!    </ol>
  !! </ul>
  !
  !> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
  !! number of "unique" grid decompositions.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE ERROR_HANDLING_DCMPLX(decomp_id)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: decomp_id
    INTEGER :: error_code
    INTEGER :: ierror
    LOGICAL :: error_flag !< Error flag to stop all processors.

    error_flag = .FALSE.


    ! Check if decomp_id outside of (1, ..., decomp_count)
    IF ((decomp_id .LT. 1) .OR. (decomp_id .GT. decomp_count)) THEN
       IF (proc_id .EQ. 0) THEN
          PRINT *, "ERROR: Attempted to communicate sub-grid boundaries of ", &
               " sub-grids of grid with decomposition ID ", decomp_ID, &
               " which is outside of 1, ..., decomp_count. decomp_count: ", &
               decomp_count, "." 
       END IF
       error_flag = .TRUE.
    END IF

    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF

  END SUBROUTINE ERROR_HANDLING_DCMPLX

END SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DCMPLX

