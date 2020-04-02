!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The EZ_PARALLEL grid decomposition subroutine.
!> Given the column count of the grid, returns the column count of the sub-grid,
!! including overlap (for use of the vertical slab decomposition). Also stores
!! the data in the grid_decomps list for use in other subroutines.
!
!> \param[in] row_count Number of rows in the grid.
!> \param[inout] col_count Number of columns in the grid, changes to number of
!! columns in the sub-grid including overlap, so that the user's time-stepping
!! code may remain unchanged.
!> \param[in] overlap Number of extra columns needed by each sub-grid to
!! successfully step forward in time.
!> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
!! number of "unique" grid decompositions.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE DECOMP_GRID_EZP(row_count, col_count, overlap, decomp_id)

  USE MPI

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: row_count, &
       overlap, &
       decomp_id
  INTEGER, INTENT(INOUT) :: col_count

  ! Check for errors in user input.
  CALL ERROR_HANDLING(row_count, col_count, overlap, decomp_id)

  IF (proc_count .EQ. 1) THEN
     ! Run in serial is there is only one processor.
     CALL SERIAL
  ELSE
     ! Run in parallel if there is more than one processor.
     CALL PARALLEL(row_count, col_count, overlap, decomp_id)
  END IF

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The serial EZ_PARALLEL grid decomposition subroutine.
  !> When running in serial, no grid decomposition occurs. This subroutine is
  !! included for file consistency, and may be depreciated in the future.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE SERIAL

    USE MPI

    IMPLICIT NONE

    CONTINUE

  END SUBROUTINE SERIAL

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The parallel EZ_PARALLEL grid decomposition subroutine.
  !> Given the column count of the grid, returns the column count of the
  !! sub-grid, including overlap (for use of the vertical slab decomposition).
  !! Also stores the data in the grid_decomps list for use in other subroutines.
  !
  !> \param[in] row_count Number of rows in the grid.
  !> \param[inout] col_count Number of columns in the grid, changes to number of
  !! columns in the sub-grid including overlap, so that the user's time-stepping
  !! code may remain unchanged.
  !> \param[in] overlap Number of extra columns needed by each sub-grid to
  !! successfully step forward in time.
  !> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
  !! number of "unique" grid decompositions.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE PARALLEL(row_count, col_count, overlap, decomp_id)

    USE MPI

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: row_count, &
         overlap, &
         decomp_id
    INTEGER, INTENT(INOUT) :: col_count
    INTEGER :: row_count_s, & !< Row count for the horizontal slab decomposition
         !! for various processors.
         col_count_s, & !< Column count for the vertical slab decomposition for
         !! various processors.
         i !< Iterator variable for DO loops.
    TYPE(GRID_DECOMPOSITION) :: grid_decomp !< The grid decomposition.

    ! Store grid data in the grid.
    grid_decomp%row_count_g = row_count
    grid_decomp%col_count_g = col_count
    grid_decomp%overlap = overlap
    grid_decomp%decomp_id = decomp_id
    
    ! Generate the horzontal-slab decomposition.
    ALLOCATE(grid_decomp%row_decomp(proc_count))
    ! Divide row_count as evenly as possible among all the processors, ignoring
    ! overlap.
    row_count_s = row_count / proc_count
    DO i = 1, proc_count
       IF (i .LE. MODULO(row_count, proc_count)) THEN
          grid_decomp%row_decomp(i) = row_count_s + 1
       ELSE
          grid_decomp%row_decomp(i) = row_count_s
       END IF
    END DO

    ! Generate the vertical-slab decomposition without overlap.
    ALLOCATE(grid_decomp%col_decomp(proc_count))
    ! Divide col_count as evenly as possible among all the processors, ignoring
    ! overlap.
    col_count_s = col_count/proc_count
    DO i = 1, proc_count
       IF (i .LE. MODULO(col_count, proc_count)) THEN
          grid_decomp%col_decomp(i) = col_count_s + 1
       ELSE
          grid_decomp%col_decomp(i) = col_count_s
       END IF
    END DO

    ! Generate the vertical-slab decomposition with overlap.
    ALLOCATE(grid_decomp%col_decomp_ovlp(proc_count))
    ! Different from without overlap only if overlap .NE. 0
    IF (overlap .EQ. 0) THEN
       grid_decomp%col_decomp_ovlp = grid_decomp%col_decomp
    ELSE
       ! All sub-grids with one neighbor get overlap additional columns.
       grid_decomp%col_decomp_ovlp(1) = grid_decomp%col_decomp(1) + overlap
       grid_decomp%col_decomp_ovlp(proc_count) = &
            grid_decomp%col_decomp(proc_count) + overlap
       ! All sub-grids with two neighbors get 2*overlap additional columns.
       grid_decomp%col_decomp_ovlp(2:proc_count-1) = &
            grid_decomp%col_decomp(2:proc_count-1) + 2*overlap
    END IF

    ! Put grid_decomp into the list of unique grid_decomps
    grid_decomps(decomp_id) = grid_decomp

    ! Free allocated arrays
    DEALLOCATE(grid_decomp%row_decomp)
    DEALLOCATE(grid_decomp%col_decomp)
    DEALLOCATE(grid_decomp%col_decomp_ovlp)

  END SUBROUTINE PARALLEL

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The EZ_PARALLEL grid decomposition error handling subroutine.
  !> Aborts the program if the grid or grid decomposition parameters trigger any
  !! of the following errors:
  !! <ul>
  !! <li> row_count...
  !!    <ol>
  !!    <li> is not positive
  !!    </ol>
  !! <li> col_count...
  !!    <ol>
  !!    <li> is not positive
  !!    <li> is too small to be divided among the processors, ignoring the overlap
  !!    <li> is too small to be divided among the processors, including overlap
  !!    </ol>
  !! <li> overlap...
  !!    <ol>
  !!    <li> is negative
  !!    </ol>
  !! <li> decomp_id...
  !!    <ol>
  !!    <li> is outside of (1, ..., decomp_count)
  !!    </ol>
  !! </ul>
  !
  !> \param[in] row_count Number of rows in the grid.
  !> \param[in] col_count Number of columns in the grid, changes to number of
  !! columns in the sub-grid including overlap, so that the user's time-stepping
  !! code may remain unchanged.
  !> \param[in] overlap Number of extra columns needed by each sub-grid to
  !! successfully step forward in time.
  !> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
  !! number of "unique" grid decompositions.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE ERROR_HANDLING(row_count, col_count, overlap, decomp_id)

    USE MPI

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: row_count, &
         col_count, &
         overlap, &
         decomp_id
    INTEGER :: error_code, &
         ierror
    LOGICAL :: error_flag !< Error flag to stop all processors.

    error_flag = .FALSE.

    ! Check if row_count is not positive.
    IF (row_count .LT. 1) THEN
       IF (proc_id .EQ. 0) THEN
          PRINT *, "ERROR: Non-positive number of rows in grid with ", &
               " decomposition ID ", decomp_ID, ". Number of rows: ", &
               row_count, "." 
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if col_count is not positive.
    IF (col_count .LT. 1) THEN
       IF (proc_id .EQ. 0) THEN
          PRINT *, "ERROR: Non-positive number of columns in grid with ", &
               " decomposition ID ", decomp_ID, ". Number of columns: ", &
               col_count, "." 
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if col_count is too small to be divided among processors, excluding overlap.
    IF (col_count .LT. proc_count) THEN
       IF (proc_id .EQ. 0) THEN
          PRINT *, "ERROR: Number of columns in grid with decomposition ID ", &
               decomp_ID, " fewer than number of processors. Number of columns: ", &
               col_count, ". Number of processors: ", proc_count, "." 
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if col_count is too small to be divided among processors, including overlap.
    IF (col_count .LT. overlap*proc_count) THEN
       IF (proc_id .EQ. 0) THEN
          PRINT *, "ERROR: Number of columns too small to decompose grid ", &
               "without sub-grid boundaries overlapping for grid with ", &
               "decomposition ID ", decomp_ID, ". Number of columns: ", &
               col_count, ". Overlap: ", overlap , "." 
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if overlap is negative.
    IF (overlap .LT. 0) THEN
       IF (proc_id .EQ. 0) THEN
          PRINT *, "ERROR: Negative overlap in grid with decomposition ID ", &
               decomp_ID, ". overlap: ", overlap, "." 
       END IF
       error_flag = .TRUE.
    END IF

    ! Check if decomp_id outside of (1, ..., decomp_count)
    IF ((decomp_id .LT. 1) .OR. (decomp_id .GT. decomp_count)) THEN
       IF (proc_id .EQ. 0) THEN
          PRINT *, "ERROR: Decomposition ID in grid with decomposition ID ", &
               decomp_ID, " is outside of 1, ..., decomp_count. decomp_count: ", &
               decomp_count, "." 
       END IF
       error_flag = .TRUE.
    END IF

    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF

  END SUBROUTINE ERROR_HANDLING
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE DECOMP_GRID_EZP
