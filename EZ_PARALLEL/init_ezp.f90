!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The EZ_PARALLEL initialization subroutine.
!> Initializes MPI, allocates memory for the grid decomposition list, and stores
!! frequenctly used variables.
!
!> \param[in] decomp_count_L Number of unique grid decompositions,
!! "unique" being in terms of size and overlap (local to subroutine).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE INIT(decomp_count_L)
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: decomp_count_L
  INTEGER :: ierror
  
  CALL MPI_INIT(ierror)

  ! Check for errors in user input.
  CALL ERROR_HANDLING(decomp_count_L)
  
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, proc_count, ierror)
  
  decomp_count = decomp_count_L
  ALLOCATE(grid_decomps(decomp_count))

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The EZ_PARALLEL initialization error handling subroutine.
  !> Aborts the program if the grid decompositions parameters trigger any
  !! of the following errors:
  !! <ul>
  !! <li> decomp_count_L...
  !!    <ol>
  !!    <li> is not positive
  !!    </ol>
  !! </ul>
  !
  !> \param[in] decomp_count_L Number of unique grid decompositions (local to
  !! subroutine).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE ERROR_HANDLING(decomp_count_L)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: decomp_count_L
    INTEGER :: error_code, &
         ierror
    LOGICAL :: error_flag !< Error flag to stop all processors.

    error_flag = .FALSE.

    ! Check if decomp_count is not positive.
    IF (decomp_count_L .LT. 1) THEN
       IF (proc_id .EQ. 0) THEN
          PRINT *, "ERROR: Non-positive number of grid decompositions: ", &
               decomp_count_L, "." 
       END IF
       error_flag = .TRUE.
    END IF

    IF (error_flag) THEN
       CALL MPI_ABORT(MPI_COMM_WORLD, error_code, ierror)
    END IF

  END SUBROUTINE ERROR_HANDLING
 
END SUBROUTINE INIT


