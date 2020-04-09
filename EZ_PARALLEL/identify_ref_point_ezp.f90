!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The EZ_PARALLEL reference point identification subroutine for DOUBLE
!! PRECISION column spacing and column reference point.
!> Assuming a linear spacing of grid points, adjusts the physical position of
!! the reference point in the physical direction corresponding to the second
!! dimension of the grid.
!
!> \param[in] col_spc Physical spacing of columns in the grid, e.g., if the
!! second dimension of the grid represents the z-direction, then col_spc
!! corresponds to the spacing of grid points in the z-direction.
!> \param[inout] col_ref Physical position of the reference point in the
!! second dimension of the grid, e.g., the z-coordinate.
!> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
!! number of "unique" grid decompositions.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE IDENTIFY_REF_POINT_DBLE(col_spc, col_ref, decomp_id)

  USE MPI
  USE EZ_PARALLEL_STRUCTS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: decomp_id
  DOUBLE PRECISION, INTENT(IN) :: col_spc
  DOUBLE PRECISION, INTENT(INOUT) :: col_ref

  ! Check for errors in user input.
  CALL ERROR_HANDLING_DBLE(decomp_id)

  IF (proc_count .EQ. 1) THEN
     ! Run in serial is there is only one processor.
     CALL SERIAL_DBLE
  ELSE
     ! Run in parallel if there is more than one processor.
     CALL PARALLEL_DBLE(col_spc, col_ref, decomp_id)
  END IF

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The serial EZ_PARALLEL reference point identification subroutine for DOUBLE
  !! PRECISION column spacing and column reference point.
  !> When running in serial, no change in reference point occurs. This
  !! subroutine is included for file consistency, and may be depreciated in the
  !! future.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE SERIAL_DBLE

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    CONTINUE

  END SUBROUTINE SERIAL_DBLE

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The parallel EZ_PARALLEL reference point identification subroutine for DOUBLE
  !! PRECISION column spacing and column reference point.
  !> Assuming a linear spacing of grid points, calculates the physical position
  !! of the sub-grid reference point in the physical direction corresponding to
  !! the first dimension of the grid.
  !
  !> \param[in] col_spc Physical spacing of columns in the grid, e.g., if the
  !! second dimension of the grid represents the z-direction, then col_spc
  !! corresponds to the spacing of grid points in the z-direction.
  !> \param[inout] col_ref Physical position of the reference point in the
  !! second dimension of the grid, e.g., the z-coordinate.
  !> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
  !! number of "unique" grid decompositions.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE PARALLEL_DBLE(col_spc, col_ref, decomp_id)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: decomp_id
    DOUBLE PRECISION, INTENT(IN) :: col_spc
    DOUBLE PRECISION, INTENT(INOUT) :: col_ref

    IF (proc_id .EQ. 0) THEN
       col_ref = col_ref + col_spc * &
            (SUM(grid_decomps(decomp_id)%col_decomp(1:proc_id+1)) &
            - grid_decomps(decomp_id)%col_decomp(proc_id+1))
    ELSE
       col_ref = col_ref + col_spc * &
            (SUM(grid_decomps(decomp_id)%col_decomp(1:proc_id+1)) &
            - grid_decomps(decomp_id)%col_decomp(proc_id+1) &
            - grid_decomps(decomp_id)%overlap)
    END IF

  END SUBROUTINE PARALLEL_DBLE

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The EZ_PARALLEL reference point identification error handling subroutine
  !! for DOUBLE PRECISION column spacing and column reference point.
  !> Aborts the program if the grid or grid decomposition parameters trigger any
  !! of the following errors:
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
          PRINT *, "ERROR: Attempted to identify reference points of ", &
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

END SUBROUTINE IDENTIFY_REF_POINT_DBLE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The EZ_PARALLEL reference point identification subroutine for DOUBLE
!! COMPLEX column spacing and column reference point.
!> Assuming a linear spacing of grid points, adjusts the physical position of
!! the reference point in the physical direction corresponding to the second
!! dimension of the grid.
!
!> \param[in] col_spc Physical spacing of columns in the grid, e.g., if the
!! second dimension of the grid represents the z-direction, then col_spc
!! corresponds to the spacing of grid points in the z-direction.
!> \param[inout] col_ref Physical position of the reference point in the
!! second dimension of the grid, e.g., the z-coordinate.
!> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
!! number of "unique" grid decompositions.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE IDENTIFY_REF_POINT_DCMPLX(col_spc, col_ref, decomp_id)

  USE MPI
  USE EZ_PARALLEL_STRUCTS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: decomp_id
  DOUBLE COMPLEX, INTENT(IN) :: col_spc
  DOUBLE COMPLEX, INTENT(INOUT) :: col_ref

  ! Check for errors in user input.
  CALL ERROR_HANDLING_DCMPLX(decomp_id)

  IF (proc_count .EQ. 1) THEN
     ! Run in serial is there is only one processor.
     CALL SERIAL_DCMPLX
  ELSE
     ! Run in parallel if there is more than one processor.
     CALL PARALLEL_DCMPLX(col_spc, col_ref, decomp_id)
  END IF

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The serial EZ_PARALLEL reference point identification subroutine for DOUBLE
  !! COMPLEX column spacing and column reference point.
  !> When running in serial, no change in reference point occurs. This
  !! subroutine is included for file consistency, and may be depreciated in the
  !! future.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE SERIAL_DCMPLX

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    CONTINUE

  END SUBROUTINE SERIAL_DCMPLX

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The parallel EZ_PARALLEL reference point identification subroutine for DOUBLE
  !! COMPLEX column spacing and column reference point.
  !> Assuming a linear spacing of grid points, calculates the physical position
  !! of the sub-grid reference point in the physical direction corresponding to
  !! the first dimension of the grid.
  !
  !> \param[in] col_spc Physical spacing of columns in the grid, e.g., if the
  !! second dimension of the grid represents the z-direction, then col_spc
  !! corresponds to the spacing of grid points in the z-direction.
  !> \param[inout] col_ref Physical position of the reference point in the
  !! second dimension of the grid, e.g., the z-coordinate.
  !> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
  !! number of "unique" grid decompositions.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE PARALLEL_DCMPLX(col_spc, col_ref, decomp_id)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: decomp_id
    DOUBLE COMPLEX, INTENT(IN) :: col_spc
    DOUBLE COMPLEX, INTENT(INOUT) :: col_ref

    IF (proc_id .EQ. 0) THEN
       col_ref = col_ref + col_spc * &
            (SUM(grid_decomps(decomp_id)%col_decomp(1:proc_id+1)) &
            - grid_decomps(decomp_id)%col_decomp(proc_id+1))
    ELSE
       col_ref = col_ref + col_spc * &
            (SUM(grid_decomps(decomp_id)%col_decomp(1:proc_id+1)) &
            - grid_decomps(decomp_id)%col_decomp(proc_id+1) &
            - grid_decomps(decomp_id)%overlap)
    END IF

  END SUBROUTINE PARALLEL_DCMPLX

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The EZ_PARALLEL reference point identification error handling subroutine
  !! for DOUBLE COMPLEX column spacing and column reference point.
  !> Aborts the program if the grid or grid decomposition parameters trigger any
  !! of the following errors:
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
          PRINT *, "ERROR: Attempted to identify reference points of ", &
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

END SUBROUTINE IDENTIFY_REF_POINT_DCMPLX



