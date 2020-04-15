!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wiscons-Madison
!> @brief
!> The <tt>EZ_PARALLEL</tt> first-dimension wavenumber array generator
!! subroutine.
!> Generates an array containing the wavenumbers along the first dimension
!! corresponding to the local processor sub-grid in spectral space, for use in
!! finding the derivative along the first dimension in spectral space.
!> \remark We use the convention of wavenumbers going from 0, ...,
!! <tt>row_count</tt>/2, -<tt>row_count</tt>/2+1, ..., -1.
!
!> \param[in] order The order of derivative desired.
!> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
!! number of "unique" grid decompositions.
!> \param[inout] spec_dim1_deriv The array to hold the wavenumbers along the
!! first dimension.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE SPECTRAL_DIM1_DERIVATIVE(order, decomp_id, spec_dim1_deriv)

  USE MPI
  USE EZ_PARALLEL_STRUCTS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: order
  INTEGER, INTENT(IN) :: decomp_id
  DOUBLE COMPLEX, INTENT(INOUT) :: spec_dim1_deriv(:,:)

  ! Check for errors in user input.
  CALL ERROR_HANDLING(decomp_id, order)

  PRINT *, "FLAG A"

  ! Run in serial if there is only one processor.
  IF (proc_count .EQ. 1) THEN
     PRINT *, "FLAG B_S"
     CALL SERIAL(order, decomp_id, spec_dim1_deriv)
  ELSE
     PRINT *, "FLAG B_P"
     CALL PARALLEL(order, decomp_id, spec_dim1_deriv)
  END IF

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The serial <tt>EZ_PARALLEL</tt> first-dimension wavenumber array generator
  !! subroutine.
  !> When running in serial, there is no grid decomposition, so we are able to
  !! generate the first-dimension wavenumber array without anything special.
  !
  !> \param[in] order The order of derivative desired.
  !> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
  !! number of "unique" grid decompositions.
  !> \param[inout] spec_dim1_deriv The array to hold the wavenumbers along the
  !! first dimension.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE SERIAL(order, decomp_id, spec_dim1_deriv)

    USE MPI
    USE EZ_PARALLEL_STRUCTS
  
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: order
    INTEGER, INTENT(IN) :: decomp_id
    DOUBLE COMPLEX, INTENT(INOUT) :: spec_dim1_deriv(:,:)
    INTEGER :: row_count_L !< Row count of grid, local to subroutine.
    INTEGER :: col_count_L !< Col count of grid, local to subroutine.
    INTEGER :: i !< Counter for DO loops.
    DOUBLE PRECISION :: min_wavenum !< The minimum wavenumber.
    DOUBLE PRECISION, ALLOCATABLE :: dim1_wavenums(:) !< Array for wavenumbers.

    PRINT *, "FLAG 0"
    col_count_L = grid_decomps(decomp_id)%col_decomp_ovlp(proc_id+1)
    row_count_L = grid_decomps(decomp_id)%row_count_g
    ALLOCATE(dim1_wavenums(row_count_L))

    PRINT *, "FLAG 1"
    
    ! Calculate dim1 wavenumbers.
    ! Positive wavenumbers.
    DO i = 1, row_count_L/2+1
       dim1_wavenums(i) = DBLE(i - 1)
    END DO
    PRINT *, "FLAG 2"
    ! Negative wavenumbers.
    min_wavenum = -(row_count_L+1)/2 ! Minimum wavenumber -1, since DO loop
    ! starts at 1.
    DO i = 1, (row_count_L-1)/2
       dim1_wavenums(i+row_count_L/2+1) = min_wavenum + i
    END DO
    PRINT *, "FLAG 3"

    ! If row_count even and order odd, we have to zero out highest wavenumber
    ! for derivative.
    IF ((MOD(row_count_L, 2) .EQ. 0) .AND. (MOD(order, 2) .EQ. 1)) THEN
       dim1_wavenums(row_count_L/2+1) = 0.0
    END IF
    PRINT *, "FLAG 4"

    ! Fill in the 2-D array. Since the sub-grid is the entire grid, it contains
    ! all of the columns.
    !ALLOCATE(spec_dim1_deriv(row_count_L,col_count_L))
    PRINT *, "FLAG 5"
    DO i = 1, col_count_L
       PRINT *, "i: ", i
       spec_dim1_deriv(:,i) = (DCMPLX(0.0, dim1_wavenums(:)))**(DBLE(order))
    END DO

    ! Deallocate the temporary arrays.
    DEALLOCATE(dim1_wavenums)

  END SUBROUTINE SERIAL

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The parallel <tt>EZ_PARALLEL</tt> first-dimension wavenumber array
  !! generator subroutine.
  !> When running in parallel, outputs the array containing the wavenumbers
  !! corresponding to the sub-grid of the processor calling the subroutine.
  !! Since the sub-grids contain the entire column, this is very similar to
  !! to the serial version of the subroutine.
  !
  !
  !> \param[in] order The order of derivative desired.
  !> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
  !! number of "unique" grid decompositions.
  !> \param[inout] spec_dim1_deriv The array to hold the wavenumbers along the
  !! first dimension.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE PARALLEL(order, decomp_id, spec_dim1_deriv)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: order
    INTEGER, INTENT(IN) :: decomp_id
    DOUBLE COMPLEX, INTENT(INOUT) :: spec_dim1_deriv(:,:)
    INTEGER :: row_count_L !< Row count of grid, local to subroutine.
    INTEGER :: col_count_L !< Col count of grid, local to subroutine.
    INTEGER :: i !< Counter for DO loops.
    DOUBLE PRECISION :: min_wavenum !< The minimum wavenumber.
    DOUBLE PRECISION, ALLOCATABLE :: dim1_wavenums(:)

    col_count_L = grid_decomps(decomp_id)%col_decomp_ovlp(proc_id+1)
    row_count_L = grid_decomps(decomp_id)%row_count_g
    ALLOCATE(dim1_wavenums(row_count_L))

    ! Calculate dim1 wavenumbers.
    ! Positive wavenumbers.
    DO i = 1, row_count_L/2+1
       dim1_wavenums(i) = DBLE(i - 1)
    END DO
    ! Negative wavenumbers.
    min_wavenum = -(row_count_L+1)/2 ! Minimum wavenumber -1, since DO loop
    ! starts at 1.
    DO i = 1, (row_count_L-1)/2
       dim1_wavenums(i+row_count_L/2+1) = min_wavenum + i
    END DO

    ! If row_count even and order odd, we have to zero out highest wavenumber
    ! for derivative.
    IF ((MOD(row_count_L, 2) .EQ. 0) .AND. (MOD(order, 2) .EQ. 1)) THEN
       dim1_wavenums(row_count_L/2+1) = 0.0
    END IF

    ! Fill in the 2-D array.
    !ALLOCATE(spec_dim1_deriv(row_count_L,col_count_L))
    DO i = 1, col_count_L
       spec_dim1_deriv(:,i) = (DCMPLX(0.0, dim1_wavenums(:)))**(DBLE(order))
    END DO

    ! Deallocate the temporary arrays.
    DEALLOCATE(dim1_wavenums)

  END SUBROUTINE PARALLEL

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief
  !> The <tt>EZ_PARALLEL</tt> first-dimension wavenumber array generator
  !! error handling subroutine.
  !> Aborts the program if the user input trigger any of the following errors:
  !! <ul>
  !! <li> order...
  !!    <ol>
  !!    <li> is negative
  !!    </ol>
  !! <li> decomp_id...
  !!    <ol>
  !!    <li> is outside of (1, ..., <tt>decomp_count</tt>)
  !!    </ol>
  !! </ul>
  !
  !> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
  !! number of "unique" grid decompositions.
  !> \param[in] order The order of derivative desired.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE ERROR_HANDLING(order, decomp_id)

    USE MPI
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: order
    INTEGER, INTENT(IN) :: decomp_id
    INTEGER :: error_code !< Error code for MPI_ABORT.
    INTEGER :: ierror
    LOGICAL :: error_flag !< Error flag to stop all processors.

    error_flag = .FALSE.

    ! Check if order is negative.
    IF (order .LT. 1) THEN
       IF (proc_id .EQ. 0) THEN
          PRINT *, "ERROR: Order of spectral derviative in first dimension ", &
               " in grid with decomposition ID ", decomp_ID, " is negative. ", &
               " order: ", order, "." 
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

END SUBROUTINE SPECTRAL_DIM1_DERIVATIVE

