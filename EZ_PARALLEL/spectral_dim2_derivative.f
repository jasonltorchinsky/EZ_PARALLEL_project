!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SPECTRAL_DIM2_DERIVATIVE. Produces the matrix used in calculating the
! derivative of the numerical solution along the second dimension (dim2) of the
! grid.
!
! ARGUMENTS: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - spec_dim2_deriv: The array that will store the matrix used in calculating
! the derivative of the numerical solution along the second dimension of the
! grid (DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len)).
! - order: The order of the derviate desired (INTEGER).
!
! NOTES: - When calculating the dim2 wavenumbers, we use the convention
! [0, ..., dim2_len/2, -dim2_len/2 + 1, ... -1], with FLOOR(dim2_len/2), as
! is assumed in Fortran integer division.
!
! Written By: Jason Turner
! Last Updated: January 12, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE SPECTRAL_DIM2_DERIVATIVE_EZP(dim1_len, dim2_len, overlap, &
& spec_dim2_deriv, order)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& overlap, &
& order, &
& proc_id, &
& proc_count, &
& ierror
DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len) :: spec_dim2_deriv

CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, proc_count, ierror)

! Check for errors in user input.
CALL ERROR_HANDLING(order, proc_id)

! Run in serial if there is only one processor.
IF (proc_count .EQ. 1) THEN
  CALL SERIAL(dim1_len, dim2_len, overlap, spec_dim2_deriv, order)
ELSE
  CALL PARALLEL(dim1_len, dim2_len, overlap, spec_dim2_deriv, order, &
  & proc_id, proc_count)
END IF


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SERIAL. Handles the execution if there is only one processor.
!
! STRUCTURE: 1) Calculate the dim2 wavenumbers.
! 2) If dim2_len even and order odd, zero out the dim2_len/2 wavenumber.
! 3) Populate spec_dim2_deriv array with (sqrt(-1)*dim2 wavenumber)**order.
!
! VARIABLES: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - spec_dim2_deriv: The array that will store the matrix used in calculating
! the derivative of the numerical solution along the second dimension of the
! grid (DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len)).
! - order: The order of the derviate desired (INTEGER).
! - i: Counting index used in DO loops (INTEGER).
! - dim2_wavenums: Array used to store the dim2 wavenumbers (DOUBLE PRECISION,
! DIMENSION(dim2_len)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE SERIAL(dim1_len, dim2_len, overlap, spec_dim2_deriv, order)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& overlap, &
& order, &
& i
DOUBLE PRECISION, DIMENSION(dim2_len) :: dim2_wavenums
DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len) :: spec_dim2_deriv

! Calculate dim2 wavenumbers.
dim2_wavenums = 0.0
dim2_wavenums(2) = 1.0
dim2_wavenums(dim2_len) = -1.0
DO i = 1, dim2_len/2 - 1
  dim2_wavenums(dim2_len - i) = DBLE(-i - 1)
  dim2_wavenums(i + 2) = DBLE(i + 1)
END DO
! If dim2_len even, we have to zero out highest wavenumber for odd-order
! derivatives.
IF ((MOD(dim2_len, 2) .EQ. 0) .AND. (MOD(order, 2) .EQ. 1)) THEN
  dim2_wavenums(dim2_len/2 + 1) = 0.0
END IF

! Populate spec_dim2_deriv array.
DO i = 1, dim2_len
  spec_dim2_deriv(:,i) = (DCMPLX(0.0, dim2_wavenums(i)))**(DBLE(order))
END DO

END SUBROUTINE SERIAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PARALLEL. Handles the execution if there is more than one processor.
!
! STRUCTURE: 1) Recover the size of the (original) grid.
! 2) Find where the processor's sub-grid (referred to simply as "sub-grid"
!    henceforth) is in the grid.
! 3) Calculate the dim2 wavenumbers for the sub-grid. If dim2_len_total even
! and order odd, zero out the dim2_len/2 wavenumber.
! 4) Populate spec_dim2_deriv array with (sqrt(-1)*dim2 wavenumber)**order.
!
! VARIABLES: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - spec_dim2_deriv: The array that will store the matrix used in calculating
! the derivative of the numerical solution along the second dimension of the
! grid (DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len)).
! - order: The order of the derviate desired (INTEGER).
! - proc_id: Processor ID (INTEGER).
! - proc_count: Total number of processors (INTEGER).
! - dim2_len_raw: Dim2 length of the interior of the sub-grid (INTEGER).
! - dim2_len_total: Dim2 length of the grid (INTEGER).
! - dim2_grid_ind_low, dim2_grid_ind_high: The dim2 indices of the sub-grid in
! the grid (if the sub-grid is indexed locally, these are the global dim2
! indices) (INTEGER).
! - dim2_len_pre: Dim2 length of all processors with lesser ID, i.e., the global
! dim2 index for the non-overlap part of the sub-grid (INTEGER).
! - ierror: Integer for holding MPI error flag IDs (INTEGER).
! - i: Counting index used in DO loops (INTEGER).
! - dim2_len_list: Array of the dim2_len for all processors (INTEGER,
! DIMENSION(proc_count)).
! - dim2_wavenums: Array used to store the dim2 wavenumbers (DOUBLE PRECISION,
! DIMENSION(dim2_len)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE PARALLEL(dim1_len, dim2_len, overlap, spec_dim2_deriv, order, &
& proc_id, proc_count)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& overlap, &
& order, &
& dim2_len_raw, &
& dim2_len_total, &
& dim2_grid_ind_low, &
& dim2_grid_ind_high, &
& dim2_len_pre, &
& proc_id, &
& proc_count, &
& ierror, &
& i
INTEGER, DIMENSION(proc_count) :: dim2_len_list
DOUBLE PRECISION, DIMENSION(dim2_len) :: dim2_wavenums
DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len) :: spec_dim2_deriv

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

! Calculate the dim2 grid indices of the sub-grid (the global coordinates
! of the sub-grid).
dim2_len_pre = SUM(dim2_len_list(1:proc_id))
IF (proc_id .EQ. 0) THEN
  dim2_grid_ind_low = dim2_len_pre + 1
  dim2_grid_ind_high = dim2_len_pre + dim2_len_raw + overlap
ELSE IF (proc_id .EQ. proc_count - 1) THEN
  dim2_grid_ind_low = dim2_len_pre + 1 - overlap
  dim2_grid_ind_high = dim2_len_pre + dim2_len_raw
ELSE
  dim2_grid_ind_low = dim2_len_pre + 1 - overlap
  dim2_grid_ind_high = dim2_len_pre + dim2_len_raw + overlap
END IF


! Calculate dim2 wavenumbers.
! If dim2_len_total even, we have to zero out highest wavenumber for odd-order
! derivatives.
IF (MOD(dim2_len_total, 2) .EQ. 0) THEN
  DO i = 1, dim2_len
    ! Non-negative wavenumbers.
    IF (i+dim2_grid_ind_low-1 .LT. dim2_len_total/2+1) THEN
      dim2_wavenums(i) = DBLE(i+dim2_grid_ind_low-2)
    ! Negative wavenumbers.
    ELSE IF (i+dim2_grid_ind_low-1 .GT. dim2_len_total/2+1) THEN
      dim2_wavenums(i) = DBLE(i-dim2_len_total+dim2_grid_ind_low-2)
    ELSE
      IF (MOD(order, 2) .EQ. 1) THEN
        dim2_wavenums(i) = 0.0
      END IF
    END IF
  END DO
ELSE
  DO i = 1, dim2_len
    ! Non-negative wavenumbers.
    IF (i+dim2_grid_ind_low-1 .LE. dim2_len_total/2+1) THEN
      dim2_wavenums(i) = DBLE(i+dim2_grid_ind_low-2)
    ! Negative wavenumbers.
    ELSE IF (i+dim2_grid_ind_low-1 .GT. dim2_len_total/2+1) THEN
      dim2_wavenums(i) = DBLE(i-dim2_len_total+dim2_grid_ind_low-2)
    END IF
  END DO
END IF

! Populate spec_dim2_deriv array with (sqrt(-1)*dim2 wavenumber)**order.
DO i = 1, dim2_len
  spec_dim2_deriv(:,i) = (DCMPLX(0.0, dim2_wavenums(i)))**(DBLE(order))
END DO

END SUBROUTINE PARALLEL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ERROR_HANDLING. Checks for input error and other possible complications.
!
! ERRORS TO CHECK: - order < 1.
!
! VARIABLES: - order: The order of the derviate desired (INTEGER).
! - proc_id: Processor ID (INTEGER).
! - proc_err_flag - Flag is triggered if processor encounters an error
! (LOGICAL).
! - global_err_flag - Flag is triggered if any processor encouters an error
! (LOGICAL).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE ERROR_HANDLING(order, proc_id)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: order, &
& proc_id
LOGICAL :: proc_err_flag, &
& global_err_flag

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

proc_err_flag = .FALSE.
global_err_flag = .FALSE.

! order < 1.
IF (order .LT. 1) THEN
  proc_err_flag = .TRUE.
  PRINT *, 'EZ_PARALLEL: Issue with processor ', proc_id, ' in subroutine ', &
  & 'call SPECTRAL_DIM2_DERIVATIVE_DBLE_CMPLX. order: ', order, ' is less ', &
  & 'than 1. Please increase the order.'
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

END SUBROUTINE SPECTRAL_DIM2_DERIVATIVE_EZP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
