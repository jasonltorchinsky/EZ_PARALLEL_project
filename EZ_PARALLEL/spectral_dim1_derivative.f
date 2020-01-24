!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SPECTRAL_DIM1_DERIVATIVE. Produces the matrix used in calculating the
! derivative of the numerical solution along the first dimension (dim1) of the
! grid.
!
! ARGUMENTS: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - spec_dim1_deriv: The array that will store the matrix used in calculating
! the derivative of the numerical solution along the second dimension of the
! grid (DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len)).
! - order: The order of the derviate desired (INTEGER).
!
! NOTES: - When calculating the dim1 wavenumbers, we use the convention
! [0, ..., dim1_len/2, -dim1_len/1 + 1, ... -1], with FLOOR(dim2_len/2), as
! is assumed in Fortran integer division.
!
! Written By: Jason Turner
! Last Updated: January 12, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE SPECTRAL_DIM1_DERIVATIVE_EZP(dim1_len, dim2_len, overlap, &
& spec_dim1_deriv, order)
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
DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len) :: spec_dim1_deriv

CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, proc_count, ierror)

! Check for errors in user input.
CALL ERROR_HANDLING(order, proc_id)

! Run in serial if there is only one processor.
IF (proc_count .EQ. 1) THEN
  CALL SERIAL(dim1_len, dim2_len, overlap, spec_dim1_deriv, order)
ELSE
  CALL PARALLEL(dim1_len, dim2_len, overlap, spec_dim1_deriv, order)
END IF


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SERIAL. Handles the execution if there is only one processor.
! (SAME AS PARALLEL)
!
! STRUCTURE: 1) Calculate the dim1 wavenumbers.
! 2) If dim1_len even and order odd, zero out the dim1_len/2 wavenumber.
! 3) Populate spec_dim1_deriv array with (sqrt(-1)*dim1 wavenumber)**order.
!
! VARIABLES: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - spec_dim1_deriv: The array that will store the matrix used in calculating
! the derivative of the numerical solution along the first dimension of the
! grid (DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len)).
! - order: The order of the derviate desired (INTEGER).
! - i: Counting index used in DO loops (INTEGER).
! - dim1_wavenums: Array used to store the dim1 wavenumbers (DOUBLE PRECISION,
! DIMENSION(dim1_len)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE SERIAL(dim1_len, dim2_len, overlap, spec_dim1_deriv, order)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& overlap, &
& order, &
& i
DOUBLE PRECISION, DIMENSION(dim1_len) :: dim1_wavenums
DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len) :: spec_dim1_deriv

! Calculate dim1 wavenumbers.
dim1_wavenums = 0.0
dim1_wavenums(2) = 1.0
dim1_wavenums(dim1_len) = -1.0
DO i = 1, dim1_len/2 - 1
  dim1_wavenums(dim1_len - i) = DBLE(-i - 1)
  dim1_wavenums(i + 2) = DBLE(i + 1)
END DO
! If dim1_len even and order odd, we have to zero out highest wavenumber for
! derivative.
IF ((MOD(dim1_len, 2) .EQ. 0) .AND. (MOD(order, 2) .EQ. 1)) THEN
  dim1_wavenums(dim1_len/2 + 1) = 0.0
END IF

DO i = 1, dim2_len
  spec_dim1_deriv(:,i) = (DCMPLX(0.0, dim1_wavenums(:)))**(DBLE(order))
END DO

END SUBROUTINE SERIAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PARALLEL. Handles the execution if there is more than one processor.
! (SAME AS SERIAL)
!
! STRUCTURE: 1) Calculate the dim1 wavenumbers.
! 2) If dim1_len even and order odd, zero out the dim1_len/2 wavenumber.
! 3) Populate spec_dim1_deriv array with (sqrt(-1)*dim1 wavenumber)**order.
!
! VARIABLES: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - spec_dim1_deriv: The array that will store the matrix used in calculating
! the derivative of the numerical solution along the first dimension of the
! grid (DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len)).
! - order: The order of the derviate desired (INTEGER).
! - i: Counting index used in DO loops (INTEGER).
! - dim1_wavenums: Array used to store the dim1 wavenumbers (DOUBLE PRECISION,
! DIMENSION(dim1_len)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE PARALLEL(dim1_len, dim2_len, overlap, spec_dim1_deriv, order)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& overlap, &
& order, &
& i
DOUBLE PRECISION, DIMENSION(dim1_len) :: dim1_wavenums
DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len) :: spec_dim1_deriv

! Calculate dim1 wavenumbers.
dim1_wavenums = 0.0
dim1_wavenums(2) = 1.0
dim1_wavenums(dim1_len) = -1.0
DO i = 1, dim1_len/2 - 1
  dim1_wavenums(dim1_len - i) = DBLE(-i - 1)
  dim1_wavenums(i + 2) = DBLE(i + 1)
END DO
! If dim1_len even and order odd, we have to zero out highest wavenumber for
! derivative.
IF ((MOD(dim1_len, 2) .EQ. 0) .AND. (MOD(order, 2) .EQ. 1)) THEN
  dim1_wavenums(dim1_len/2 + 1) = 0.0
END IF

DO i = 1, dim2_len
    spec_dim1_deriv(:,i) = (DCMPLX(0.0, dim1_wavenums(:)))**(DBLE(order))
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
& proc_id, &
& ierror
LOGICAL :: proc_err_flag, &
& global_err_flag

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

proc_err_flag = .FALSE.
global_err_flag = .FALSE.

! order < 1.
IF (order .LT. 1) THEN
  proc_err_flag = .TRUE.
  PRINT *, 'EZ_PARALLEL: Issue with processor ', proc_id, ' in subroutine ', &
  & 'call SPECTRAL_DIM1_DERIVATIVE_DBLE_CMPLX. order: ', order, ' is less ', &
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

END SUBROUTINE SPECTRAL_DIM1_DERIVATIVE_EZP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
