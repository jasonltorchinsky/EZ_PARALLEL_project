!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! WRITE_OUTPUT. Contains the output subroutine. For general output management.
!
! VARIABLES: - WRITE_OUTPUT: The subroutine which writes temp_grid to file.
!
! Written By: Jason Turner
! Last Updated: January 17, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE OUTPUT

IMPLICIT NONE

PRIVATE

! Defines standard integer-, real-precision types.
INCLUDE 'integer_types.h'
INCLUDE 'real_types.h'

PUBLIC :: WRITE_OUTPUT

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The output routine for the temperature grid.
!
! STRUCTURE: N/A.
!
! VARIABLES: - step: The current time step number (INTEGER(qb)).
! - step_parity: Which index of temp_grid should be output (either (:,:,1) or
! (:,:,2)) (INTEGER(qb)).
! - i, j: Countind indices for DO loops (INTEGER(qb)).
! - fname: Output file name (CHARACTER(LEN = 30))
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE WRITE_OUTPUT(step, step_parity)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE INITIALIZE

IMPLICIT NONE

INTEGER(qb) :: step, &
& step_parity, &
& i, &
& j
CHARACTER(LEN = 30) :: fname

WRITE(fname,'(A,I0.8,A)') '.\output_data\out_', step, '.csv'

OPEN(100, file = fname, form = 'formatted')

DO j = 1 , y_len
  DO i = 1, x_len
    WRITE(100,'(E32.16,A,1x)',  ADVANCE = 'NO') temp_grid(i, j, step_parity), &
    & ','
  END DO
  WRITE(100,'(1x)')
END DO
CLOSE(100)

PRINT *, 'Wrote grid to ', fname, '.'

END SUBROUTINE WRITE_OUTPUT
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END MODULE

