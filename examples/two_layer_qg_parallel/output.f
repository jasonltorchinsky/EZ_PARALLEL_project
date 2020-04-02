!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! OUTPUT. Contains the output subroutine. For general output management.
!
! VARIABLES: - WRITE_OUTPUT: The subroutine which writes phys_pot_vort_grid to
! file (PUBLIC).
!
! Written By: Jason Turner
! Last Updated: January 28, 2020
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
! The output routine for the physical pot vort grid.
!
! STRUCTURE: N/A.
!
! VARIABLES: - step: The current time step number (INTEGER(qb)).
! - time: The current time in the simulation (REAL(dp)).
! - dt: The current timestep size, adapted throughout simulation (REAL(dp)).
! - i, j: Counting indices for DO loops (INTEGER(qb)).
! - layer1_file_name, layer2_file_name: Output file name for layer 1, layer 2
! (CHARACTER(LEN = 33))
! - timestep_info_file_name: Output file name for timestep information
! (CHARACTER(LEN = 39))
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE WRITE_OUTPUT(step, time, dt)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI ! ADDED TO PARALLEL

USE INITIALIZE

IMPLICIT NONE

INTEGER(qb) :: step, &
& i, &
& j, &
& proc_id ! ADDED TO PARALLEL
REAL(dp) :: time, &
& dt
CHARACTER(LEN = 38) :: layer1_file_name, & ! CHANGED FOR PARALLEL
& layer2_file_name
CHARACTER(LEN = 44) :: timestep_info_file_name ! CHANGED FOR PARALLEL

CALL GET_ID_EZP(proc_id)

WRITE(layer1_file_name,'(A,I0.8,A,I0.4,A)') './output_data/layer1_', step, &
& '_', proc_id, '.csv' ! CHANGED FOR PARALLEL
WRITE(layer2_file_name,'(A,I0.8,A,I0.4,A)') './output_data/layer2_', step, &
& '_', proc_id, '.csv' ! CHANGED FOR PARALLEL
WRITE(timestep_info_file_name,'(A,I0.8,A,I0.4,A)') './output_data/out_', step, &
& '_', proc_id, '_info.txt' ! CHANGED FOR PARALLEL

OPEN(1001, file = layer1_file_name, form = 'formatted')
OPEN(1002, file = layer2_file_name, form = 'formatted')

DO j = 1 , y_len
  DO i = 1, x_len
    WRITE(1001, '(E32.16)', ADVANCE = 'NO') &
    REAL(phys_pot_vort_grid(i,j,1), dp)
    WRITE(1001, '(A)', ADVANCE = 'NO') ','

    WRITE(1002, '(E32.16)', ADVANCE = 'NO') &
    REAL(phys_pot_vort_grid(i,j,2), dp)
    WRITE(1002, '(A)', ADVANCE = 'NO') ','
  END DO
  WRITE(1001, '(1x)')
  WRITE(1002, '(1x)')
END DO

CLOSE(1001)
CLOSE(1002)


OPEN(1005, file = timestep_info_file_name, form = 'formatted')
WRITE(1005,'(A,I0.8,1x)') 'step = ', step
WRITE(1005,'(A,E32.16,1x)') 'time = ', time
WRITE(1005,'(A,E32.16,1x)') 'dt = ', dt
CLOSE(1005)

PRINT *,  'Processor ', proc_id, ' wrote outputs for step ', &
& step, '.' ! CHANGED FOR PARALLEL

END SUBROUTINE WRITE_OUTPUT
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END MODULE
