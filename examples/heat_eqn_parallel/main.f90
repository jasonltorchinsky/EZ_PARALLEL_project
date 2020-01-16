!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  This is an example of a serial code that our package will be able to
!  paralellize. Its structure is akin to that found in the UCLA LES.
!
!  Variables:
!    - start_time, end_time : Beginning, end time of simulation.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM PARALLEL_HEAT_EQN_SOLVER

IMPLICIT NONE

! Real precision types, a la Metcal et. al (2004, p. 71).
INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6, 37)
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(33, 4931)

REAL(dp) :: start_time, end_time

CALL CPU_TIME(start_time)
CALL MAIN
CALL CPU_TIME(end_time)

WRITE(*,"(A,F10.5,A)") "Execution time: ", end_time - start_time, "."
WRITE(*,"(A)") "Simulation complete. Normal termination..."


CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE MAIN
! This is the main program driver. It calls routines to read the model
! initialization file, configure memory and points. It also initializes the
! file and timesteps.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE INITIALIZE
USE TIME_STEPPER

! ADDED TO PARALLEL.
USE MPI

IMPLICIT NONE

! ADDED TO PARALLEL.
CALL INIT_MPI

CALL INITIALIZE_PARAMETERS

! ADDED TO PARALLEL.
CALL DECOMP_GRID(y_len, 1)

! ADDED TO PARALLEL.
CALL IDENTIFY_REF_POINT(y_len, y_ref, dy, 1)

CALL INITIALIZE_GRID

CALL TIME_STEP

CALL FIN_MPI

RETURN

END SUBROUTINE MAIN

END PROGRAM
