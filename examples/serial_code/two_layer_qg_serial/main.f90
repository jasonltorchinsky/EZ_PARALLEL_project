!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! A recreation of the non-dimensional, two-layer quasi-geostrophic (QG)
! equations solver originally supplied to us by Di Qi of NYU. Their code
! solves the equation in the case of equal layers, a rigid lid, and a doubly-
! periodic domain. Their code was originally written in MATLAB.
!
! AUTHOR: Jason Turner, under the advision of Sam Stechmann.
! LAST EDITTED: June 21, 2019.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM SERIAL_TWO_LAYER_QG_EQN_SOLVER

IMPLICIT NONE

! Defines standard integer-, real-precision types.
INCLUDE 'integer_types.h'
INCLUDE 'real_types.h'

REAL(dp) :: t1, t2

CALL CPU_TIME(t1)
CALL MAIN
CALL CPU_TIME(t2)

WRITE(*,"(A,F18.6,A)") "Execution time: ", t2 - t1, "."
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

IMPLICIT NONE

CALL INITIALIZE_PARAMETERS

CALL INITIALIZE_GRID

CALL TIME_STEP

RETURN

END SUBROUTINE MAIN

END PROGRAM
