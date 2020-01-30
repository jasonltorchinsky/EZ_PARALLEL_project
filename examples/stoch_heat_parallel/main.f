!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STOCH_HEAT_SOLVER_PARALLEL. A finite-difference code for solving the stochastic
! heat equation on a two-dimensional grid.
!
! NOTES: - All output files are written as .csv files to the output_data
! subdirectory.
! - For this simulation, we consider the first dimension (dim1) of the grid to
! be in the x-direction (despite Fortran being a comlumn-major language). This
! is because the dim1 is faster to index through, and thus writing the grid to
! file using this indexing results in dim1 being written as a row.
!
! Written By: Jason Turner
! Last Updated: January 30, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM STOCH_HEAT_SOLVER_PARALLEL

IMPLICIT NONE

! Defines standard integer-, real-precision types.
INCLUDE 'integer_types.h'
INCLUDE 'real_types.h'

REAL(dp) :: start_time, &
& end_time

CALL CPU_TIME(start_time)
CALL MAIN
CALL CPU_TIME(end_time)

WRITE(*,"(A,F10.5,A)") "Execution time: ", end_time - start_time, "."
WRITE(*,*) "STOCH_HEAT_SOLVER_PARALLEL execution complete. ', &
& 'Normal termination..."

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Main program.
!
! STRUCTURE: 1) Calls the subroutines that initialize parameters, intializes the
! grid, and evolves the grid forward in time.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE MAIN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI ! ADDED TO PARALLEL

USE INITIALIZE
USE TIME_STEPPER

IMPLICIT NONE

CALL INIT_MPI_EZP ! ADDED TO PARALLEL

CALL INITIALIZE_PARAMETERS

CALL DECOMP_GRID_EZP(y_len, 1_qb) ! ADDED TO PARALLEL
CALL IDENTIFY_REF_POINT_EZP(y_len, y_ref, dy, 1_qb) ! ADDED TO PARALLEL

CALL INITIALIZE_GRID

CALL TIME_STEP

CALL FIN_MPI_EZP ! ADDED TO PARALLEL

RETURN

END SUBROUTINE MAIN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END PROGRAM
