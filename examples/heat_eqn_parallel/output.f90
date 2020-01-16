!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Definition of grid, allocation of memory and I/O management.
!
!  Variables:
!    - nxp, nyp : Number of x-, y-grid points.
!    - dx, dy : x-, y- grid spacing.
!    - timmin, timmax : Initial and final time of simulation.
!    - dt : Time step size.
!    - CX, CY : Thermal diffusivity in x-, y-directions.
!    - u : Array for grid.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE OUTPUT

IMPLICIT NONE

PRIVATE

PUBLIC :: WRITE_OUTPUT

CONTAINS


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE WRITE_OUTPUT(step, time_index)
!  Writes the given array to a .dat file
!
!  Variables:
!    - step : Step in simulation.
!    - time_index : Use in (:,:,1) or (:,:,2) of temperature_grid.
!    - i, j : x-, y-grid point index.
!    - fname : File name.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE INITIALIZE

IMPLICIT NONE

INTEGER :: step, time_index, i, j, proc_id
CHARACTER*35 :: fname

! ADDED TO PARALLEL.
CALL GET_ID(proc_id)

WRITE(fname,"(A,I0.8,A,I0.4,A)") ".\output_data\out_", step, "_", proc_id, &
&  ".dat"

OPEN(100, file = fname, form = "formatted")

DO j = 1 , y_len
  DO i = 1, x_len
    WRITE(100, "(E32.16, 1x)",  ADVANCE = "NO") &
    &  temperature_grid(i, j, time_index)
  END DO
  WRITE(100, "(1x)")
END DO
CLOSE(100)

WRITE(*,"(A, A, A)") "Wrote grid to ", fname, "."

END SUBROUTINE WRITE_OUTPUT

END MODULE

