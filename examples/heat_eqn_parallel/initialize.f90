!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Initializaes the simultion, including reading the input file and setting
!  the initial conditions for the grid.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE INITIALIZE

IMPLICIT NONE

PRIVATE

! Real precision types, a la Metcal et. al (2004, p. 71).
INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6, 37)
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(33, 4931)

INTEGER, PUBLIC :: x_len, y_len, output_freq
REAL(dp), PUBLIC :: x_ref, y_ref, dx, dy, dt, init_time, fin_time, &
  & x_diffus, y_diffus
REAL(dp), PUBLIC, ALLOCATABLE :: temperature_grid(:, :, :)
PUBLIC :: INITIALIZE_PARAMETERS, INITIALIZE_GRID

CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE INITIALIZE_PARAMETERS
!  Reads the namelist containing the parameters for the simulation.
!
!  Variables:
!    - x_len, y_len : Number of x-, y-grid points.
!    - dx, dy : x-, y- grid spacing.
!    - dt : Time step size.
!    - x_diffus, y_diffus : Thermal diffusivity in x-, y-directions.
!    - x_ref, y_ref : Position of reference grid point (bottom-left).
!    - init_time, fin_time : Initial and final time of simulation.
!    - output_freq : Output frequency of the simulation.
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IMPLICIT NONE

NAMELIST /model/ x_len, y_len, dx, dy, dt, x_diffus, y_diffus, &
  & x_ref, y_ref, init_time, fin_time, output_freq

OPEN(1000, file = "NAMELIST")
READ(1000, nml = model)
CLOSE(1000)

END SUBROUTINE INITIALIZE_PARAMETERS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE INITIALIZE_GRID
!  Sets the initial conditions for the grid using a function. Note that this
!  function must be defined here.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IMPLICIT NONE

INTEGER :: i, j

! Have temperature grid contain two time steps (last index).
ALLOCATE(temperature_grid(x_len, y_len, 2))
temperature_grid(:,:,1) = 0.0_dp

! Fill in the initial condition for the temperature grid.
DO j = 1, y_len
  DO i = 1, x_len
    temperature_grid(i, j, 1) = &
      & INITIAL_CONDITION(x_ref + i * dx, y_ref + j * dy)
  END DO
END DO

! Zeros out boundary.
! Faster to assign value to array which is contiguous in memory all at once.
temperature_grid(:, 1, 1) = 0.0
temperature_grid(:, y_len, 1) = 0.0
! Faster to assign value to array which is not contiguous in memory with DO.
DO j = 1, y_len
  temperature_grid(1, j, 1) = 0.0
  temperature_grid(x_len, j, 1) = 0.0
END DO

! ADDED TO PARALLEL.
CALL SHARE_SUBGRID_BOUNDARIES(x_len, y_len, 1, temperature_grid(:,:,1))

temperature_grid(:,:,2) = temperature_grid(:,:,1)

END SUBROUTINE INITIALIZE_GRID

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REAL(dp) FUNCTION INITIAL_CONDITION(x_pos, y_pos) RESULT(output)
!  Function to set the initial condition based on grid point position.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IMPLICIT NONE

REAL(dp), INTENT(in) :: x_pos, y_pos

! f(x, y) = x + 2 * y.
output = x_pos + 2.0_dp * y_pos

END FUNCTION INITIAL_CONDITION

END MODULE
