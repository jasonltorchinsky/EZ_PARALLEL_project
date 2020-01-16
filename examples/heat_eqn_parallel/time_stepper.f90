!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Progresses the simulation forward in time, outputting the temperature grid
!  desired frequency.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE TIME_STEPPER

IMPLICIT NONE

PRIVATE

! Real precision types, a la Metcal et. al (2004, p. 71).
INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6, 37)
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(33, 4931)

PUBLIC :: TIME_STEP

CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE TIME_STEP
!  The main driver for the model's time integration. It calls the routine
!  tstep, which updates variables. It writes output files.
!
!  Variables:
!    - step : Step in simulation.
!    - max_step : Steps required to finish simulation.
!    - time : Curent time of simulation.
!    - stability_constraint : Stability condition for scheme.
!    - step_parity : Fill in (:,:,1) or (:,:,2) of temperature_grid.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE INITIALIZE
USE OUTPUT

IMPLICIT NONE

INTEGER :: step, max_step, step_parity
REAL(dp) :: time, stability_constraint

! Stability constraint derived from our time step method, which is Forward
! Euler in time, and a 1st order centered difference for each spatial
! derivative (using Von Neumann analysis).
stability_constraint = 1.0_dp/(2.0_dp * (x_diffus/(dx**2) + y_diffus/(dy**2)))

IF ((dt .GT. stability_constraint)) THEN
  WRITE(*,"(A,A)") "WARNING: Time step is too large numerical ", &
  & "scheme is unstable."
END IF

step_parity = 1
step = 0
time = init_time
max_step = (fin_time - init_time)/dt

IF (max_step .LT. output_freq) THEN
  ERROR STOP "Output frequency too large. Please reduce."
END IF

! Output the initial condition.
CALL WRITE_OUTPUT(step, step_parity)

! Step simulation forward in time.
DO WHILE (time .LE. fin_time)
  step_parity = 3 - step_parity
  CALL TIME_STEP_SCHEME(step_parity)
  time = time + dt
  step = step + 1

  ! Write output at desired frequency.
  IF (MOD(step, output_freq) .EQ. 0) THEN
    CALL WRITE_OUTPUT(step, step_parity)
  END IF
END DO

END SUBROUTINE TIME_STEP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE TIME_STEP_SCHEME(time_index)
!  The numerical scheme for the time step, using a Forward Euler scheme in time
!  and a first-order center difference in space.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE INITIALIZE

IMPLICIT NONE

INTEGER :: time_index, prev_time_index

prev_time_index = 3 - time_index

temperature_grid(2:x_len-1, 2:y_len-1, time_index) = &
&  temperature_grid(2:x_len-1, 2:y_len-1, prev_time_index) &
&  + (dt * x_diffus / (dx**2)) &
&  * (temperature_grid(3:x_len, 2:y_len-1, prev_time_index) &
&     - 2 * temperature_grid(2:x_len-1, 2:y_len-1, prev_time_index) &
&     + temperature_grid(1:x_len-2, 2:y_len-1, prev_time_index)) &
&  + (dt * y_diffus / (dy**2)) &
&  * (temperature_grid(2:x_len-1, 3:y_len, prev_time_index) &
&     - 2 * temperature_grid(2:x_len-1, 2:y_len-1, prev_time_index) &
&     + temperature_grid(2:x_len-1, 1:y_len-2, prev_time_index))

! ADDED TO PARALLEL.
CALL SHARE_SUBGRID_BOUNDARIES(x_len, y_len, 1, &
&  temperature_grid(:,:,time_index))

END SUBROUTINE TIME_STEP_SCHEME

END MODULE
