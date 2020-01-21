!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TIMESTEPPER. Contains the time stepping subroutine to progress the simulation
! forward in time.
!
! VARIABLES: - TIME_STEP: The subroutine which drives the model's time
! integration (PUBLIC).
!
! Written By: Jason Turner
! Last Updated: January 17, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE TIME_STEPPER

IMPLICIT NONE

PRIVATE

! Defines standard integer-, real-precision types.
INCLUDE 'integer_types.h'
INCLUDE 'real_types.h'

PUBLIC :: TIME_STEP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The main driver for the model's time integration. It calls the subroutine
! TIME_STEP_SCHEME which updates variables and the subroutine WRITE_OUTPUT
! which writes the temperature grid to file.
!
! STRUCTURE: 1) Checks if the simulation parameters violate the CFL condition.
! 2) Checks if output_freq is too large.
! 3) Writes the initial condition to file.
! 4) Steps the simulation forward in time, writing the output at the desired
! frequency.
!
! VARIABLES: - step: The current time step number (INTEGER(qb)).
! - max_step: The maximum number of steps the simulation can take, given
! the value of fin_time (INTEGER(qb)).
! - step_parity: Which index of temp_grid should be updated (either (:,:,1) or
! (:,:,2)) (INTEGER(qb)).
! - time: The current time in the simulation (REAL(dp)).
! - stability_constraint: The CFL condition for the time-stepping scheme
! (REAL(dp)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE TIME_STEP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE INITIALIZE
USE OUTPUT

IMPLICIT NONE

INTEGER(qb) :: step, &
& max_step, &
& step_parity
REAL(dp) :: time, &
& stability_constraint

! Stability constraint derived from our time step method, which is Forward
! Euler in time, and a 1st order centered difference for each spatial
! derivative (using Von Neumann analysis).
stability_constraint = 1.0_dp/(2.0_dp &
  & * (x_diffus/(dx**2.0_dp) + y_diffus/(dy**2.0_dp)))

IF ((dt .GT. stability_constraint)) THEN
  PRINT *, 'WARNING: Time step size (dt) exceeds the CFL parameter. ', &
  & 'Simulation may be unstable.'
END IF

step_parity = 1_qb
step = 0_qb
time = init_time
max_step = (fin_time - init_time)/dt

IF (max_step .LT. output_freq) THEN
  ERROR STOP 'Output frequency too large. Please reduce.'
END IF

! Output the initial condition.
CALL WRITE_OUTPUT(step, step_parity)

! Step simulation forward in time.
DO WHILE (time .LE. fin_time)
  step_parity = 3_qb - step_parity
  CALL TIME_STEP_SCHEME(step_parity)
  time = time + dt
  step = step + 1_qb

  ! Write output at desired frequency.
  IF (MOD(step, output_freq) .EQ. 0) THEN
    CALL WRITE_OUTPUT(step, step_parity)
  END IF
END DO

END SUBROUTINE TIME_STEP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The numerical scheme for the time step, using a Forward Euler scheme in time
! and a first-order center difference in space. We also include a deterministic
! forcing, a relaxation term, and white noise (see Equation (2) of Hottovy, S.,
! & Stechmann, S. (2015)).
!
! STRUCTURE: 1) Gets a random number for the white noise.
! 2) Changes the step parity.
! 3) Steps forward in time.
!
! VARIABLES: - step_parity: Which index of temp_grid should be updated
! (either (:,:,1) or (:,:,2)) (INTEGER(qb)).
! - prev_step_parity: Parity of the previous time step (INTEGER(qb)).
! - white_noise: Random number between 0 and 1 for the stochastic variable
! (REAL(dp), DIMENSION(x_len-1,y_len-1)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE TIME_STEP_SCHEME(step_parity)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE INITIALIZE

IMPLICIT NONE

INTEGER(qb) :: step_parity, &
& prev_step_parity
REAL(dp), DIMENSION(x_len-1,y_len-1) :: white_noise

! Get a random number for the white noise.
CALL RANDOM_NUMBER(white_noise)
white_noise = 2.0_dp * (white_noise - 0.5_dp)

! Update the step parity.
prev_step_parity = 3_qb - step_parity

! Step forward in time.
temp_grid(2:x_len-1, 2:y_len-1, step_parity) = &
& temp_grid(2:x_len-1, 2:y_len-1, prev_step_parity) &
& + (dt * x_diffus / (dx**2.0_dp)) &
& * (temp_grid(3:x_len, 2:y_len-1, prev_step_parity) &
  & - 2.0_dp * temp_grid(2:x_len-1, 2:y_len-1, prev_step_parity) &
  & + temp_grid(1:x_len-2, 2:y_len-1, prev_step_parity)) &
& + (dt * y_diffus / (dy**2.0_dp)) &
& * (temp_grid(2:x_len-1, 3:y_len, prev_step_parity) &
  & - 2.0_dp * temp_grid(2:x_len-1, 2:y_len-1, prev_step_parity) &
  & + temp_grid(2:x_len-1, 1:y_len-2, prev_step_parity)) &
& + deter_force - (1.0_dp/damp_coeff) &
& * (temp_grid(2:x_len-1, 2:y_len-1, step_parity) - relax_trgt) &
& + stoch_mag * white_noise

END SUBROUTINE TIME_STEP_SCHEME
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END MODULE
