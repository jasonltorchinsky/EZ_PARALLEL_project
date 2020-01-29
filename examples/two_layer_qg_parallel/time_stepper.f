!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TIMESTEPPER. Contains the time stepping subroutine to progress the simulation
! forward in time.
!
! VARIABLES: - TIME_STEP: The subroutine which drives the model's time
! integration (PUBLIC).
!
! Written By: Jason Turner
! Last Updated: January 28, 2020
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
! STRUCTURE: 1) Checks if output_frew is too large.
! 2) Writes the initial condition to file.
! 3) FFTs the initial condition, sets up the spectral biharmonic viscosity
! operator.
! 4) Steps the simulation forward in time, writing the output at the desired
! frequency.
!
! VARIABLES: - step: The current time step number (INTEGER(qb)).
! - i, j: Counting index used in DO loops (INTEGER(qb)).
! - time: The current time in the simulation (REAL(dp)).
! - dt: The current timestep size, adapted throughout simulation (REAL(dp)).
! - error_toler_0: The error tolerance for the current step of the
! simulation (REAL(dp)).
! - esdirk_coeff: The coefficients for the explicit singly diagonally implicit
! Runge-Kutta (ESDIRK) portion of the additive Runge-Kutta method, used for the
! 'fast' phenomena (COMPLEX(dp), DIMENSION(5,5)).
! - erk_coeff: The coefficients for the explicit Runge-Kutta (ERK) portion of
! the additive Runge-Kutta method, used for the 'slow' phenomena (COMPLEX(dp),
! DIMENSION(6,7)).
! - spec_x_deriv, spec_y_deriv: The array that will store the matrix used in
! calculating the derivative of the numerical solution along the x-, y-direction
! of the grid (COMPLEX(dp), DIMENSION(x_len, y_len)).
! - spec_biharm: The spectral biharmonic operator used in the hyperviscosity
! term of the QG equations (COMPLEX(dp), DIMENSION(x_len, y_len, 2)).
! - freq_pot_vort_grid: The potential vorticity grid in Fourier space
! (COMPLEX(dp), DIMENSION(x_len, y_len, 2)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE TIME_STEP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI ! ADDED TO PARALLEL

USE INITIALIZE
USE JACOBIAN_EKMAN_SHEAR_SOLVE
USE OUTPUT

IMPLICIT NONE

INTEGER(qb) :: step, &
& i, &
& j, &
& proc_id ! ADDED TO PARALLEL
REAL(dp) :: time, &
& dt, &
& error_toler_0
COMPLEX(dp), DIMENSION(5, 5) :: esdirk_coeff
COMPLEX(dp), DIMENSION(6, 7) :: erk_coeff
COMPLEX(dp), DIMENSION(x_len, y_len) :: spec_x_deriv, &
& spec_y_deriv
COMPLEX(dp), DIMENSION(x_len, y_len, 2):: spec_biharm, &
& freq_pot_vort_grid

CALL GET_ID_EZP(proc_id) ! ADDED TO PARALLEL

! Abort if output_freq is too big.
!/* CHANGED FOR PARALLEL
!IF (num_timesteps .LT. output_freq) THEN
!  ERROR STOP 'Output frequency too large. Please reduce.'
!END IF
!*/

!/* ADDED TO PARALLEL
IF ((proc_id .EQ. 0) .AND. (num_timesteps .LT. output_freq)) THEN
  ERROR STOP 'Output frequency too large. Please reduce.'
ELSE IF (num_timesteps .LT. output_freq) THEN
  STOP
END IF
!*/
! Output the initial condition.
time = init_time
dt = init_dt
step = 0_qb
CALL WRITE_OUTPUT(step, time, dt)

! FFT the physical potential vorticity to frequency space for timestepping.
freq_pot_vort_grid = phys_pot_vort_grid
!/* CHANGED FOR PARALLEL
!CALL CFFT2DF(x_len, y_len, freq_pot_vort_grid(:,:,1))
!CALL CFFT2DF(x_len, y_len, freq_pot_vort_grid(:,:,2))
!/*

CALL CFFT2DF_EZP(x_len, y_len, 0_qb, &
  & freq_pot_vort_grid(:,:,1)) ! ADDED TO PARALLEL
CALL CFFT2DF_EZP(x_len, y_len, 0_qb, &
  & freq_pot_vort_grid(:,:,2)) ! ADDED TO PARALLEL

! Set up hyperviscosity for time-stepping.
!/* CHANGED FOR PARALLEL
!CALL SPECTRAL_X_DERIVATIVE(x_len, y_len, spec_x_deriv, 2_qb * biharm_order)
!CALL SPECTRAL_Y_DERIVATIVE(x_len, y_len, spec_y_deriv, 2_qb * biharm_order)
!*/

CALL SPECTRAL_DIM1_DERIVATIVE_EZP(x_len, y_len, 0_qb, spec_x_deriv, &
  & 2_qb * biharm_order) ! ADDED TO PARALLEL
CALL SPECTRAL_DIM2_DERIVATIVE_EZP(x_len, y_len, 0_qb, spec_y_deriv, &
  & 2_qb * biharm_order) ! ADDED TO PARALLEL
spec_biharm(:,:,1) = (-1.0_dp, 0.0_dp)**(REAL(biharm_order + 1_qb, dp)) &
& * biharm_visc_coeff * (spec_x_deriv + spec_y_deriv)
spec_biharm(:,:,2) = spec_biharm(:,:,1)
! Calculate the ARK coefficients.
CALL CALCULATE_ARK_COEFF(erk_coeff, esdirk_coeff)
! Step simulation forward in time.
DO step = 1, num_timesteps
  error_toler_0 = 0.8_dp * error_toler
  CALL TIME_STEP_SCHEME(freq_pot_vort_grid, spec_biharm, time, &
  & error_toler_0, dt, erk_coeff, esdirk_coeff)
  ! Write output at desired frequency.
  IF (MOD(step, output_freq) .EQ. 0_qb) THEN
    ! Output the current step.
    CALL WRITE_OUTPUT(step, time, dt)
  END IF
END DO

END SUBROUTINE TIME_STEP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The numerical scheme for the time step, using an additive Runge-Kutta method
! described by Kennedy, C. A., et. al in "Additive Runge-Kutta Schemes for
! Convection-Diffusion-Reaction Equations" (July 2001). We use
! ARK4(3)6L[2]SA - ERK for the Jacobian, Ekamn friction, and vertical shear
! terms (which act on slow time scales) and ARK4(3)6L[2]SA - ESDIRK for
! the hyperviscosity term (which acts on fast time scales). Since the latter
! is diagonally implicit, the equation for each stage has been rearranged to
! isolate the new value of potential vorticity at that stage. The coefficient
! that arises from this rearrangement is store in stage_coeff.
!
! STRUCTURE: 1) Calculates each of the 5 stages of the ARK method.
! 2) Checks that the one-step error is within the error bound.
! 3) Updates the physical pot vort grid, and adjusts the timestep size.
!
! VARIABLES: - time: The current time in the simulation (REAL(dp)).
! - dt: The current timestep size (REAL(dp)).
! - error_toler_0, error_toler_1: Tolerance parameter for the current, next step
! for adaptive time-stepping (REAL(dp)).
! - max_pot_vort: The potential vorticity bound (REAL(dp)).
! - esdirk_coeff: The coefficients for the explicit singly diagonally implicit
! Runge-Kutta (ESDIRK) portion of the additive Runge-Kutta method, used for the
! 'fast' phenomena (COMPLEX(dp), DIMENSION(5,5)).
! - erk_coeff: The coefficients for the explicit Runge-Kutta (ERK) portion of
! the additive Runge-Kutta method, used for the 'slow' phenomena (COMPLEX(dp),
! DIMENSION(6,7)).
! - freq_pot_vort_grid: The potential vorticity in Fourier space (COMPLEX(dp),
! DIMENSION(x_len, y_len, 2)).
! - spec_biharm: The spectral biharmonic operator used in the hyperviscosity
! term of the QG equations (COMPLEX(dp), DIMENSION(x_len, y_len, 2)).
! - stage_coeff: The multiplicative coefficient common across all stages of the
! ARK method.
! - jacobian_ekman_shear_0, jacobian_ekman_shear_1, jacobian_ekman_shear_2,
! jacobian_ekman_shear_3, jacobian_ekman_shear_4, jacobian_ekman_shear_5: The
! change in potential vorticity at stage 0, 1, 2, 3, 4, 5 due the the Jacobian,
! Ekman friction, and vertical shear terms of the QG equations (COMPLEX(dp),
! DIMENSION(x_len, y_len, 2)).
! - biharm_visc_0, biharm_visc_1, biharm_visc_2, biharm_visc_3, biharm_visc_4,
! biharm_visc_5: The hyperviscosity at stage 0, 1, 2, 3, 4, 5 (COMPLEX(dp),
! DIMENSION(x_len, y_len, 2)).
! - freq_pot_vort_1, freq_pot_vort_2, freq_pot_vort_3, freq_pot_vort_4,
! freq_pot_vort_5: The potential vorticity in Fourier space at stage 1, 2, 3, 4,
! 5 (COMPLEX(dp), DIMENSION(x_len, y_len, 2)).
! - error_control: The array storing the one-step error at each grid point
! (COMPLEX(dp), DIMENSION(x_len, y_len, 2)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RECURSIVE SUBROUTINE TIME_STEP_SCHEME(freq_pot_vort_grid, spec_biharm, time, &
& error_toler_0, dt, erk_coeff, esdirk_coeff)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI ! ADDED TO PARALLEL

USE INITIALIZE
USE JACOBIAN_EKMAN_SHEAR_SOLVE

IMPLICIT NONE

INTEGER(qb) :: proc_id ! ADDED TO PARALLEL
REAL(dp) :: time, &
& dt, &
& error_toler_0, &
!& error_toler_1, & ! CHANGED FOR PARALLEL
!/* ADDED TO PARALLEL
& error_toler_1, &
& error_toler_11, &
& error_toler_12, &
!*/
& max_pot_vort
COMPLEX(dp), DIMENSION(5,5) :: esdirk_coeff
COMPLEX(dp), DIMENSION(6,7) :: erk_coeff
COMPLEX(dp), DIMENSION(x_len, y_len, 2) :: freq_pot_vort_grid, &
& spec_biharm, &
& stage_coeff, &
& jacobian_ekman_shear_0, &
& biharm_visc_0, &
& freq_pot_vort_1, &
& jacobian_ekman_shear_1, &
& biharm_visc_1, &
& freq_pot_vort_2, &
& jacobian_ekman_shear_2, &
& biharm_visc_2, &
& freq_pot_vort_3, &
& jacobian_ekman_shear_3, &
& biharm_visc_3, &
& freq_pot_vort_4, &
& jacobian_ekman_shear_4, &
& biharm_visc_4, &
& freq_pot_vort_5, &
& jacobian_ekman_shear_5, &
& biharm_visc_5, &
& error_control

stage_coeff = 1.0_dp/(1.0_dp - 0.25_dp*dt*spec_biharm)

! First stage of additive RK method.
  ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
  ! at this stage.
jacobian_ekman_shear_0 = (0.0_dp, 0.0_dp)
biharm_visc_0 = (0.0_dp, 0.0_dp)
CALL JACOBIAN_EKMAN_SHEAR(freq_pot_vort_grid, x_len, y_len, dt, &
  & deform_wavenum, rotat_wavenum, vert_shear, ekman_fric_coeff, &
  & jacobian_ekman_shear_0)
biharm_visc_0 = spec_biharm * freq_pot_vort_grid
  ! Calculate potential vorticity in frequency space at this stage.
freq_pot_vort_1 = (0.0_dp, 0.0_dp)
freq_pot_vort_1 = stage_coeff * (freq_pot_vort_grid + CMPLX(dt, 0.0_dp, dp) &
  & * (erk_coeff(1,1) * jacobian_ekman_shear_0 &
    & + esdirk_coeff(1,1) * biharm_visc_0))

! Second stage of additive RK method.
  ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
  ! at this stage.
jacobian_ekman_shear_1 = (0.0_dp, 0.0_dp)
biharm_visc_1 = (0.0_dp, 0.0_dp)
CALL JACOBIAN_EKMAN_SHEAR(freq_pot_vort_1, x_len, y_len, dt, &
  & deform_wavenum, rotat_wavenum, vert_shear, ekman_fric_coeff, &
  & jacobian_ekman_shear_1)
biharm_visc_1 = spec_biharm * freq_pot_vort_1
  ! Calculate potential vorticity in frequency space at this stage.
freq_pot_vort_2 = (0.0_dp, 0.0_dp)
freq_pot_vort_2 = stage_coeff * (freq_pot_vort_grid + CMPLX(dt, 0.0_dp, dp) &
  & * (erk_coeff(1,2) * jacobian_ekman_shear_0 &
    & + erk_coeff(2,2) * jacobian_ekman_shear_1 &
    & + esdirk_coeff(1,2) * biharm_visc_0 &
    & + esdirk_coeff(2,2) * biharm_visc_1))

! Third stage of additive RK method.
  ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
  ! at this stage.
jacobian_ekman_shear_2 = (0.0_dp, 0.0_dp)
biharm_visc_2 = (0.0_dp, 0.0_dp)
CALL JACOBIAN_EKMAN_SHEAR(freq_pot_vort_2, x_len, y_len, dt, &
  & deform_wavenum, rotat_wavenum, vert_shear, ekman_fric_coeff, &
  & jacobian_ekman_shear_2)
biharm_visc_2 = spec_biharm * freq_pot_vort_2
  ! Calculate potential vorticity in frequency space at this stage.
freq_pot_vort_3 = (0.0_dp, 0.0_dp)
freq_pot_vort_3 = stage_coeff * (freq_pot_vort_grid + CMPLX(dt, 0.0_dp, dp) &
  & * (erk_coeff(1,3) * jacobian_ekman_shear_0 &
    & + erk_coeff(2,3) * jacobian_ekman_shear_1 &
    & + erk_coeff(3,3) * jacobian_ekman_shear_2 &
    & + esdirk_coeff(1,3) * biharm_visc_0 &
    & + esdirk_coeff(2,3) * biharm_visc_1 &
    & + esdirk_coeff(3,3) * biharm_visc_2))

! Fourth stage of additive RK method.
  ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
  ! at this stage.
jacobian_ekman_shear_3 = (0.0_dp, 0.0_dp)
biharm_visc_3 = (0.0_dp, 0.0_dp)
CALL JACOBIAN_EKMAN_SHEAR(freq_pot_vort_3, x_len, y_len, dt, &
  & deform_wavenum, rotat_wavenum, vert_shear, ekman_fric_coeff, &
  & jacobian_ekman_shear_3)
biharm_visc_3 = spec_biharm * freq_pot_vort_3
  ! Calculate potential vorticity in frequency space at this stage.
freq_pot_vort_4 = (0.0_dp, 0.0_dp)
freq_pot_vort_4 = stage_coeff * (freq_pot_vort_grid + CMPLX(dt, 0.0_dp, dp) &
  & * (erk_coeff(1,4) * jacobian_ekman_shear_0 &
    & + erk_coeff(2,4) * jacobian_ekman_shear_1 &
    & + erk_coeff(3,4) * jacobian_ekman_shear_2 &
    & + erk_coeff(4,4) * jacobian_ekman_shear_3 &
    & + esdirk_coeff(1,4) * biharm_visc_0 &
    & + esdirk_coeff(2,4) * biharm_visc_1 &
    & + esdirk_coeff(3,4) * biharm_visc_2 &
    & + esdirk_coeff(4,4) * biharm_visc_3))

! Fifth stage of additive RK method.
  ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
  ! at this stage.
jacobian_ekman_shear_4 = (0.0_dp, 0.0_dp)
biharm_visc_4 = (0.0_dp, 0.0_dp)
CALL JACOBIAN_EKMAN_SHEAR(freq_pot_vort_4, x_len, y_len, dt, &
  & deform_wavenum, rotat_wavenum, vert_shear, ekman_fric_coeff, &
  & jacobian_ekman_shear_4)
biharm_visc_4 = spec_biharm * freq_pot_vort_4
  ! Calculate potential vorticity in frequency space at this stage.
freq_pot_vort_5 = (0.0_dp, 0.0_dp)
freq_pot_vort_5 = stage_coeff * (freq_pot_vort_grid + CMPLX(dt, 0.0_dp, dp) &
  & * (erk_coeff(1,5) * jacobian_ekman_shear_0 &
    & + erk_coeff(2,5) * jacobian_ekman_shear_1 &
    & + erk_coeff(3,5) * jacobian_ekman_shear_2 &
    & + erk_coeff(4,5) * jacobian_ekman_shear_3 &
    & + erk_coeff(5,5) * jacobian_ekman_shear_4 &
    & + esdirk_coeff(1,5) * biharm_visc_0 &
    & + esdirk_coeff(3,5) * biharm_visc_2 &
    & + esdirk_coeff(4,5) * biharm_visc_3 &
    & + esdirk_coeff(5,5) * biharm_visc_4))

! Sixth stage of additive RK method.
  ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
  ! at this stage.
jacobian_ekman_shear_5 = (0.0_dp, 0.0_dp)
biharm_visc_5 = (0.0_dp, 0.0_dp)
CALL JACOBIAN_EKMAN_SHEAR(freq_pot_vort_5, x_len, y_len, dt, &
  & deform_wavenum, rotat_wavenum, vert_shear, ekman_fric_coeff, &
  & jacobian_ekman_shear_5)
biharm_visc_5 = spec_biharm * freq_pot_vort_5

! Error control, see the final two rows of the Butcher Tableau for each RK
! method.
  ! Calculate the error for adaptive time stepping.
error_control = (0.0_dp, 0.0_dp)
error_control = erk_coeff(1,6) * (jacobian_ekman_shear_0 + biharm_visc_0) &
& + erk_coeff(3,6) * (jacobian_ekman_shear_2 + biharm_visc_2) &
& + erk_coeff(4,6) * (jacobian_ekman_shear_3 + biharm_visc_3) &
& + erk_coeff(5,6) * (jacobian_ekman_shear_4 + biharm_visc_4) &
& + erk_coeff(6,6) * (jacobian_ekman_shear_5 + biharm_visc_5)

! If the one-step error exceeds the tolerated value, decrease dt by 3/4 and
! try again.
!/* CHANGED FOR PARALLEL
!CALL CFFT2DB(x_len, y_len, error_control(:,:,1))
!CALL CFFT2DB(x_len, y_len, error_control(:,:,2))
!error_toler_1 = dt * MAXVAL(ABS(error_control))
!*/
!*/ ADDED TO PARALLEL
CALL CFFT2DB_EZP(x_len, y_len, 0_qb, error_control(:,:,1))
CALL CFFT2DB_EZP(x_len, y_len, 0_qb, error_control(:,:,2))
CALL MAXVAL_EZP(x_len, y_len, ABS(error_control(:,:,1)), error_toler_11)
CALL MAXVAL_EZP(x_len, y_len, ABS(error_control(:,:,2)), error_toler_12)
error_toler_1 = dt * MAXVAL((/ error_toler_11, error_toler_12 /))
!*/

CALL GET_ID_EZP(proc_id) ! ADDED TO PARALLEL
IF (error_toler_1 .GT. error_toler) THEN
  dt = 0.75_dp * dt
  !/* CHANGED FOR PARALLEL
  !PRINT *, 'Step error too large, retrying...'
  !*/
  !/* ADDDED TO PARALLEL
  IF (proc_id .EQ. 0_qb) THEN
    PRINT *, 'Step error too large, retrying...'
  END IF
  !*/
  CALL TIME_STEP_SCHEME(freq_pot_vort_grid, spec_biharm, time, &
  & error_toler_0, dt, erk_coeff, esdirk_coeff)

END IF

! Successful step, proceed to evaluation.
time = time + dt
!/* CHANGED FOR PARALLEL
!PRINT *, 'Time at current step: ', time, '.'
!PRINT *, 'dt at current step: ', dt, '.'
!PRINT *, 'Error at current step: ', error_toler_1, '.'
!*/
!*/ ADDED TO PARALLEL
IF (proc_id .EQ. 0_qb) THEN
  PRINT *, 'Time at current step: ', time, '.'
  PRINT *, 'dt at current step: ', dt, '.'
  PRINT *, 'Error at current step: ', error_toler_1, '.'
END IF
!*/

freq_pot_vort_grid = freq_pot_vort_grid + CMPLX(dt, 0.0_dp, dp) &
& * (erk_coeff(1,7) * (jacobian_ekman_shear_0 + biharm_visc_0) &
  & + erk_coeff(3,7) * (jacobian_ekman_shear_2 + biharm_visc_2) &
  & + erk_coeff(4,7) * (jacobian_ekman_shear_3 + biharm_visc_3) &
  & + erk_coeff(5,7) * (jacobian_ekman_shear_4 + biharm_visc_4) &
  & + erk_coeff(6,7) * (jacobian_ekman_shear_5 + biharm_visc_5))

! Tranform solution to physical space to check against upper bound of physical
! potential vorticity.
phys_pot_vort_grid = freq_pot_vort_grid
!/* CHANGED FOR PARALLEL
!CALL CFFT2DB(x_len, y_len, phys_pot_vort_grid(:,:,1))
!CALL CFFT2DB(x_len, y_len, phys_pot_vort_grid(:,:,2))
!*/
CALL CFFT2DB_EZP(x_len, y_len, 0_qb, &
  & phys_pot_vort_grid(:,:,1)) ! ADDED TO PARALLEL
CALL CFFT2DB_EZP(x_len, y_len, 0_qb, &
  & phys_pot_vort_grid(:,:,2)) ! ADDED TO PARALLEL
! Zero-out complex part artifacts leftover from inverse FFT.
phys_pot_vort_grid = REAL(phys_pot_vort_grid)
  ! Transform back to frequency space now that complex artifacts are gone.
freq_pot_vort_grid = phys_pot_vort_grid
!/* CHANGED FOR PARALLEL
!CALL CFFT2DF(x_len, y_len, freq_pot_vort_grid(:,:,1))
!CALL CFFT2DF(x_len, y_len, freq_pot_vort_grid(:,:,2))
!*/
CALL CFFT2DF_EZP(x_len, y_len, 0_qb, &
  & freq_pot_vort_grid(:,:,1)) ! ADDED TO PARALLEL
CALL CFFT2DF_EZP(x_len, y_len, 0_qb, &
  & freq_pot_vort_grid(:,:,2)) ! ADDED TO PARALLEL
! Step size adjustment: EPS, PI.3.4.
dt = ((0.75_dp * error_toler)/error_toler_1)**(0.075_dp) &
* (error_toler_0/error_toler_1)**(0.1_dp) * dt
error_toler_0 = error_toler_1
!/* CHANGED FOR PARALLEL
!PRINT *, 'dt for next step: ', dt, '.'
!*/
!*/ ADDED TO PARALLEL
IF (proc_id .EQ. 0_qb) THEN
  PRINT *, 'dt for next step: ', dt, '.'
END IF
!*/

max_pot_vort = MAXVAL(REAL(ABS(phys_pot_vort_grid), dp))
!/* CHANGED FOR PARALLEL
!IF (max_pot_vort .GT. pot_vort_bound) THEN
!  WRITE(*,'(A,F10.3,A,F10.3,A)') 'ERROR: Max phys_pot_vort_grid = ', &
!    & max_pot_vort, ' exceeds pot_vort_bound = ', pot_vort_bound, '.'
!  ERROR STOP
!END IF
!*/
!/* ADDED TO PARALLEL
IF ((proc_id .EQ. 0_qb) .AND. (max_pot_vort .GT. pot_vort_bound)) THEN
  WRITE(*,'(A,F10.3,A,F10.3,A)') 'ERROR: Max phys_pot_vort_grid = ', &
    & max_pot_vort, ' exceeds pot_vort_bound = ', pot_vort_bound, '.'
  ERROR STOP
ELSE IF (max_pot_vort .GT. pot_vort_bound) THEN
  STOP
END IF
!*/

END SUBROUTINE TIME_STEP_SCHEME
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Calculates the coefficients for an additive Runge-Kutta method
! described by Kennedy, C. A., et. al in "Additive Runge-Kutta Schemes for
! Convection-Diffusion-Reaction Equations" (July 2001). We use
! ARK4(3)6L[2]SA - ERK for the Jacobian, Ekman friction, and vertical shear
! terms (which act on slow time scales) and ARK4(3)6L[2]SA - ESDIRK for
! the hyperviscosity term (which acts on fast time scales). Since the latter
! is diagonally implicit, the equation for each stage has been rearranged to
! isolate the new value of potential vorticity at that stage.
!
! STRUCTURE: N/A.
!
! VARIABLES: - erk_coeff: The coefficients for the explicit Runge-Kutta (ERK)
! portion of the additive Runge-Kutta method, used for the 'slow' phenomena
! (COMPLEX(dp), DIMENSION(6,7)).
! - esdirk_coeff: The coefficients for the explicit singly diagonally implicit
! Runge-Kutta (ESDIRK) portion of the additive Runge-Kutta method, used for the
! 'fast' phenomena (COMPLEX(dp), DIMENSION(5,5)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE CALCULATE_ARK_COEFF(erk_coeff, esdirk_coeff)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IMPLICIT NONE

COMPLEX(dp), DIMENSION(5,5) :: esdirk_coeff
COMPLEX(dp), DIMENSION(6,7) :: erk_coeff

! Zero-out the arrays before filling them with the proper coefficients.
erk_coeff(:,:) = CMPLX(0.0_dp, 0.0_dp, dp)
esdirk_coeff(:,:) = CMPLX(0.0_dp, 0.0_dp, dp)

! ARK coefficients for calculating the first stage.
erk_coeff(1,1) = CMPLX((1.0_dp/2.0_dp), 0.0_dp, dp)
esdirk_coeff(1,1) = CMPLX((1.0_dp/4.0_dp), 0.0_dp, dp)

! ARK coefficients for calculating the second stage.
erk_coeff(1,2) = CMPLX((13861.0_dp/62500.0_dp), 0.0_dp, dp)
erk_coeff(2,2) = CMPLX((6889.0_dp/62500.0_dp), 0.0_dp, dp)
esdirk_coeff(1,2) = CMPLX((8611.0_dp/62500.0_dp), 0.0_dp, dp)
esdirk_coeff(2,2) = CMPLX((-1743.0_dp/31250.0_dp), 0.0_dp, dp)

! ARK coefficients for calculating the third stage.
erk_coeff(1,3) = CMPLX((-116923316275.0_dp/2393684061468.0_dp), 0.0_dp, dp)
erk_coeff(2,3) = CMPLX((-2731218467317.0_dp/15368042101831.0_dp), 0.0_dp, dp)
erk_coeff(3,3) = CMPLX((9408046702089.0_dp/11113171139209.0_dp), 0.0_dp, dp)
esdirk_coeff(1,3) = CMPLX((5012029.0_dp/34652500.0_dp), 0.0_dp, dp)
esdirk_coeff(2,3) = CMPLX((-654441.0_dp/2922500.0_dp), 0.0_dp, dp)
esdirk_coeff(3,3) = CMPLX((174375.0_dp/388108.0_dp), 0.0_dp, dp)

! ARK coefficients for calculating the fourth stage.
erk_coeff(1,4) = CMPLX((-451086348788.0_dp/2902428689909.0_dp), 0.0_dp, dp)
erk_coeff(2,4) = CMPLX((-2682348792572.0_dp/7519795681897.0_dp), 0.0_dp, dp)
erk_coeff(3,4) = CMPLX((12662868775082.0_dp/11960479115383.0_dp), 0.0_dp, dp)
erk_coeff(4,4) = CMPLX((3355817975965.0_dp/11060851509271.0_dp), 0.0_dp, dp)
esdirk_coeff(1,4) = CMPLX((15267082809.0_dp/155376265600.0_dp), 0.0_dp, dp)
esdirk_coeff(2,4) = CMPLX((-71443401.0_dp/120774400.0_dp), 0.0_dp, dp)
esdirk_coeff(3,4) = CMPLX((730878875.0_dp/902184768.0_dp), 0.0_dp, dp)
esdirk_coeff(4,4) = CMPLX((2285395.0_dp/8070912.0_dp), 0.0_dp, dp)

! ARK coefficients for calculating the fifth stage.
erk_coeff(1,5) = CMPLX((647845179188.0_dp/3216320057751.0_dp), 0.0_dp, dp)
erk_coeff(2,5) = CMPLX((73281519250.0_dp/8382639484533.0_dp), 0.0_dp, dp)
erk_coeff(3,5) = CMPLX((552539513391.0_dp/3454668386233.0_dp), 0.0_dp, dp)
erk_coeff(4,5) = CMPLX((3354512671639.0_dp/8306763924573.0_dp), 0.0_dp, dp)
erk_coeff(5,5) = CMPLX((4040.0_dp/17871.0_dp), 0.0_dp, dp)
esdirk_coeff(1,5) = CMPLX((82889.0_dp/524892.0_dp), 0.0_dp, dp)
esdirk_coeff(3,5) = CMPLX((15625.0_dp/83664.0_dp), 0.0_dp, dp)
esdirk_coeff(4,5) = CMPLX((69875.0_dp/102672.0_dp), 0.0_dp, dp)
esdirk_coeff(5,5) = CMPLX((-2260.0_dp/8211.0_dp), 0.0_dp, dp)

! ARK coefficients for calculating the error control.
erk_coeff(1,6) = &
CMPLX((82889.0_dp/524892.0_dp) - (4586570599.0_dp/29645900160.0_dp), 0.0_dp, dp)
erk_coeff(3,6) = &
CMPLX((15625.0_dp/83664.0_dp) - (178811875.0_dp/945068544.0_dp), 0.0_dp, dp)
erk_coeff(4,6) = &
CMPLX((69875.0_dp/102672.0_dp) - (814220225.0_dp/1159782912.0_dp), 0.0_dp, dp)
erk_coeff(5,6) = &
CMPLX((-2260.0_dp/8211.0_dp) - (-3700637.0_dp/11593932.0_dp), 0.0_dp, dp)
erk_coeff(6,6) = &
CMPLX((1.0_dp/4.0_dp) - (61727.0_dp/225920.0_dp), 0.0_dp, dp)

! ARK coefficients for calculating the next time step.
erk_coeff(1,7) = CMPLX((82889.0_dp/524892.0_dp), 0.0_dp, dp)
erk_coeff(3,7) = CMPLX((15625.0_dp/83664.0_dp), 0.0_dp, dp)
erk_coeff(4,7) = CMPLX((69875.0_dp/102672.0_dp), 0.0_dp, dp)
erk_coeff(5,7) = CMPLX((-2260.0_dp/8211.0_dp), 0.0_dp, dp)
erk_coeff(6,7) = CMPLX((1.0_dp/4.0_dp), 0.0_dp, dp)

END SUBROUTINE CALCULATE_ARK_COEFF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END MODULE
