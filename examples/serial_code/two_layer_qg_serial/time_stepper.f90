!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Progresses the simulation forward in time, outputting the potential
!  vorticity grid at desired frequency.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE TIME_STEPPER

IMPLICIT NONE

PRIVATE

! Defines standard integer-, real-precision types.
INCLUDE 'integer_types.h'
INCLUDE 'real_types.h'

PUBLIC :: TIME_STEP

CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE TIME_STEP
!  The main driver for the model's time integration. It calls the routine
!  tstep, which updates variables. It writes output files.
!
!  Variables:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE INITIALIZE
USE OUTPUT

IMPLICIT NONE

INTEGER(qb) :: timestep, i, j
REAL(dp) :: time, dt, error_toler_0
REAL(dp), DIMENSION(:), ALLOCATABLE :: wavenumbers_1d
REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: spectral_biharm
COMPLEX(dp), DIMENSION(xy_len, xy_len, 2):: frequency_pot_vorticity_grid
COMPLEX(dp), DIMENSION(5,5) :: esdirk_coeff
COMPLEX(dp), DIMENSION(6,7) :: erk_coeff

! Abort if output_freq is too big.
IF (num_timesteps .LT. output_freq) THEN
  ERROR STOP "Output frequency too large. Please reduce."
END IF

! Output the initial condition.
time = init_time
dt = init_dt
timestep = 0_qb
CALL WRITE_OUTPUT(timestep, time, init_dt)

! FFT the physical potential vorticity to frequency space for timestepping.
frequency_pot_vorticity_grid = physical_pot_vorticity_grid
CALL CFFT2DF(xy_len, xy_len, frequency_pot_vorticity_grid(:,:,1))
CALL CFFT2DF(xy_len, xy_len, frequency_pot_vorticity_grid(:,:,2))

! Set up hyperviscosity for time-stepping.
! Store 1D wavenumbers to make hypervisocsity operator.
ALLOCATE(wavenumbers_1d(xy_len))
wavenumbers_1d(:) = 0.0_dp
wavenumbers_1d(xy_len) = -1.0_dp
DO i = 1, xy_len/2
  wavenumbers_1d(i + 1) = wavenumbers_1d(i) + 1.
  wavenumbers_1d(xy_len - i) = wavenumbers_1d(xy_len - i + 1) - 1.
END DO
! Create the spectral biharmonic operator for hypervisocsity.
ALLOCATE(spectral_biharm(xy_len, xy_len, 2))
spectral_biharm(:,:,:) = 0.0_dp
DO j = 1, xy_len
  DO i = 1, xy_len
    spectral_biharm(i, j, :) = -1.0_dp * biharm_visc_coeff * &
    (wavenumbers_1d(i)**2 + wavenumbers_1d(j)**2)**biharm_order
  END DO
END DO
DEALLOCATE(wavenumbers_1d)

! Calculate the ARK coefficients.
CALL CALCULATE_ARK_COEFF(erk_coeff, esdirk_coeff)

! Step simulation forward in time.
DO timestep = 1, num_timesteps
  error_toler_0 = 0.8_dp * error_toler
  CALL TIME_STEP_SCHEME(frequency_pot_vorticity_grid, spectral_biharm, time, &
  error_toler_0, dt, erk_coeff, esdirk_coeff)

  ! Write output at desired frequency.
  IF (MOD(timestep, output_freq) .EQ. 0) THEN
    ! Output the current timestep.
    CALL WRITE_OUTPUT(timestep, time, dt)
  END IF
END DO

END SUBROUTINE TIME_STEP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE TIME_STEP_SCHEME(freq_pot_vort_grid, spectral_biharm, &
time, error_toler_0, dt, erk_coeff, esdirk_coeff)
! The numerical scheme for the time step, using an additive Runge-Kutta method
! described by Kennedy, C. A., et. al in "Additive Runge-Kutta Schemes for
! Convection-Diffusion-Reaction Equations" (July 2001). We use
! ARK4(3)6L[2]SA - ERK for the Jacobian, Ekamn friction, and vertical shear
! terms (which act on slow time scales) and ARK4(3)6L[2]SA - ESDIRK for
! the hyperviscosity term (which acts on fast time scales). Since the latter
! is diagonally implicit, the equation for each stage has been rearranged to
! isolate the new value of potential vorticity at that stage. The coefficient
! that arises from this rearrangement is store in stage_coeff.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE INITIALIZE
USE JACOBIAN_EKMAN_SHEAR_SOLVE

IMPLICIT NONE

REAL(dp) :: time, dt, error_toler_0, error_toler_1, max_pot_vorticity
REAL(dp), DIMENSION(xy_len, xy_len, 2) :: spectral_biharm
COMPLEX(dp), DIMENSION(xy_len, xy_len, 2) :: freq_pot_vort_grid, &
stage_coeff, jacobian_ekman_shear_0, biharm_visc_0, frequency_pot_vorticity_1, &
jacobian_ekman_shear_1, biharm_visc_1, frequency_pot_vorticity_2, &
jacobian_ekman_shear_2, biharm_visc_2, frequency_pot_vorticity_3, &
jacobian_ekman_shear_3, biharm_visc_3, frequency_pot_vorticity_4, &
jacobian_ekman_shear_4, biharm_visc_4, frequency_pot_vorticity_5, &
jacobian_ekman_shear_5, biharm_visc_5, error_control
COMPLEX(dp), DIMENSION(5,5) :: esdirk_coeff
COMPLEX(dp), DIMENSION(6,7) :: erk_coeff


101 stage_coeff = CMPLX(1.0_dp/(1.0_dp - 0.25_dp*dt*spectral_biharm), 0.0_dp, dp)

! First stage of additive RK method.
  ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
  ! at this stage.
jacobian_ekman_shear_0 = (0.0_dp, 0.0_dp)
biharm_visc_0 = (0.0_dp, 0.0_dp)
CALL JACOBIAN_EKMAN_SHEAR(freq_pot_vort_grid, xy_len, dt, &
deform_wavenumber, rotat_wavenumber, vertical_shear, ekman_friction_coeff, &
jacobian_ekman_shear_0)
biharm_visc_0 = CMPLX(spectral_biharm, 0.0_dp, dp) * freq_pot_vort_grid
  ! Store RK coefficients for this stage.
!erk_coeff_0 = CMPLX((1.0_dp/2.0_dp), 0.0_dp, dp)
!esdirk_coeff_0 = CMPLX((1.0_dp/4.0_dp), 0.0_dp, dp)
  ! Calculate potential vorticity in frequency space at this stage.
frequency_pot_vorticity_1 = (0.0_dp, 0.0_dp)
frequency_pot_vorticity_1 = stage_coeff &
* (freq_pot_vort_grid + CMPLX(dt, 0.0_dp, dp) &
* (erk_coeff(1,1) * jacobian_ekman_shear_0 &
+ esdirk_coeff(1,1) * biharm_visc_0) )

! Second stage of additive RK method.
  ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
  ! at this stage.
jacobian_ekman_shear_1 = (0.0_dp, 0.0_dp)
biharm_visc_1 = (0.0_dp, 0.0_dp)
CALL JACOBIAN_EKMAN_SHEAR(frequency_pot_vorticity_1, xy_len, dt, &
deform_wavenumber, rotat_wavenumber, vertical_shear, ekman_friction_coeff, &
jacobian_ekman_shear_1)
biharm_visc_1 = CMPLX(spectral_biharm, 0.0_dp, dp) &
* frequency_pot_vorticity_1
  ! Store RK coefficients for this stage.
!erk_coeff_0 = CMPLX((13861.0_dp/62500.0_dp), 0.0_dp, dp)
!erk_coeff_1 = CMPLX((6889.0_dp/62500.0_dp), 0.0_dp, dp)
!esdirk_coeff_0 = CMPLX((8611.0_dp/62500.0_dp), 0.0_dp, dp)
!esdirk_coeff_1 = CMPLX((-1743.0_dp/31250.0_dp), 0.0_dp, dp)
  ! Calculate potential vorticity in frequency space at this stage.
frequency_pot_vorticity_2 = (0.0_dp, 0.0_dp)
frequency_pot_vorticity_2 = stage_coeff &
* (freq_pot_vort_grid + CMPLX(dt, 0.0_dp, dp) &
* (erk_coeff(1,2) * jacobian_ekman_shear_0 &
+ erk_coeff(2,2) * jacobian_ekman_shear_1 &
+ esdirk_coeff(1,2) * biharm_visc_0 &
+ esdirk_coeff(2,2) * biharm_visc_1))

! Third stage of additive RK method.
  ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
  ! at this stage.
jacobian_ekman_shear_2 = (0.0_dp, 0.0_dp)
biharm_visc_2 = (0.0_dp, 0.0_dp)
CALL JACOBIAN_EKMAN_SHEAR(frequency_pot_vorticity_2, xy_len, dt, &
deform_wavenumber, rotat_wavenumber, vertical_shear, ekman_friction_coeff, &
jacobian_ekman_shear_2)
biharm_visc_2 = CMPLX(spectral_biharm, 0.0_dp, dp) &
* frequency_pot_vorticity_2
  ! Store RK coefficients for this stage.
!erk_coeff_0 = CMPLX((-116923316275.0_dp/2393684061468.0_dp), 0.0_dp, dp)
!erk_coeff_1 = CMPLX((-2731218467317.0_dp/15368042101831.0_dp), 0.0_dp, dp)
!erk_coeff_2 = CMPLX((9408046702089.0_dp/11113171139209.0_dp), 0.0_dp, dp)
!esdirk_coeff_0 = CMPLX((5012029.0_dp/34652500.0_dp), 0.0_dp, dp)
!esdirk_coeff_1 = CMPLX((-654441.0_dp/2922500.0_dp), 0.0_dp, dp)
!esdirk_coeff_2 = CMPLX((174375.0_dp/388108.0_dp), 0.0_dp, dp)
  ! Calculate potential vorticity in frequency space at this stage.
frequency_pot_vorticity_3 = (0.0_dp, 0.0_dp)
frequency_pot_vorticity_3 = stage_coeff &
* (freq_pot_vort_grid + CMPLX(dt, 0.0_dp, dp) &
* (erk_coeff(1,3) * jacobian_ekman_shear_0 &
+ erk_coeff(2,3) * jacobian_ekman_shear_1 &
+ erk_coeff(3,3) * jacobian_ekman_shear_2 &
+ esdirk_coeff(1,3) * biharm_visc_0 &
+ esdirk_coeff(2,3) * biharm_visc_1 &
+ esdirk_coeff(3,3) * biharm_visc_2))

! Fourth stage of additive RK method.
  ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
  ! at this stage.
jacobian_ekman_shear_3 = (0.0_dp, 0.0_dp)
biharm_visc_3 = (0.0_dp, 0.0_dp)
CALL JACOBIAN_EKMAN_SHEAR(frequency_pot_vorticity_3, xy_len, dt, &
deform_wavenumber, rotat_wavenumber, vertical_shear, ekman_friction_coeff, &
jacobian_ekman_shear_3)
biharm_visc_3 = CMPLX(spectral_biharm, 0.0_dp, dp) &
* frequency_pot_vorticity_3
  ! Store RK coefficients for this stage.
!erk_coeff_0 = CMPLX((-451086348788.0_dp/2902428689909.0_dp), 0.0_dp, dp)
!erk_coeff_1 = CMPLX((-2682348792572.0_dp/7519795681897.0_dp), 0.0_dp, dp)
!erk_coeff_2 = CMPLX((12662868775082.0_dp/11960479115383.0_dp), 0.0_dp, dp)
!erk_coeff_3 = CMPLX((3355817975965.0_dp/11060851509271.0_dp), 0.0_dp, dp)
!esdirk_coeff_0 = CMPLX((15267082809.0_dp/155376265600.0_dp), 0.0_dp, dp)
!esdirk_coeff_1 = CMPLX((-71443401.0_dp/120774400.0_dp), 0.0_dp, dp)
!esdirk_coeff_2 = CMPLX((730878875.0_dp/902184768.0_dp), 0.0_dp, dp)
!esdirk_coeff_3 = CMPLX((2285395.0_dp/8070912.0_dp), 0.0_dp, dp)
  ! Calculate potential vorticity in frequency space at this stage.
frequency_pot_vorticity_4 = (0.0_dp, 0.0_dp)
frequency_pot_vorticity_4 = stage_coeff &
* (freq_pot_vort_grid + CMPLX(dt, 0.0_dp, dp) &
* (erk_coeff(1,4) * jacobian_ekman_shear_0 &
+ erk_coeff(2,4) * jacobian_ekman_shear_1 &
+ erk_coeff(3,4) * jacobian_ekman_shear_2 &
+ erk_coeff(4,4) * jacobian_ekman_shear_3 &
+ esdirk_coeff(1,4) * biharm_visc_0 &
+ esdirk_coeff(2,4) * biharm_visc_1 &
+ esdirk_coeff(3,4) * biharm_visc_2 &
+ esdirk_coeff(4,4) * biharm_visc_3))

! Fifth stage of additive RK method.
  ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
  ! at this stage.
jacobian_ekman_shear_4 = (0.0_dp, 0.0_dp)
biharm_visc_4 = (0.0_dp, 0.0_dp)
CALL JACOBIAN_EKMAN_SHEAR(frequency_pot_vorticity_4, xy_len, dt, &
deform_wavenumber, rotat_wavenumber, vertical_shear, ekman_friction_coeff, &
jacobian_ekman_shear_4)
biharm_visc_4 = CMPLX(spectral_biharm, 0.0_dp, dp) &
* frequency_pot_vorticity_4
  ! Store RK coefficients for this stage.
!erk_coeff_0 = CMPLX((647845179188.0_dp/3216320057751.0_dp), 0.0_dp, dp)
!erk_coeff_1 = CMPLX((73281519250.0_dp/8382639484533.0_dp), 0.0_dp, dp)
!erk_coeff_2 = CMPLX((552539513391.0_dp/3454668386233.0_dp), 0.0_dp, dp)
!erk_coeff_3 = CMPLX((3354512671639.0_dp/8306763924573.0_dp), 0.0_dp, dp)
!erk_coeff_4 = CMPLX((4040.0_dp/17871.0_dp), 0.0_dp, dp)
!esdirk_coeff_0 = CMPLX((82889.0_dp/524892.0_dp), 0.0_dp, dp)
!esdirk_coeff_2 = CMPLX((15625.0_dp/83664.0_dp), 0.0_dp, dp)
!esdirk_coeff_3 = CMPLX((69875.0_dp/102672.0_dp), 0.0_dp, dp)
!esdirk_coeff_4 = CMPLX((-2260.0_dp/8211.0_dp), 0.0_dp, dp)
  ! Calculate potential vorticity in frequency space at this stage.
frequency_pot_vorticity_5 = (0.0_dp, 0.0_dp)
frequency_pot_vorticity_5 = stage_coeff &
* (freq_pot_vort_grid + CMPLX(dt, 0.0_dp, dp) &
* (erk_coeff(1,5) * jacobian_ekman_shear_0 &
+ erk_coeff(2,5) * jacobian_ekman_shear_1 &
+ erk_coeff(3,5) * jacobian_ekman_shear_2 &
+ erk_coeff(4,5) * jacobian_ekman_shear_3 &
+ erk_coeff(5,5) * jacobian_ekman_shear_4 &
+ esdirk_coeff(1,5) * biharm_visc_0 &
+ esdirk_coeff(3,5) * biharm_visc_2 &
+ esdirk_coeff(4,5) * biharm_visc_3 &
+ esdirk_coeff(5,5) * biharm_visc_4))

! Sixth stage of additive RK method.
  ! Calculate Jacobian, Ekman friction, vertical shear, and biharm. viscosity
  ! at this stage.
jacobian_ekman_shear_5 = (0.0_dp, 0.0_dp)
biharm_visc_5 = (0.0_dp, 0.0_dp)
CALL JACOBIAN_EKMAN_SHEAR(frequency_pot_vorticity_5, xy_len, dt, &
deform_wavenumber, rotat_wavenumber, vertical_shear, ekman_friction_coeff, &
jacobian_ekman_shear_5)
biharm_visc_5 = CMPLX(spectral_biharm, 0.0_dp, dp) &
* frequency_pot_vorticity_5

! Error control, see the final two rows of the Butcher Tableau for each RK
! method.
  ! Store RK coefficients for checking the error. Since the ERK and ESDIRK
  ! coefficients match, we only need to store the ERK ones
!erk_coeff_0 = &
!CMPLX((82889.0_dp/524892.0_dp) - (4586570599.0_dp/29645900160.0_dp), 0.0_dp, dp)
!erk_coeff_2 = &
!CMPLX((15625.0_dp/83664.0_dp) - (178811875.0_dp/945068544.0_dp), 0.0_dp, dp)
!erk_coeff_3 = &
!CMPLX((69875.0_dp/102672.0_dp) - (814220225.0_dp/1159782912.0_dp), 0.0_dp, dp)
!erk_coeff_4 = &
!CMPLX((-2260.0_dp/8211.0_dp) - (-3700637.0_dp/11593932.0_dp), 0.0_dp, dp)
!erk_coeff_5 = &
!CMPLX((1.0_dp/4.0_dp) - (61727.0_dp/225920.0_dp), 0.0_dp, dp)
  ! Calculate the error for adaptive time stepping.
error_control = (0.0_dp, 0.0_dp)
error_control = &
erk_coeff(1,6) * (jacobian_ekman_shear_0 + biharm_visc_0) &
+ erk_coeff(3,6) * (jacobian_ekman_shear_2 + biharm_visc_2) &
+ erk_coeff(4,6) * (jacobian_ekman_shear_3 + biharm_visc_3) &
+ erk_coeff(5,6) * (jacobian_ekman_shear_4 + biharm_visc_4) &
+ erk_coeff(6,6) * (jacobian_ekman_shear_5 + biharm_visc_5)

! If the one-step error exceeds the tolerated value, decrease dt by 3/4 and
! try again.
CALL CFFT2DB(xy_len, xy_len, error_control(:,:,1))
CALL CFFT2DB(xy_len, xy_len, error_control(:,:,2))
PRINT *, "Serial error parameter: ", MAXVAL(ABS(error_control))
error_toler_1 = dt * MAXVAL(ABS(error_control))

IF (error_toler_1 .GT. error_toler) THEN
  dt = 0.75_dp * dt
  PRINT *, "Step error too large, retrying..."
  GO TO 101
END IF

! Successful step, proceed to evaluation.
time = time + dt
PRINT *, "Time at current step: ", time, "."
PRINT *, "dt at current step: ", dt, "."
PRINT *, "Error at current step: ", error_toler_1, "."

  ! Store RK coefficients for calculating the next time step. Since the ERK and
  ! ESDIRK coefficients match, we only need to store the ERK ones
!erk_coeff_0 = CMPLX((82889.0_dp/524892.0_dp), 0.0_dp, dp)
!erk_coeff_2 = CMPLX((15625.0_dp/83664.0_dp), 0.0_dp, dp)
!erk_coeff_3 = CMPLX((69875.0_dp/102672.0_dp), 0.0_dp, dp)
!erk_coeff_4 = CMPLX((-2260.0_dp/8211.0_dp), 0.0_dp, dp)
!erk_coeff_5 = CMPLX((1.0_dp/4.0_dp), 0.0_dp, dp)
freq_pot_vort_grid = &
freq_pot_vort_grid + CMPLX(dt, 0.0_dp, dp) &
* (erk_coeff(1,7) * (jacobian_ekman_shear_0 + biharm_visc_0) &
+ erk_coeff(3,7) * (jacobian_ekman_shear_2 + biharm_visc_2) &
+ erk_coeff(4,7) * (jacobian_ekman_shear_3 + biharm_visc_3) &
+ erk_coeff(5,7) * (jacobian_ekman_shear_4 + biharm_visc_4) &
+ erk_coeff(6,7) * (jacobian_ekman_shear_5 + biharm_visc_5))

! Tranform solution to physical space to check against upper bound of physical
! potential vorticity.
physical_pot_vorticity_grid = freq_pot_vort_grid
CALL CFFT2DB(xy_len, xy_len, physical_pot_vorticity_grid(:,:,1))
CALL CFFT2DB(xy_len, xy_len, physical_pot_vorticity_grid(:,:,2))

! Zero-out complex part artifacts leftover from inverse FFT.
physical_pot_vorticity_grid = REAL(physical_pot_vorticity_grid)
  ! Transform back to frequency space now that complex artifacts are gone.
freq_pot_vort_grid = physical_pot_vorticity_grid
CALL CFFT2DF(xy_len, xy_len, freq_pot_vort_grid(:,:,1))
CALL CFFT2DF(xy_len, xy_len, freq_pot_vort_grid(:,:,2))

! Step size adjustment: EPS, PI.3.4.
dt = ((0.75_dp * error_toler)/error_toler_1)**(0.075_dp) &
* (error_toler_0/error_toler_1)**(0.1_dp) * dt
error_toler_0 = error_toler_1
PRINT *, "dt for next step: ", dt, "."

max_pot_vorticity = MAXVAL(REAL(ABS(physical_pot_vorticity_grid), dp))
IF (max_pot_vorticity .GT. pot_vorticity_bound) THEN
  WRITE(*,"(A,F10.3,A,F10.3,A)") "ERROR: Max physical_pot_vorticity_grid = ", &
    max_pot_vorticity, " exceeds pot_vorticity_bound = ", pot_vorticity_bound, "."
  ERROR STOP
END IF

END SUBROUTINE TIME_STEP_SCHEME

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE CALCULATE_ARK_COEFF(erk_coeff, esdirk_coeff)
! Calculates the coefficients for an additive Runge-Kutta method
! described by Kennedy, C. A., et. al in "Additive Runge-Kutta Schemes for
! Convection-Diffusion-Reaction Equations" (July 2001). We use
! ARK4(3)6L[2]SA - ERK for the Jacobian, Ekamn friction, and vertical shear
! terms (which act on slow time scales) and ARK4(3)6L[2]SA - ESDIRK for
! the hyperviscosity term (which acts on fast time scales). Since the latter
! is diagonally implicit, the equation for each stage has been rearranged to
! isolate the new value of potential vorticity at that stage.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

END MODULE
