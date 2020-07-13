!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : STOCH_HEAT_SERIAL TIME_STEPPER
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
!
!> @brief Module containing the time-stepping subroutine needed for the example
!! serial stochastic heat equation code.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE TIME_STEPPER

IMPLICIT NONE

PRIVATE

! Defines standard integer-, real-precision types.
INCLUDE 'integer_types.h'
INCLUDE 'real_types.h'

PUBLIC :: TIME_STEP

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Driver for the time-stepping. First ensures that the CFL condition is met
  !! and that the outputFreq is not too large.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE TIME_STEP

    USE INITIALIZE
    USE OUTPUT

    IMPLICIT NONE

    INTEGER(qb) :: step !< Current time-step number.
    INTEGER(qb) :: stepParity !< Parity of the time-step, to determine if the
    !! step is stored in colWtrVpr(:,:,0) or colWtrVpr(:,:,1).
    REAL(dp) :: cflCondition !< CFL condition parameter.

    ! Stability constraint derived from our time step method, which is forward
    ! Euler in time, and a 1st order centered difference for each spatial
    ! derivative (using Von Neumann analysis). We check it here even though this
    ! code uses a different differential equation just for safety.
    cflCondition = 1.0_dp/(2.0_dp * (diffCoeff/(dx**2.0_dp) &
         + diffCoeff/(dy**2.0_dp)))

    IF ((dt .GT. cflCondition)) THEN
       PRINT *, "WARNING: Time step size ", dt, " exceeds the CFL parameter ", &
            cflCondition, ". Simulation may be unstable."
    END IF

    stepParity = 0_qb
    step = 0_qb

    IF (numSteps .LT. outputFreq) THEN
       ERROR STOP 'Output frequency too large. Please reduce.'
    END IF

    ! Output the initial condition.
    IF (outputFreq .NE. 0) THEN
       CALL WRITE_OUTPUT(step, stepParity)
    END IF

    ! Step simulation forward in time.
    DO WHILE (step .LT. numSteps)
       stepParity = 2_qb - (stepParity + 1_qb)
       CALL TIME_STEP_SCHEME(stepParity)
       step = step + 1_qb

       ! Write output at desired frequency.
       IF (outputFreq .NE. 0) THEN
          IF (MOD(step, outputFreq) .EQ. 0) THEN
             CALL WRITE_OUTPUT(step, stepParity)
          END IF
       END IF
    END DO

  END SUBROUTINE TIME_STEP

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Executes the numerical scheme for the time-step, using a Forward Euler
  !! scheme in time and a first-order center difference in space. We also
  !! include a deterministic forcing, a relaxation term, and white noise
  !! (see Equation (2) of Hottovy, S., Stechmann, S. (2015)).
  !
  !> @param[in] stepParity Parity of the time-step.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE TIME_STEP_SCHEME(stepParity)

    USE INITIALIZE

    IMPLICIT NONE

    INTEGER(qb), INTENT(IN) :: stepParity
    INTEGER(qb) :: prevStepParity !< Parity of the previous time-step.
    REAL(dp), PARAMETER :: pi_dp = 4.0_dp * ATAN(1.0_dp) !< Double precision pi.
    REAL(dp) :: whiteNoise1(1:numPts-2,1:numPts-2), whiteNoise2(1:numPts-2,1:numPts-2)
    !< White noise for the stochasticity of the time-step.

    ! Get a random number for the Gaussian white noise.
    CALL RANDOM_NUMBER(whiteNoise1)
    CALL RANDOM_NUMBER(whiteNoise2)
    whiteNoise1 = stochMag * SQRT(dt)* SQRT(-2.0_dp * LOG(whiteNoise1)) * COS(2.0_dp * pi_dp * whiteNoise2) 
    !< Gaussian white noise using Box-Muller transform.

    ! Update the step parity.
    prevStepParity = 2_qb - (stepParity + 1_qb)

    ! Step forward in time.
    colWtrVpr(1:numPts-2, 1:numPts-2, stepParity) = &
         colWtrVpr(1:numPts-2, 1:numPts-2, prevStepParity) &
         + dt * deterForce &
         + dt * interCoeff * (colWtrVpr(0:numPts-3, 1:numPts-2, prevStepParity) &
         + colWtrVpr(2:numPts-1, 1:numPts-2, prevStepParity) &
         + colWtrVpr(1:numPts-2, 0:numPts-3, prevStepParity) &
         + colWtrVpr(1:numPts-2, 2:numPts-1, prevStepParity) &
         - 4.0_dp * colWtrVpr(1:numPts-2, 1:numPts-2, prevStepParity)) &
         - (dt / dampCoeff) * (colWtrVpr(1:numPts-2, 1:numPts-2, prevStepParity) &
         - relaxTrgt) + whiteNoise1

  END SUBROUTINE TIME_STEP_SCHEME

END MODULE TIME_STEPPER
