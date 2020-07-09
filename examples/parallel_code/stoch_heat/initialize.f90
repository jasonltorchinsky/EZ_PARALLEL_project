!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : STOCH_HEAT_PARALLEL INITIALIZE
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
!
!> @brief Module containing the initialization subroutines needed for the
!! example parallel stochastic heat equation code.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE INITIALIZE

  USE MPI !< ADDED TO PARALLEL.
  USE EZ_PARALLEL_STRUCTS !< ADDED TO PARALLEL.
  USE EZ_PARALLEL !< ADDED TO PARALLEL.

  IMPLICIT NONE

  PRIVATE

  ! Defines standard integer-, real-precision types.
  INCLUDE 'integer_types.h'
  INCLUDE 'real_types.h'

  INTEGER(qb), PUBLIC :: numPts !< Number of lattice points in one direction (N)
  INTEGER(qb), PUBLIC :: numPtsY !< ADDED TO PARALLEL. Number of points in the
  !! y-direction (along the rows).
  REAL(dp), PUBLIC :: dx !< Size of one grid cell in the x-direction (Delta x)
  REAL(dp), PUBLIC :: dy !< Size of one grid cell in the y-direction (Delta y)
  REAL(dp), PUBLIC :: deterForce !< Magnitude of the deterministic forcing (F)
  REAL(dp), PUBLIC :: dampCoeff !< Dampening coefficient (tau)
  REAL(dp), PUBLIC :: diffCoeff !< Diffusion coefficient (b)
  REAL(dp), PUBLIC :: stochMag !< Magnitude of stochastic forcing (D_*)
  REAL(dp), PUBLIC :: relaxTrgt !< Relaxation target (q^*)
  REAL(dp), PUBLIC :: dt !< Time-step size.
  INTEGER(dp), PUBLIC :: numSteps !< Number of time-steps to execute.
  INTEGER(qb), PUBLIC :: outputFreq !< Output frequency for the simulation.
  !! 1 = output a file each time-step. 0 = no output.

  REAL(dp), PUBLIC :: sysLen !< System size in one direction (L)
  REAL(dp), PUBLIC :: interCoeff !< Interation coefficient (b_0)
  REAL(dp), ALLOCATABLE, PUBLIC :: colWtrVpr(:,:,:) !< The integrated column water
  !! vapor. We store both the current and the previous time-step for the forward
  !! Euler time-integration method.

  TYPE(SCHEME), PUBLIC :: sch !< <tt>SCHEME</tt> For use of EZ_PARALLEL. ADDED
  !! TO PARALLEL.
  
  PUBLIC :: INITIALIZE_PARAMETERS
  PUBLIC :: INITIALIZE_GRID
  
CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Reads the NAMELIST for input parameters, and calculate some other constant
  !! parameters.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE INITIALIZE_PARAMETERS

    IMPLICIT NONE

    NAMELIST /model/ numPts, dx, dy, deterForce, dampCoeff, diffCoeff, &
         stochMag, relaxTrgt, dt, numSteps, outputFreq

    OPEN(1000, file = "NAMELIST")
    READ(1000, nml = model)
    CLOSE(1000)

    IF (dx .NE. dy) THEN
       PRINT *, "Lattice spacing is not the same in both directions, ", &
            " defaulting to dx for dy."
    END IF
    
    sysLen = numPts * dx
    interCoeff = diffCoeff / (dx ** 2.0_dp)
    numPtsY = numPts

  END SUBROUTINE INITIALIZE_PARAMETERS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Allocates memory to colWtrVpr  and fills in the initial condition.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE INITIALIZE_GRID

    USE MPI !< ADDED TO PARALLEL.
    USE EZ_PARALLEL_STRUCTS !< ADDED TO PARALLEL.
    USE EZ_PARALLEL !< ADDED TO PARALLEL.

    IMPLICIT NONE

    
    COMPLEX(dp), ALLOCATABLE :: initCond(:,:)
    REAL(dp) :: yRef, yRefInit !< ADDED TO PARALLEL. Y-position of reference
    !! point on grid and initial grid.
    INTEGER(qb) :: i, j !< Counters for DO loops.
    INTEGER(qb) :: numPtsYInit !< Number of points in the y-direction along the
    !! initial grid.
    TYPE(SCHEME) :: schInit !< ADDED TO PARALLEL. Scheme for initial condition grid.

    yRef = 0.0_dp
    yRefInit = 0.0_dp
    numPtsYInit = numPtsY
    
    CALL CREATE_SCHEME(numPts, numPtsYInit, 1.0_dp, yRef, MPI_COMM_WORLD, &
         MPI_DOUBLE_COMPLEX, 1_qb, schInit) !< ADDED TO PARALLEL.
    ! Since we initialize the grid in Fourier space, we set the grid spacing
    ! in the y-direction to be 1 and the reference position in the y-direction
    ! to be zero.
    CALL CREATE_SCHEME_FFT(schInit) !< ADDED TO PARALLEL.
    CALL CREATE_SCHEME(numPts, numPtsY, 1.0_dp, yRef, MPI_COMM_WORLD, &
         MPI_DOUBLE_PRECISION, 1_qb, sch) !< ADDED TO PARALLEL.
    
    ! colWtrVpr contains two time steps (third index).
    ALLOCATE(colWtrVpr(0:numPts-1,0:numPtsY-1,0:1)) !< CHANGED FOR PARALLEL.
    colWtrVpr(:,:,0) = 0.0_dp

    ! We initialze a complex grid for the inverse FFT.
    ALLOCATE(initCond(0:numPts-1,0:numPtsY-1)) !< CHANGED FOR PARALLEL.
    initCond = (0.0_dp, 0.0_dp)
    
    ! Fill in the initial condition, which is uncorrelated Gaussian noise in
    ! Fourier space.
    DO j = 0, numPtsY-1 !< CHANGED FOR PARALLEL.
       DO i = 0, numPts-1 !< CHANGED FOR PARALLEL.
          initCond(i,j) = INITIAL_CONDITION(i,INT(yRef,qb) + j) !< CHANGED FOR PARALLEL.
       END DO
    END DO

    ! Inverse FFT the initial condition.
    !CALL CFFT2DB(numPts-2_qb, numPts-2_qb, initCond(1:numPts-2,1:numPts-2))
    CALL EXECUTE_SCHEME_IFFT(initCond, FFT_2D, schInit) !< ADDED TO PARALLEL.
    colWtrVpr(:,:,0) = REAL(initCond, dp)
    colWtrVpr = (numPts-1_qb) * (numPts-1_qb) * colWtrVpr ! Inverse FFT of
    ! EZ_PARALLEL is normalized, whereas the one in the paper is not.
    
    ! Set boundary value to the relax target. We use Dirichlet boundary
    ! conditions for compatibility with EZ_PARALLEL.
    ! Faster to assign value to array which is contiguous in memory all at once.
    colWtrVpr(:,0,0) = relaxTrgt
    colWtrVpr(:,numPtsY-1,0) = relaxTrgt !< CHANGED FOR PARALLEL.
    ! Faster to assign value to array which is not contiguous in memory with DO.
    DO j = 0, numPtsY-1 !< CHANGED FOR PARALLEL.
       colWtrVpr(0,j,0) = relaxTrgt
       colWtrVpr(numPts-1,j,0) = relaxTrgt
    END DO

    CALL DESTROY_SCHEME(schInit)

    CALL SHARE_SUBGRID_BDRY(colWtrVpr(:,:,0), sch) !< ADDED TO PARALLEL.

    colWtrVpr(:,:,1) = colWtrVpr(:,:,0)

  END SUBROUTINE INITIALIZE_GRID

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> The initial condition function. We use Gaussian noise in Fourier space.
  !
  !> @param[in] xWvNum The wavenumber in the spectral x-direction.
  !> @param[in] yWvNum The wavenumber in the spectral y-direction.
  !> @param[out] output The value of the initial condition at (xWvNum, yWvNum).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION INITIAL_CONDITION(xWvNum, yWvNum) RESULT(output)

    IMPLICIT NONE

    INTEGER(qb), INTENT(IN) :: xWvNum
    INTEGER(qb), INTENT(IN) :: yWvNum
    REAL(dp) :: mean !< Mean of the Gaussian noise at the wavenumber.
    REAL(dp) :: var !< Variance of the Gaussian noise at the wavenumber.
    REAL(dp) :: pi_dp !< Double-precision pi.
    REAL(dp) :: temp1_dp, temp2_dp !< Temporary double-precision value.
    REAL(dp) :: rand1, rand2 !< Random double precision in [0,1) for
    !! generating random Gaussian variable via Box-Mueller transform.

    ! Caluclate mean and variance of Gaussian at wavenumber.
    IF ((xWvNum .EQ. 0_qb) .AND. (yWvNum .EQ. 0_qb)) THEN
       mean = 0.0_dp * deterForce * dampCoeff + 1.5_dp *  relaxTrgt
    ELSE
       mean = 0.0_dp
    END IF

    pi_dp = 4.0_dp * ATAN(1.0_dp)
    temp1_dp = interCoeff * (4.0_dp + (1.0_dp / (interCoeff * dampCoeff)) &
         - 2.0_dp * COS(2.0_dp * pi_dp * REAL(xWvNum, dp) * dx / sysLen) &
         - 2.0_dp * COS(2.0_dp * pi_dp * REAL(yWvNum, dp) * dy / sysLen))
    var = 1.0e-5 * stochMag * stochMag / (2.0_dp * temp1_dp)
    
    ! Random complex initial condition.
    CALL RANDOM_NUMBER(rand1)
    CALL RANDOM_NUMBER(rand2)
    temp1_dp = SQRT(-2.0_dp * LOG(rand1)) * COS(2.0_dp * pi_dp * rand2)
    temp2_dp = SQRT(-2.0_dp * LOG(rand1)) * SIN(2.0_dp * pi_dp * rand2)
    temp1_dp = mean + temp1_dp * SQRT(var)
    temp2_dp = mean + temp2_dp * SQRT(var)
    
    output = 1.0_dp/SQRT(2.0_dp) * CMPLX(temp1_dp, temp2_dp, dp)

  END FUNCTION INITIAL_CONDITION

END MODULE INITIALIZE
