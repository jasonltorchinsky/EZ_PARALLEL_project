!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! INITIALIZE. Contains PUBLIC parameters about the temperature grid (temp_grid),
! and subroutines for reading the NAMELIST parameters and filling in the
! initial condition.
!
! VARIABLES: - x_len, y_len: The size (dimension 1, dimension 2) of the grid for
! simulation (INTEGER(qb), PUBLIC).
! - output_freq: The output frequency for the simulation, at every output_freq
! timestep, the grid will be written to file (INTEGER(qb), PUBLIC).
! - x_ref, y_ref: The coordinates of the reference point (lower-left corner of
! the grid) for setting the initial condition (REAL(dp), PUBLIC).
! - dx, dy: The grid spacing in the x- and y-directions (REAL(dp), PUBLIC).
! - dt: The timestep size (REAL(dp), PUBLIC).
! - x_diffus, y_diffus: The heat diffusivity in the x- and y-directions
! (REAL(dp), PUBLIC).
! - deter_force: Magnitude of the deterministic forcing (REAL(dp), PUBLIC).
! - damp_coeff: Magnitude of the dampening coefficient (REAL(dp), PUBLIC).
! - stoch_mag: Magnitude of the stochastic noise (REAL(dp), PUBLIC).
! - relax_trgt: Relaxation target (REAL(dp), PUBLIC).
! - dir_bdy_val: Value used for the Dirichlet boundary condition.
! - init_time, fin_time: The starting/(max) ending time for the simulation
! (REAL(dp), PUBLIC). Note, the simulation will never go further than one time
! step past the ending time.
! - temp_grid: The temperature grid that will be evolved in time (REAL(dp),
! DIMENSION(:,:,:), ALLOCATABLE, PUBLIC). Note, this grid has a third dimension
! which is used to hold two time steps.
! - INITIALIZE_PARAMETERS: The subroutine which reads the NAMELIST (PUBLIC).
! - INITIALIZE_GRID: The subroutine which allocates memory for temp_grid and
! fills in the initial condition.
!
! Written By: Jason Turner
! Last Updated: January 17, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE INITIALIZE

IMPLICIT NONE

PRIVATE

! Defines standard integer-, real-precision types.
INCLUDE 'integer_types.h'
INCLUDE 'real_types.h'

INTEGER(qb), PUBLIC :: x_len, &
& y_len, &
& output_freq
REAL(dp), PUBLIC :: x_ref, &
& y_ref, &
& dx, &
& dy, &
& dt, &
& x_diffus, &
& y_diffus, &
& deter_force, &
& damp_coeff, &
& stoch_mag, &
& relax_trgt, &
& dir_bdy_val, &
& init_time, &
& fin_time

REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: temp_grid
PUBLIC :: INITIALIZE_PARAMETERS, &
& INITIALIZE_GRID

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Reads the NAMELIST.
!
! STRUCTURE: 1) Reads the NAMELIST.
!
! VARIABLES: N/A.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE INITIALIZE_PARAMETERS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IMPLICIT NONE

NAMELIST /model/ x_len, &
& y_len, &
& dx, &
& dy, &
& dt, &
& x_diffus, &
& y_diffus, &
& deter_force, &
& damp_coeff, &
& stoch_mag, &
& relax_trgt, &
& dir_bdy_val, &
& x_ref, &
& y_ref, &
& init_time, &
& fin_time, &
& output_freq

OPEN(1000, file = "NAMELIST")
READ(1000, nml = model)
CLOSE(1000)

END SUBROUTINE INITIALIZE_PARAMETERS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Allocate memory and fills in the initial condition for temp_grid.
!
! STRUCTURE: 1) Allocates memory to temp_grid.
! 2) Sets the initial condition on the temp_grid.
! 3) Sets a homogeneous Dirichlet boundary conition.
!
! VARIABLES: - i, j: Counting indices for DO loops (INTEGER(qb)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE INITIALIZE_GRID
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IMPLICIT NONE

INTEGER(qb) :: i, &
& j

! temp_grid contains two time steps (third index).
ALLOCATE(temp_grid(x_len, y_len, 2))
temp_grid(:,:,1) = 0.0_dp

! Fill in the initial condition, based on the INITIAL_CONDITION function defined
! below.
DO j = 1, y_len
  DO i = 1, x_len
    temp_grid(i, j, 1) = &
      & INITIAL_CONDITION(x_ref + (i-1) * dx, y_ref + (j-1) * dy)
  END DO
END DO

! Set boundary value to a fixed value.
! Faster to assign value to array which is contiguous in memory all at once.
temp_grid(:, 1, 1) = dir_bdy_val
temp_grid(:, y_len, 1) = dir_bdy_val
! Faster to assign value to array which is not contiguous in memory with DO.
DO j = 1, y_len
  temp_grid(1, j, 1) = dir_bdy_val
  temp_grid(x_len, j, 1) = dir_bdy_val
END DO

temp_grid(:,:,2) = temp_grid(:,:,1)

END SUBROUTINE INITIALIZE_GRID
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The initial temperature distribution.
!
! STRUCTURE: 1) Calculates pi.
! 2) Calculates the value of the distribution at this input position.
!
! VARIABLES: - pi_dp: Constant for pi (3.1415...) (REAL(dp)).
! - x_pos, y_pos: The spatial coordinates of the input (REAL(dp),
! INTENT(IN)).
! - output: The value of the distribution at x_pos, y_pos (REAL(dp)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REAL(dp) FUNCTION INITIAL_CONDITION(x_pos, y_pos) RESULT(output)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IMPLICIT NONE

REAL(dp) :: pi_dp
REAL(dp), INTENT(in) :: x_pos, &
& y_pos

pi_dp = 4.0_dp * ATAN(1.0_dp)

! Random Initial Condition.
CALL RANDOM_NUMBER(output)
output = relax_trgt * (output + 0.5)

END FUNCTION INITIAL_CONDITION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END MODULE
