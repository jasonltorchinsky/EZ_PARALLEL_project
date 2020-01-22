!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! INITIALIZE. Contains PUBLIC parameters about the temperature grid (temp_grid),
! and subroutines for reading the NAMELIST parameters and filling in the
! initial condition.
!
! VARIABLES: - x_len, y_len: The size (dimension 1, dimension 2) of the grid for
! simulation (INTEGER(qb), PUBLIC).
! - num_timesteps: Number of timesteps to be performed (INTEGER(qb), PUBLIC).
! - biharm_order: Order of the biharmonic operator for the hyperviscosity
! (INTEGER(qb), PUBLIC).
! - output_freq: The output frequency for the simulation, at every output_freq
! timestep, the grid will be written to file (INTEGER(qb), PUBLIC).
! - x_ref, y_ref: The coordinates of the reference point (lower-left corner of
! the grid) for setting the initial condition (REAL(dp), PUBLIC).
! - dx, dy: The grid spacing in the x- and y-directions (REAL(dp), PUBLIC).
! - init_dt: The initial timestep size (REAL(dp), PUBLIC).
! - pot_vort_bound: Program aborts if the potential vorticity (pot_vort)
! exceeds this bound at any timestep (REAL(dp), PUBLIC).
! - deform_wavenum: The baroclinic deformation wavenumber corrseponding
! to the Rossby radius of deformation (REAL(dp), PUBLIC).
! - rotat_wavenum: The wavenumber corresponding to the rotation
! coefficient, which controls the advection of streamfunctions (REAL(dp),
! PUBLIC).
! - vert_shear: Large-scale vertical shear, opposite in each direction in
! background to induce baroclinic instability (REAL(dp), PUBLIC).
! - ekman_fric_coeff: The coefficient for Ekman friction in layer 2
! (REAL(dp), PUBLIC).
! - biharm_visc_coeff : The coefficient for hyperviscosity in the system
! (REAL(dp), PUBLIC).
! - init_time: The starting/ending time for the simulation (REAL(dp), PUBLIC).
! - error_toler: Tolerance parameter for adaptive time-stepping (REAL(dp),
! PUBLIC).
! - phys_pot_vort_grid: The pot_vort grid in physical space that will be evolved
! in time (COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC).
! - INITIALIZE_PARAMETERS: The subroutine which reads the NAMELIST (PUBLIC).
! - INITIALIZE_GRID: The subroutine which allocates memory for temp_grid and
! fills in the initial condition (PUBLIC).
!
! Written By: Jason Turner
! Last Updated: January 21, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE INITIALIZE

IMPLICIT NONE

PRIVATE

! Defines standard integer-, real-precision types.
INCLUDE 'integer_types.h'
INCLUDE 'real_types.h'

INTEGER(qb), PUBLIC :: x_len, &
& y_len, &
& num_timesteps, &
& biharm_order, &
& output_freq
REAL(dp), PUBLIC :: x_ref, &
& y_ref, &
& dx, &
& dy, &
& init_dt, &
& pot_vort_bound, &
& deform_wavenum, &
& rotat_wavenum, &
& vert_shear, &
& ekman_fric_coeff, &
& biharm_visc_coeff, &
& init_time, &
& error_toler
COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: phys_pot_vort_grid
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
& num_timesteps, &
& init_dt, &
& pot_vort_bound, &
& deform_wavenum, &
& rotat_wavenum, &
& vert_shear, &
& ekman_fric_coeff, &
& biharm_visc_coeff, &
& biharm_order, &
& x_ref, &
& y_ref, &
& init_time, &
& error_toler, &
& output_freq

OPEN(1000, file = "NAMELIST")
READ(1000, nml = model)
CLOSE(1000)

END SUBROUTINE INITIALIZE_PARAMETERS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Allocate memory and fills in the initial condition for temp_grid.
!
! STRUCTURE: 1) Allocates memory to phys_pot_vort_grid.
! 2) Sets the initial condition on the phys_pot_vort_grid.
!
! VARIABLES: - i, j: Counting indices for DO loops (INTEGER(qb)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE INITIALIZE_GRID
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IMPLICIT NONE

INTEGER(qb) :: i, &
& j

! The physical potential vorticity grid contains layer 1 (:,:,1) and
! layer 2 (:, :, 2).
ALLOCATE(phys_pot_vort_grid(x_len, y_len, 2))
phys_pot_vort_grid(:,:,:) = (0.0_dp, 0.0_dp)

! Fill in the initial condition for the physical potential vorticity grid.
DO j = 1, y_len
  DO i = 1, x_len
    phys_pot_vort_grid(i, j, 1) = &
    & CMPLX(INITIAL_CONDITION_1(x_ref + (i-1) * dx, y_ref + (j-1) * dy), &
      & 0.0_dp, dp)
    phys_pot_vort_grid(i, j, 2) = &
    & CMPLX(INITIAL_CONDITION_2(x_ref + (i-1) * dx, y_ref + (j-1) * dy), &
      & 0.0_dp, dp)
  END DO
END DO


END SUBROUTINE INITIALIZE_GRID
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The initial potential vorticity distribution for layer 1.
!
! STRUCTURE: 1) Calculates pi.
! 2) Calculates the value of the distribution at this input position.
!
! VARIABLES: - pi_dp: Constant for pi (3.1415...) (REAL(dp)).
! - x_pos, y_pos: The spatial coordinates of the input (REAL(dp),
! INTENT(IN)).
! - output: The value of the distribution at x_pos, y_pos (REAL(dp)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REAL(dp) FUNCTION INITIAL_CONDITION_1(x_pos, y_pos) RESULT(output)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IMPLICIT NONE

REAL(dp) :: pi_dp
REAL(dp), INTENT(in) :: x_pos, &
& y_pos

pi_dp = 4.0_dp * ATAN(1.0_dp)

! f(x, y) = sin(x/x_len * 2 * pi)
output = SIN(x_pos/x_len * 2.0_dp * pi_dp)

END FUNCTION INITIAL_CONDITION_1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The initial potential vorticity distribution for layer 2.
!
! STRUCTURE: 1) Calculates pi.
! 2) Calculates the value of the distribution at this input position.
!
! VARIABLES: - pi_dp: Constant for pi (3.1415...) (REAL(dp)).
! - x_pos, y_pos: The spatial coordinates of the input (REAL(dp),
! INTENT(IN)).
! - output: The value of the distribution at x_pos, y_pos (REAL(dp)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REAL(dp) FUNCTION INITIAL_CONDITION_2(x_pos, y_pos) RESULT(output)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IMPLICIT NONE

REAL(dp) :: pi_dp
REAL(dp), INTENT(in) :: x_pos, &
& y_pos

pi_dp = 4.0_dp * ATAN(1.0_dp)

! f(x, y) = COS(y/y_len * 2 * pi)
output = COS(y_pos/y_len * 2.0_dp * pi_dp)

END FUNCTION INITIAL_CONDITION_2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END MODULE
