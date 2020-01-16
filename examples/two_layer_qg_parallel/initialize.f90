!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Initializaes the simultion, including reading the input file and setting
!  the initial conditions for the grid.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE INITIALIZE

IMPLICIT NONE

PRIVATE

! Defines standard integer-, real-precision types.
INCLUDE 'integer_types.h'
INCLUDE 'real_types.h'

INTEGER(qb) :: x_len, y_len, num_timesteps, biharm_order, output_freq
REAL(dp) :: x_ref, y_ref, dx, dy, init_dt, pot_vorticity_bound, &
& deform_wavenumber, rotat_wavenumber, vertical_shear, ekman_friction_coeff, &
& biharm_visc_coeff, init_time, error_toler
COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: physical_pot_vorticity_grid

PUBLIC :: x_len, y_len, num_timesteps, biharm_order, output_freq
PUBLIC :: x_ref, y_ref, dx, dy, init_dt, pot_vorticity_bound, &
& deform_wavenumber, rotat_wavenumber, vertical_shear, ekman_friction_coeff, &
& biharm_visc_coeff, init_time, error_toler
PUBLIC :: physical_pot_vorticity_grid
PUBLIC :: INITIALIZE_PARAMETERS, INITIALIZE_GRID

CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE INITIALIZE_PARAMETERS
!  Reads the namelist containing the parameters for the simulation.
!
!  Variables:
!  - xy_len : Number of grid points in the x-, y- direction.
!  - dx, dy : Grid spacing in the x-, y-directions, for the intial condition.
!  - num_timesteps : Number of timesteps to go through
!  - init_dt : The initial timestep size.
!  - pot_vorticity_bound : Program aborts if the potential vorticity exceeds
!    this bound at any point at any time.
!  - deform_wavenumber : The baroclinic deformation wavenumber corrseponding
!    to the Rossby radius of deformation.
!  - rotat_wavenumber : The wavenumber corresponding to the rotation
!    coefficient, which controls the advection of streamfunctions.
!  - vertical_shear : Large-scale vertical shear, opposite in each direction in
!    background to induce baroclinic instability.
!  - biharm_visc_coeff : The coefficient for hyperviscosity in the system.
!  - x_ref, y_ref : The x-, y-coordinates of the reference point, for the
!    initial condition.
!  - error_toler : Tolerance parameter for adaptive time-stepping.
!  - output_freq : Output frequency of the simulation.
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IMPLICIT NONE

NAMELIST /model/ x_len, y_len, dx, dy, num_timesteps, init_dt, pot_vorticity_bound, &
& deform_wavenumber, rotat_wavenumber, vertical_shear, ekman_friction_coeff, &
& biharm_visc_coeff, biharm_order, x_ref, y_ref, init_time, error_toler, &
& output_freq

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

INTEGER(qb) :: i, j

! The physical potential vorticity grid contains layer 1 (:,:,1) and
! layer 2 (:, :, 2).
ALLOCATE(physical_pot_vorticity_grid(x_len, y_len, 2))
physical_pot_vorticity_grid(:,:,1) = (0.0_dp, 0.0_dp)

! Fill in the initial condition for the physical potential vorticity grid.
DO j = 1, y_len
  DO i = 1, x_len
    physical_pot_vorticity_grid(i, j, 1) = &
      & CMPLX(INITIAL_CONDITION_1(x_ref + (i-1) * dx, y_ref + (j-1) * dy), 0.0_dp, dp)
    physical_pot_vorticity_grid(i, j, 2) = &
      & CMPLX(INITIAL_CONDITION_2(x_ref + (i-1) * dx, y_ref + (j-1) * dy), 0.0_dp, dp)
  END DO
END DO


END SUBROUTINE INITIALIZE_GRID

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REAL(dp) FUNCTION INITIAL_CONDITION_1(x_pos, y_pos) RESULT(output)
!  Function to set the initial condition for layer 1 based on grid point
!  position.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IMPLICIT NONE

REAL(dp), INTENT(in) :: x_pos, y_pos
REAL(dp), PARAMETER :: pi = 4.0_dp * ATAN(1.0_dp)

! f(x, y) = cos(2*pi*x) * sin(2*pi*y)
output = COS(2 * pi * x_pos) * SIN(2 * pi * y_pos)

END FUNCTION INITIAL_CONDITION_1

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REAL(dp) FUNCTION INITIAL_CONDITION_2(x_pos, y_pos) RESULT(output)
!  Function to set the initial condition for layer 2 based on grid point
!  position.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IMPLICIT NONE

REAL(dp), INTENT(in) :: x_pos, y_pos
REAL(dp), PARAMETER :: pi = 4.0_dp * ATAN(1.0_dp)

! f(x, y) = cos(2*pi*x) * sin(2*pi*y)
output = COS(2 * pi * x_pos) * SIN(2 * pi * y_pos)

END FUNCTION INITIAL_CONDITION_2

END MODULE
