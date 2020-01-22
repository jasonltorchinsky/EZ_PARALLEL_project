!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! JACOBIAN_EKMAN_SHEAR_SOLVE. Calculates the Jacobian term, the Ekman friction
! term, and the vertical shear term in the QG equation for use in progressing
! the simulation forward in time.
!
! VARIABLES: - spec_laplacian: The numerical Laplacian operator in Fourier
! space (COMPLEX(dp), DIMENSION(x_len_sub, y_len_sub)).
! - spec_inv_barotropic, spec_inv_baroclinic: The numerical inverse
! barotropic, baroclinic operators in Fourier space (see two-layer QG equations
! in Di, Q., & Majda, A. (2015)) (COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE).
! - spec_x_deriv, spec_y_deriv: The array that will store the matrix used in
! calculating the derivative of the numerical solution along the x-, y-direction
! of the grid (COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE).
! -scaled_spec_x_deriv, scaled_spec_y_deriv: The array that will store the
! matrix used in calculating the derivative of the numerical solution along the
! x-, y-direction of the scaled grid (the scaled grid is used to de-alias the
! Fourier transform) (COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE).
! -SPECTRAL_X_DERIVATIVE, SPECTRAL_Y_DERIVATIVE: Produces the matrix used in
! calculating the derivative of the numerical solution along the x-,
! y-directions grid (PUBLIC).
! - JACOBIAN_EKMAN_SHEAR_SOLVE: Calculates the Jacobian term, the Ekman friction
! term, and the vertical shear term in the QG equation for use in progressing
! the simulation forward in time (PUBLIC).
!
!
! NOTES: - When calculating the x-, y-wavenumbers, we use the convention
! [0, ..., x_len_sub/2, -x_len_sub/1 + 1, ... -1], with FLOOR(x_len_sub/2), as
! is assumed in Fortran integer division.
!
! Written By: Jason Turner
! Last Updated: January 21, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE JACOBIAN_EKMAN_SHEAR_SOLVE

IMPLICIT NONE

PRIVATE

! Defines standard integer-, real-precision types.
INCLUDE 'integer_types.h'
INCLUDE 'real_types.h'


COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: spec_laplacian, &
& spec_inv_barotropic, &
& spec_inv_baroclinic
COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: spec_x_deriv, &
& spec_y_deriv
COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: scaled_spec_x_deriv, &
& scaled_spec_y_deriv
LOGICAL :: ran = .FALSE.

PUBLIC :: SPECTRAL_X_DERIVATIVE, &
& SPECTRAL_Y_DERIVATIVE, &
& JACOBIAN_EKMAN_SHEAR

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Calculates the Jacobian term, the Ekman friction term, and the vertical shear
! term in the QG equation for use in progressing the simulation forward in time
!
! STRUCTURE: 1) If the subroutine hasn't been called before, calculate all of
! the differential operators needed.
! 2) Calculate the contributions from the mean shear and rotation deformation
! in layers 1 and 2, and the Ekman friction in layer 2.
! 3) Rescale the grid to de-alias it and calculate the derivatives in Fourier
! space.
! 4) Transform into physical space to perform the multiplications to calculate
! the Jacobian, and transform back into frequency space to return the final
! result.
!
! VARIABLES: - x_len_sub, y_len_sub: The size of the grid for simulation
! (INTEGER(qb)).
! - scaled_x_len, scaled_y_len: The size of the scaled-up grid used for
! de-aliasing the Jacobian (INTEGER(qb)).
! - i, j: Counting index used in DO loops (INTEGER(qb)).
! - dt: The current timestep size, adapted throughout simulation (REAL(dp)).
! - deform_wavenum_sub: The baroclinic deformation wavenumber corrseponding
! to the Rossby radius of deformation (REAL(dp)).
! - rotat_wavenum_sub: The wavenumber corresponding to the rotation
! coefficient, which controls the advection of streamfunctions (REAL(dp)).
! - vert_shear_sub: Large-scale vertical shear, opposite in each direction in
! background to induce baroclinic instability (REAL(dp)).
! - ekman_fric_coeff_sub: The coefficient for Ekman friction in layer 2
! (REAL(dp)).
! - infinity (ADDED TO PARALLEL): A very large number to check if 1.0/0.0 has
! ocurred (REAL(dp)).
! - barotropic_pot_vort, baroclinic_pot_vort: The brotropic, baroclinic modes of
! potential vorticity (COMPLEX(dp), DIMENSION(x_len_sub, y_len_sub)).
! - barotropic_pot_vort_strmfunc, baroclinic_pot_vort_strmfunc: The
! streamfunction for barotropic, baroclinic modes of potential vorticity
! (COMPLEX(dp), DIMENSION(x_len_sub, y_len_sub)).
! - freq_pot_vort_grid: The pot vort grid in Fourier space that will be evolved
! in time (COMPLEX(dp), DIMENSION(x_len_sub, y_len_sub, 2)).
! - jacobian_ekman_shear_grid: The change in potential vorticity at a timestep
! due the the Jacobian, Ekman friction, and vertical shear terms of the QG
! equations (COMPLEX(dp), DIMENSION(x_len_sub, y_len_sub, 2)).
! - freq_strmfunc: The streamfunction of the flow in Fourier space (COMPLEX(dp),
! DIMENSION(x_len_sub, y_len_sub, 2)).
! - freq_jacobian: The Jacobian term of the QG equations in Fourier space
! (COMPLEX(dp), DIMENSION(x_len_sub, y_len_sub, 2)).
! - scaled_freq_pot_vort_grid: The scaled pot vort grid in Fourier space that
! is used for de-aliasing the Fourier transforms (COMPLEX(dp), DIMENSION(:,:,:),
! ALLOCATABLE).
! - scaled_freq_strmfunc: The scaled streamfunction of the flow in Fourier space
! that is used for de-aliasing the Fourier transforms (COMPLEX(dp),
! DIMENSION(:,:,:), ALLOCATABLE).
! - scaled_strmfunc_x_deriv, scaled_strmfunc_y_deriv: The scaled x-,
! y-derivatives of the streamfunction in Fourier space, used for de-aliasing the
! Fourier transforms (COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE).
! - scaled_pot_vort_x_deriv, scaled_pot_vort_y_deriv: The scaled x-,
! y-derivatives of the pot vort in Fourier space, used for de-aliasing the
! Fourier transforms (COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE).
! - scaled_freq_jacobian: The scaled Jacobian term of the QG equations in
! Fourier space used for de-aliasing the Fourier transforms (COMPLEX(dp),
! DIMENSION(:,:,:), ALLOCATABLE).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE JACOBIAN_EKMAN_SHEAR(freq_pot_vort_grid, x_len_sub, y_len_sub, dt, &
& deform_wavenum_sub, rotat_wavenum_sub, vert_shear_sub, ekman_fric_coeff_sub, &
& jacobian_ekman_shear_grid)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE INITIALIZE

IMPLICIT NONE

INTEGER(qb) :: x_len_sub, &
& y_len_sub, &
& scaled_x_len, &
& scaled_y_len, &
& i, &
& j ! ADDED TO PARALLEL
REAL(dp) :: dt, &
& deform_wavenum_sub, &
& rotat_wavenum_sub, &
& vert_shear_sub, &
& ekman_fric_coeff_sub, &
& infinity ! ADDED TO PARALLEL
COMPLEX(dp), DIMENSION(x_len_sub, y_len_sub) :: barotropic_pot_vort, &
& baroclinic_pot_vort, &
& barotropic_pot_vort_strmfunc, &
& baroclinic_pot_vort_strmfunc
COMPLEX(dp), DIMENSION(x_len_sub, y_len_sub, 2) :: freq_pot_vort_grid, &
& jacobian_ekman_shear_grid, &
& freq_strmfunc, &
& freq_jacobian
COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: scaled_freq_pot_vort_grid, &
& scaled_freq_strmfunc, &
& scaled_strmfunc_x_deriv, &
& scaled_strmfunc_y_deriv, &
& scaled_pot_vort_x_deriv, &
& scaled_pot_vort_y_deriv, &
& scaled_freq_jacobian

! If it hasn't been called before, calculate all differential operators needed.
IF (.NOT. ran) THEN
  ! Get the spectral Laplacian operator.
  ALLOCATE(spec_x_deriv(x_len_sub, y_len_sub, 1))
  ALLOCATE(spec_y_deriv(x_len_sub, y_len_sub, 1))
  !CALL SPECTRAL_X_DERIVATIVE(x_len_sub, y_len_sub, spec_x_deriv, 2_qb)
    ! CHANGED FOR PARALLEL.
  !CALL SPECTRAL_Y_DERIVATIVE(x_len_sub, y_len_sub, spec_y_deriv, 2_qb)
    ! CHANGED FOR PARALLEL.
  CALL SPECTRAL_DIM1_DERIVATIVE_EZP(x_len_sub, y_len_sub, 0, spec_x_deriv, &
    & 2_qb) ! ADDED TO SERIAL
  CALL SPECTRAL_DIM2_DERIVATIVE_EZP(x_len_sub, y_len_sub, 0, spec_y_deriv, &
    & 2_qb) ! ADDED TO SERIAL
  ALLOCATE(spec_laplacian(x_len_sub, y_len_sub))
  spec_laplacian = spec_x_deriv(:,:,1) + spec_y_deriv(:,:,1)
  DEALLOCATE(spec_x_deriv)
  DEALLOCATE(spec_y_deriv)

  ! Get the inverse spectral barotropic and inverse baroclinic operators.
  ALLOCATE(spec_inv_barotropic(x_len_sub, y_len_sub))
  ALLOCATE(spec_inv_baroclinic(x_len_sub, y_len_sub))
  spec_inv_barotropic = 1.0_dp/spec_laplacian
  !spec_inv_barotropic(1,1) = 0.0_dp ! CHANGED FOR PARALLEL
  spec_inv_baroclinic = 1.0_dp/(spec_laplacian - deform_wavenum_sub**(2.0_dp))
  !spec_inv_baroclinic(1,1) = 0.0_dp ! CHANGED FOR PARALLEL
  !~~ ADDED BLOCK TO PARALLEL
  infinity = HUGE(infinity)
  DO j = 1_qb, y_len_sub
    DO i = 1_qb, x_len_sub
      IF (ABS(spec_inv_barotropic(i,j)) .GT. infinity) THEN
        spec_inv_barotropic(i,j) = 0.0_dp
      END IF
      IF (ABS(spec_inv_baroclinic(i,j)) .GT. infinity) THEN
        spec_inv_baroclinic(i,j) = 0.0_dp
      END IF
    END DO
  END DO
  !~~

  ! Get the first-order differential operators.
  ALLOCATE(spec_x_deriv(x_len_sub, y_len_sub, 2))
  ALLOCATE(spec_y_deriv(x_len_sub, y_len_sub, 2))
  !CALL SPECTRAL_X_DERIVATIVE(x_len_sub, y_len_sub, spec_x_deriv(:,:,1), 1_qb)
    ! CHANGED FOR PARALLEL
  CALL SPECTRAL_DIM1_DERIVATIVE_EZP(x_len_sub, y_len_sub, 0_qb, spec_x_deriv, &
    & 1_qb) ! ADDED TO SERIAL
  spec_x_deriv(:,:,2) = spec_x_deriv(:,:,1)
  !CALL SPECTRAL_Y_DERIVATIVE(x_len_sub, y_len_sub, spec_y_deriv(:,:,1), 1_qb)
    ! CHANGED FOR PARALLEL
  CALL SPECTRAL_DIM2_DERIVATIVE_EZP(x_len_sub, y_len_sub, 0_qb, spec_y_deriv, &
    & 1_qb) ! ADDED TO SERIAL
  spec_y_deriv(:,:,2) = spec_y_deriv(:,:,1)

  ALLOCATE(scaled_spec_x_deriv(3_qb * x_len_sub/2_qb, 3_qb * y_len_sub/2_qb, 2))
  ALLOCATE(scaled_spec_y_deriv(3_qb * x_len_sub/2_qb, 3_qb * y_len_sub/2_qb, 2))
  !CALL SPECTRAL_X_DERIVATIVE(3_qb * x_len_sub/2_qb, 3_qb * y_len_sub/2_qb, &
  !  & scaled_spec_x_deriv(:,:,1), 1_qb) ! CHANGED FOR PARALLEL
  CALL SPECTRAL_DIM1_DERIVATIVE_EZP(3_qb * x_len_sub/2_qb, &
    & 3_qb * y_len_sub/2_qb, 0_qb, scaled_spec_x_deriv, 1_qb) ! ADDED TO SERIAL
  scaled_spec_x_deriv(:,:,2) = scaled_spec_x_deriv(:,:,1)
  !CALL SPECTRAL_Y_DERIVATIVE(3_qb * x_len_sub/2_qb, 3_qb * y_len_sub/2_qb, &
  !  & scaled_spec_y_deriv(:,:,1), 1_qb) ! CHANGED FOR PARALLEL
  CALL SPECTRAL_DIM2_DERIVATIVE_EZP(3_qb * x_len_sub/2_qb, &
    & 3_qb * y_len_sub/2_qb, 0_qb, scaled_spec_y_deriv, 1_qb) ! ADDED TO SERIAL
  scaled_spec_y_deriv(:,:,2) = scaled_spec_y_deriv(:,:,1)

  ran = .TRUE.
END IF

! Calculate the spectral baroclinic and barotropic potential vorticities.
barotropic_pot_vort = 0.5_dp * (freq_pot_vort_grid(:,:,1) &
& + freq_pot_vort_grid(:,:,2))
baroclinic_pot_vort = 0.5_dp * (freq_pot_vort_grid(:,:,1) &
& - freq_pot_vort_grid(:,:,2))

! Calculate the streamfunctions for the spectral baroclinic and barotropic
! potential vorticities.
barotropic_pot_vort_strmfunc = spec_inv_barotropic * barotropic_pot_vort
baroclinic_pot_vort_strmfunc = spec_inv_baroclinic * baroclinic_pot_vort

! Calculate the strmfunction for the spectral potential vort.
freq_strmfunc(:,:,1) = barotropic_pot_vort_strmfunc &
& + baroclinic_pot_vort_strmfunc
freq_strmfunc(:,:,2) = barotropic_pot_vort_strmfunc &
& - baroclinic_pot_vort_strmfunc

! Calculate the contributions from the mean shear, rotation deformation, and
! the Ekman friction (the last in layer 2 only).
jacobian_ekman_shear_grid = (0.0_dp, 0.0_dp)
  ! Set layer 1.
jacobian_ekman_shear_grid(:,:,1) = CMPLX((-1.0_dp), 0.0_dp, dp) &
& * CMPLX(vert_shear_sub, 0.0_dp, dp) * spec_x_deriv(:,:,1) &
& * freq_pot_vort_grid(:,:,1) - CMPLX((rotat_wavenum_sub**(2.0_dp) &
  & + vert_shear_sub * deform_wavenum_sub**(2.0_dp)), 0.0_dp, dp) &
& * spec_x_deriv(:,:,1) * freq_strmfunc(:,:,1)
  ! Set layer 2.
jacobian_ekman_shear_grid(:,:,2) = CMPLX((1.0_dp), 0.0_dp, dp) &
& * CMPLX(vert_shear_sub, 0.0_dp, dp) * spec_x_deriv(:,:,1) &
& * freq_pot_vort_grid(:,:,2) - CMPLX((rotat_wavenum_sub**(2.0_dp) &
  & - vert_shear_sub * deform_wavenum_sub**(2.0_dp)), 0.0_dp, dp) &
&* spec_x_deriv(:,:,1) * freq_strmfunc(:,:,2) - CMPLX(ekman_fric_coeff_sub, &
  & 0.0_dp, dp) * spec_laplacian * freq_strmfunc(:,:,2)

! UH OH, SCALING CONTINUE WORKING FROM HERE.

! Must rescale the potential vort and the streamfunction grids to dealias the
! Jacobian, see Orszag, S. "On the Elimination of Aliasing in Finite-Difference 
! Schemes by Filtering High-Wavenumber Components" (1971).
scaled_x_len = 3_qb * x_len_sub/2_qb
scaled_y_len = 3_qb * y_len_sub/2_qb
  ! Rescale the frequency potential vorticity.
ALLOCATE(scaled_freq_pot_vort_grid(scaled_x_len, scaled_y_len, 2))
scaled_freq_pot_vort_grid = (0.0_dp, 0.0_dp)
scaled_freq_pot_vort_grid(1:x_len_sub/2+1, 1:y_len_sub/2+1, :) = &
& (2.25_dp, 0.0_dp) * freq_pot_vort_grid(1:x_len_sub/2+1, 1:y_len_sub/2+1, :)
scaled_freq_pot_vort_grid(1:x_len_sub/2+1, y_len_sub+2:scaled_y_len, :) = &
& (2.25_dp, 0.0_dp) * freq_pot_vort_grid(1:x_len_sub/2+1, &
  & y_len_sub/2+2:y_len_sub, :)
scaled_freq_pot_vort_grid(x_len_sub+2:scaled_x_len, 1:y_len_sub/2+1, :) = &
& (2.25_dp, 0.0_dp) * freq_pot_vort_grid(x_len_sub/2+2:x_len_sub, &
  & 1:y_len_sub/2+1, :)
scaled_freq_pot_vort_grid(x_len_sub+2:scaled_x_len, &
  & y_len_sub+2:scaled_y_len, :) = (2.25_dp, 0.0_dp) &
& * freq_pot_vort_grid(x_len_sub/2+2:x_len_sub, y_len_sub/2+2:y_len_sub, :)
  ! Rescale the potential vorticity streamfunction.
ALLOCATE(scaled_freq_strmfunc(scaled_x_len, scaled_y_len, 2))
scaled_freq_strmfunc = (0.0_dp, 0.0_dp)
scaled_freq_strmfunc(1:x_len_sub/2+1, 1:y_len_sub/2+1, :) = &
& (2.25_dp, 0.0_dp) * freq_strmfunc(1:x_len_sub/2+1, 1:y_len_sub/2+1, :)
scaled_freq_strmfunc(1:x_len_sub/2+1, y_len_sub+2:scaled_y_len, :) = &
& (2.25_dp, 0.0_dp) * freq_strmfunc(1:x_len_sub/2+1, y_len_sub/2+2:y_len_sub, :)
scaled_freq_strmfunc(x_len_sub+2:scaled_x_len, 1:y_len_sub/2+1, :) = &
& (2.25_dp, 0.0_dp) * freq_strmfunc(x_len_sub/2+2:x_len_sub, 1:y_len_sub/2+1,:)
scaled_freq_strmfunc(x_len_sub+2:scaled_x_len, y_len_sub+2:scaled_y_len, :) = &
& (2.25_dp, 0.0_dp) * freq_strmfunc(x_len_sub/2+2:x_len_sub, &
  & y_len_sub/2+2:y_len_sub, :)

! To avoid convolution, we will calculate the Jacobian in physical space, and
! then transform it back to freq space. We want
! J(streamfunction, potential vorticity).
  ! Calculate x-derivative of the potential vorticity.
ALLOCATE(scaled_strmfunc_x_deriv(scaled_x_len, scaled_y_len, 2))
scaled_strmfunc_x_deriv = scaled_spec_x_deriv * scaled_freq_strmfunc
CALL CFFT2DB(scaled_x_len, scaled_y_len, scaled_strmfunc_x_deriv(:,:,1))
CALL CFFT2DB(scaled_x_len, scaled_y_len, scaled_strmfunc_x_deriv(:,:,2))
  ! Take only the real part to get rid of machine-epsilon errors from the
  ! inverse FFT.
scaled_strmfunc_x_deriv = REAL(scaled_strmfunc_x_deriv)

  ! Calculate y-derivative of the streamfunction.
ALLOCATE(scaled_strmfunc_y_deriv(scaled_x_len, scaled_y_len, 2))
scaled_strmfunc_y_deriv = scaled_spec_y_deriv * scaled_freq_strmfunc
CALL CFFT2DB(scaled_x_len, scaled_x_len, scaled_strmfunc_y_deriv(:,:,1))
CALL CFFT2DB(scaled_y_len, scaled_y_len, scaled_strmfunc_y_deriv(:,:,2))
  ! Take only the real part to get rid of machine-epsilon errors from the
  ! inverse FFT.
scaled_strmfunc_y_deriv = REAL(scaled_strmfunc_y_deriv)

  ! Calculate x-derivative of the streamfunction.
ALLOCATE(scaled_pot_vort_x_deriv(scaled_x_len, scaled_y_len, 2))
scaled_pot_vort_x_deriv = scaled_spec_x_deriv * scaled_freq_pot_vort_grid
CALL CFFT2DB(scaled_x_len, scaled_y_len, scaled_pot_vort_x_deriv(:,:,1))
CALL CFFT2DB(scaled_x_len, scaled_y_len, scaled_pot_vort_x_deriv(:,:,2))
  ! Take only the real part to get rid of machine-epsilon errors from the
  ! inverse FFT.
scaled_pot_vort_x_deriv = REAL(scaled_pot_vort_x_deriv)

  ! Calculate y-derivative of the potential vorticity.
ALLOCATE(scaled_pot_vort_y_deriv(scaled_x_len, scaled_y_len, 2))
scaled_pot_vort_y_deriv = scaled_spec_y_deriv * scaled_freq_pot_vort_grid
CALL CFFT2DB(scaled_x_len, scaled_y_len, scaled_pot_vort_y_deriv(:,:,1))
CALL CFFT2DB(scaled_x_len, scaled_y_len, scaled_pot_vort_y_deriv(:,:,2))
  ! Take only the real part to get rid of machine-epsilon errors from the
  ! inverse FFT.
scaled_pot_vort_y_deriv = REAL(scaled_pot_vort_y_deriv)

  ! Calculate the actual Jacobian.
ALLOCATE(scaled_freq_jacobian(scaled_x_len, scaled_y_len, 2))
scaled_freq_jacobian = scaled_strmfunc_x_deriv * scaled_pot_vort_y_deriv &
& - scaled_strmfunc_y_deriv * scaled_pot_vort_x_deriv
CALL CFFT2DF(scaled_x_len, scaled_y_len, scaled_freq_jacobian(:,:,1))
CALL CFFT2DF(scaled_x_len, scaled_y_len, scaled_freq_jacobian(:,:,2))
! The larger grid size means the entries of the scaled Jacobian are 9/4 times
! larger than they should be, so we must correct that.
scaled_freq_jacobian = CMPLX((4.0_dp/9.0_dp), 0.0_dp, dp) * scaled_freq_jacobian


! Reduce the Jacobian to the original grid size.
freq_jacobian(1:x_len_sub/2+1, 1:y_len_sub/2+1, :) = &
scaled_freq_jacobian(1:x_len_sub/2+1, 1:y_len_sub/2+1, :)
freq_jacobian(1:x_len_sub/2+1, y_len_sub/2+2:y_len_sub, :) = &
scaled_freq_jacobian(1:x_len_sub/2+1, y_len_sub+2:scaled_y_len, :)
freq_jacobian(x_len_sub/2+2:x_len_sub, 1:y_len_sub/2+1, :) = &
scaled_freq_jacobian(x_len_sub+2:scaled_x_len, 1:y_len_sub/2+1, :)
freq_jacobian(x_len_sub/2+2:x_len_sub, y_len_sub/2+2:y_len_sub, :) = &
scaled_freq_jacobian(x_len_sub+2:scaled_x_len, y_len_sub+2:scaled_y_len, :)

! Deallocate all scaled arrays.
DEALLOCATE(scaled_freq_pot_vort_grid)
DEALLOCATE(scaled_freq_strmfunc)
DEALLOCATE(scaled_strmfunc_x_deriv)
DEALLOCATE(scaled_strmfunc_y_deriv)
DEALLOCATE(scaled_pot_vort_x_deriv)
DEALLOCATE(scaled_pot_vort_y_deriv)
DEALLOCATE(scaled_freq_jacobian)

! Add in the Jacobian to the output matrix.
jacobian_ekman_shear_grid = jacobian_ekman_shear_grid - freq_jacobian

END SUBROUTINE JACOBIAN_EKMAN_SHEAR
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Calculates the x-direction derivative operator in Fourier space and stores it
! in spec_x_deriv.
!
! STRUCTURE: 1) Calculate the x-wavenumbers.
! 2) If x_len_sub even and order odd, zero out the x_len_sub/2 wavenumber.
! 3) Populate spec_x_deriv array with (sqrt(-1)*x wavenumber)**order.
!
! VARIABLES: - x_len_sub, y_len_sub: The shape of the grid (INTEGER(qb)).
! - spec_x_deriv: The array that will store the matrix used in calculating
! the derivative of the numerical solution along the x-direction of the
! grid (COMPLEX(dp), DIMENSION(x_len_sub, y_len_sub)).
! - order: The order of the derviate desired (INTEGER(qb)).
! - i: Counting index used in DO loops (INTEGER(qb)).
! - x_wavenums: Array used to store the x-direction wavenumbers (REAL(dp),
! DIMENSION(x_len_sub)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE SPECTRAL_X_DERIVATIVE(x_len_sub, y_len_sub, spec_x_deriv, order)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IMPLICIT NONE

INTEGER(qb) :: x_len_sub, &
& y_len_sub, &
& order, &
& i
REAL(dp), DIMENSION(x_len_sub) :: x_wavenums
COMPLEX(dp), DIMENSION(x_len_sub, y_len_sub) :: spec_x_deriv

! Calculate the x-direction wavenumbers.
x_wavenums = 0.0_dp
x_wavenums(2) = 1.0_dp
x_wavenums(x_len_sub) = -1.0_dp
DO i = 1_qb, x_len_sub/2_qb - 1_qb
  x_wavenums(x_len_sub - i) = REAL(-i - 1_qb, dp)
  x_wavenums(i + 2_qb) = REAL(i + 1_qb, dp)
END DO
! If x_len_sub even and order odd, we have to zero out highest wavenumber for
! derivative.
IF ((MOD(x_len_sub, 2_qb) .EQ. 0_qb) .AND. (MOD(order, 2_qb) .EQ. 1_qb)) THEN
  x_wavenums(x_len_sub/2_qb + 1_qb) = 0.0_dp
END IF

DO i = 1, y_len_sub
  spec_x_deriv(:,i) = (CMPLX(0.0, x_wavenums(:), dp))**(REAL(order, dp))
END DO

END SUBROUTINE SPECTRAL_X_DERIVATIVE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Calculates the y-direction derivative operator in Fourier space and stores it
! in spec_y_deriv.
!
! STRUCTURE: 1) Calculate the y-wavenumbers.
! 2) If y_len_sub even and order odd, zero out the y_len_sub/2 wavenumber.
! 3) Populate spec_y_deriv array with (sqrt(-1)*y wavenumber)**order.
!
! VARIABLES: - x_len_sub, y_len_sub: The shape of the grid (INTEGER(qb)).
! - spec_x_deriv: The array that will store the matrix used in calculating
! the derivative of the numerical solution along the x-direction of the
! grid (COMPLEX(dp), DIMENSION(x_len_sub, y_len_sub)).
! - order: The order of the derviate desired (INTEGER(qb)).
! - i: Counting index used in DO loops (INTEGER(qb)).
! - y_wavenums: Array used to store the y-direction wavenumbers (REAL(dp),
! DIMENSION(x_len_sub)).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE SPECTRAL_Y_DERIVATIVE(x_len_sub, y_len_sub, spec_y_deriv, order)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IMPLICIT NONE

INTEGER(qb) :: x_len_sub, &
& y_len_sub, &
& order, &
& i
REAL(dp), DIMENSION(y_len_sub) :: y_wavenums
COMPLEX(dp), DIMENSION(x_len_sub, y_len_sub) :: spec_y_deriv

! Calculate the y-direction wavenumbers.
y_wavenums = 0.0_dp
y_wavenums(2) = 1.0_dp
y_wavenums(x_len_sub) = -1.0_dp
DO i = 1_qb, y_len_sub/2_qb - 1_qb
  y_wavenums(y_len_sub - i) = REAL(-i - 1_qb, dp)
  y_wavenums(i + 2_qb) = REAL(i + 1_qb, dp)
END DO
! If x_len_sub even and order odd, we have to zero out highest wavenumber for
! derivative.
IF ((MOD(y_len_sub, 2_qb) .EQ. 0_qb) .AND. (MOD(order, 2_qb) .EQ. 1_qb)) THEN
  y_wavenums(y_len_sub/2_qb + 1_qb) = 0.0_dp
END IF

DO i = 1, x_len_sub
  spec_y_deriv(i,:) = (CMPLX(0.0, y_wavenums(:), dp))**(REAL(order, dp))
END DO

END SUBROUTINE SPECTRAL_Y_DERIVATIVE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END MODULE
