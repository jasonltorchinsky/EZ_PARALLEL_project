!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Progresses the simulation forward in time, outputting the potential
!  vort grid at desired freq.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE JACOBIAN_EKMAN_SHEAR_SOLVE

IMPLICIT NONE

PRIVATE

! Defines standard integer-, real-precision types.
INCLUDE 'integer_types.h'
INCLUDE 'real_types.h'

LOGICAL :: ran = .FALSE.
COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: spectral_x_deriv, &
spectral_y_deriv, scaled_spectral_x_deriv, scaled_spectral_y_deriv
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: spectral_laplacian, &
spectral_inv_barotropic, spectral_inv_baroclinic

PUBLIC :: JACOBIAN_EKMAN_SHEAR

CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE JACOBIAN_EKMAN_SHEAR(freq_pot_vort, grid_len, dt, &
deform_wavenum, rotat_wavenum, vert_shear, ekman_fric_coeff, &
jacobian_ekman_shear_grid)
!  The main driver for the model's time integration. It calls the routine
!  tstep, which updates variables. It writes output files.
!
!  Variables:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE INITIALIZE

IMPLICIT NONE

INTEGER(qb) :: grid_len, scaled_grid_len, i
REAL(dp) :: dt, deform_wavenum, rotat_wavenum, vert_shear, &
ekman_fric_coeff
COMPLEX(dp), DIMENSION(grid_len, grid_len, 2) :: freq_pot_vort, &
jacobian_ekman_shear_grid, freq_strmfunc, freq_jacobian
COMPLEX(dp), DIMENSION(grid_len, grid_len) :: barotropic_pot_vort, &
baroclinic_pot_vort, barotropic_pot_vort_strmfunc, baroclinic_pot_vort_strmfunc
REAL(dp), DIMENSION(:), ALLOCATABLE :: wavenums_1d
COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: scaled_freq_pot_vort, &
scaled_freq_strmfunc, scaled_strmfunc_x_deriv, &
scaled_strmfunc_y_deriv, scaled_pot_vort_x_deriv, &
scaled_pot_vort_y_deriv, scaled_freq_jacobian

! Initialize the spectral derivative, spectral Laplacian, spectral inverse
! barotropic, and spectral inverse baroclinic operators.
IF (.NOT. ran) THEN
  ! Calculate 1D wavenumbers.
  ALLOCATE(wavenums_1d(grid_len))
  wavenums_1d(:) = 0.0_dp
  wavenums_1d(2) = 1.0_dp
  wavenums_1d(grid_len) = -1.0_dp
  DO i = 1, grid_len/2_qb - 1_qb
    wavenums_1d(grid_len - i) = REAL(-i, dp) - 1.0_dp
    wavenums_1d(i + 2_qb) = REAL(i, dp) + 1.0_dp
  END DO

  ! Set x-, y-wavenumbers. Note, multiplication by x-, y-wavenumbers in Fourier
  ! space corresponds to x-, y-derivative in physical space. These are the
  ! correct ones when forming the Laplacian operator, but will need to be
  ! changed later.
  ALLOCATE(spectral_x_deriv(grid_len, grid_len, 2))
  ALLOCATE(spectral_y_deriv(grid_len, grid_len, 2))
  DO i = 1, grid_len
    spectral_x_deriv(:,i,1) = (0.0_dp, 1.0_dp) &
    * CMPLX(wavenums_1d(:), 0.0_dp, dp)
    spectral_y_deriv(:,i,1) = (0.0_dp, 1.0_dp) &
    * CMPLX(wavenums_1d(i), 0.0_dp, dp)
  END DO
  spectral_x_deriv(:,:,2) = spectral_x_deriv(:,:,1)
  spectral_y_deriv(:,:,2) = spectral_y_deriv(:,:,1)

  ! Set spectral Laplacian operator.
  ALLOCATE(spectral_laplacian(grid_len, grid_len))
  spectral_laplacian = REAL((spectral_x_deriv(:,:,1))**2 &
  + (spectral_y_deriv(:,:,1))**2, dp)

  ! Set the inverse barotropic operator (barotropic potential vort is the
  ! Laplacian of the barotropic strmfunction).
  ALLOCATE(spectral_inv_barotropic(grid_len, grid_len))
  spectral_inv_barotropic = 1.0_dp/spectral_laplacian(:,:)
  spectral_inv_barotropic(1,1) = 0.0_dp

  ! Set the inverse baroclinic operator (baroclinic potential vort is the
  ! difference of the Laplacian of the baroclinic strmfunction and
  ! deform_wavenumber**2 * baroclinic strmfunction).
  ALLOCATE(spectral_inv_baroclinic(grid_len, grid_len))
  spectral_inv_baroclinic = 1.0_dp/(spectral_laplacian(:,:) &
  - deform_wavenum**2)
  spectral_inv_baroclinic(1,1) = 0.0_dp

  ! Correcting the x-, y-wavenumbers to properly correspond to differentiation
  ! in physical space.
  wavenums_1d(grid_len/2_qb + 1_qb) = 0.0_dp
  DO i = 1, grid_len
    spectral_x_deriv(:,i,1) = (0.0_dp, 1.0_dp) &
    * CMPLX(wavenums_1d(:), 0.0_dp, dp)
    spectral_y_deriv(:,i,1) = (0.0_dp, 1.0_dp) &
    * CMPLX(wavenums_1d(i), 0.0_dp, dp)
  END DO
  spectral_x_deriv(:,:,2) = spectral_x_deriv(:,:,1)
  spectral_y_deriv(:,:,2) = spectral_y_deriv(:,:,1)

  ! In order to de-alias the Jacobian, we will switch to a grid with 3/2 times
  ! the sidelength of the original grid. Thus, we will need rescaled spectral
  ! derivative operators.
  DEALLOCATE(wavenums_1d)
  scaled_grid_len = 3_qb * grid_len/2_qb
  ALLOCATE(wavenums_1d(scaled_grid_len))
  wavenums_1d(:) = 0.0_dp
  wavenums_1d(2) = 1.0_dp
  wavenums_1d(scaled_grid_len) = -1.0_dp
  DO i = 1, scaled_grid_len/2_qb - 2_qb
    wavenums_1d(i + 2_qb) = REAL(i, dp) + 1.0_dp
    wavenums_1d(scaled_grid_len - i) = REAL((-1_qb * i), dp) - 1.0_dp
  END DO

  ALLOCATE(scaled_spectral_x_deriv(scaled_grid_len, scaled_grid_len, 2))
  ALLOCATE(scaled_spectral_y_deriv(scaled_grid_len, scaled_grid_len, 2))
  DO i = 1, scaled_grid_len
    scaled_spectral_x_deriv(:,i,1) = (0.0_dp, 1.0_dp) &
    * CMPLX(wavenums_1d(:), 0.0_dp, dp)
    scaled_spectral_y_deriv(:,i,1) = (0.0_dp, 1.0_dp) &
    * CMPLX(wavenums_1d(i), 0.0_dp, dp)
  END DO
  scaled_spectral_x_deriv(:,:,2) = scaled_spectral_x_deriv(:,:,1)
  scaled_spectral_y_deriv(:,:,2) = scaled_spectral_y_deriv(:,:,1)

  ran = .TRUE.
END IF

! Calculate the spectral baroclinic and barotropic potential vorticities.
barotropic_pot_vort = 0.5_dp * (freq_pot_vort(:,:,1) &
+ freq_pot_vort(:,:,2))
baroclinic_pot_vort = 0.5_dp * (freq_pot_vort(:,:,1) &
- freq_pot_vort(:,:,2))

! Calculate the streamfunctions for the spectral baroclinic and barotropic
! potential vorticities.
barotropic_pot_vort_strmfunc = spectral_inv_barotropic &
* barotropic_pot_vort
baroclinic_pot_vort_strmfunc = spectral_inv_baroclinic &
* baroclinic_pot_vort

! Calculate the strmfunction for the spectral potential vort.
freq_strmfunc(:,:,1) = &
barotropic_pot_vort_strmfunc + baroclinic_pot_vort_strmfunc
freq_strmfunc(:,:,2) = &
barotropic_pot_vort_strmfunc - baroclinic_pot_vort_strmfunc

! Calculate the constributions from the mean shear, rotation deformation, and
! the Ekman friction (in layer 2 only).
jacobian_ekman_shear_grid = (0.0_dp, 0.0_dp)
  ! Set layer 1.
jacobian_ekman_shear_grid(:,:,1) = CMPLX((-1.0_dp), 0.0_dp, dp) &
* CMPLX(vert_shear, 0.0_dp, dp) * spectral_x_deriv(:,:,1) * freq_pot_vort(:,:,1) &
- CMPLX((rotat_wavenum**2 + vert_shear * deform_wavenum**2), 0.0_dp, dp) &
* spectral_x_deriv(:,:,1) * freq_strmfunc(:,:,1)
  ! Set layer 2.
jacobian_ekman_shear_grid(:,:,2) = CMPLX((1.0_dp), 0.0_dp, dp) &
* CMPLX(vert_shear, 0.0_dp, dp) * spectral_x_deriv(:,:,1) * freq_pot_vort(:,:,2) &
- CMPLX((rotat_wavenum**2 - vert_shear * deform_wavenum**2), 0.0_dp, dp) &
* spectral_x_deriv(:,:,1) * freq_strmfunc(:,:,2) &
- CMPLX(ekman_fric_coeff, 0.0_dp, dp) * spectral_laplacian &
* freq_strmfunc(:,:,2)

! Must rescale the potential vort and the streamfunction grids to dealias the
! Jacobian, see Orszag, S. "On the Elimination of Aliasing in Finite-Difference 
! Schemes by Filtering High-Wavenumber Components" (1971).
scaled_grid_len = 3_qb * grid_len/2_qb
  ! Rescale the frequency potential vorticity.
ALLOCATE(scaled_freq_pot_vort(scaled_grid_len, scaled_grid_len, 2))
scaled_freq_pot_vort = (0.0_dp, 0.0_dp)
scaled_freq_pot_vort(1:grid_len/2+1, 1:grid_len/2+1, :) = &
(2.25_dp, 0.0_dp) * freq_pot_vort(1:grid_len/2+1, 1:grid_len/2+1, :)
scaled_freq_pot_vort(1:grid_len/2+1, grid_len+2:scaled_grid_len, :) = &
(2.25_dp, 0.0_dp) * freq_pot_vort(1:grid_len/2+1, grid_len/2+2:grid_len, :)
scaled_freq_pot_vort(grid_len+2:scaled_grid_len, 1:grid_len/2+1, :) = &
(2.25_dp, 0.0_dp) * freq_pot_vort(grid_len/2+2:grid_len, 1:grid_len/2+1, :)
scaled_freq_pot_vort(grid_len+2:scaled_grid_len, grid_len+2:scaled_grid_len, :) = &
(2.25_dp, 0.0_dp) * freq_pot_vort(grid_len/2+2:grid_len, grid_len/2+2:grid_len, :)
  ! Rescale the potential vorticity streamfunction.
ALLOCATE(scaled_freq_strmfunc(scaled_grid_len, scaled_grid_len, 2))
scaled_freq_strmfunc = (0.0_dp, 0.0_dp)
scaled_freq_strmfunc(1:grid_len/2+1, 1:grid_len/2+1, :) = &
(2.25_dp, 0.0_dp) * freq_strmfunc(1:grid_len/2+1, 1:grid_len/2+1, :)
scaled_freq_strmfunc(1:grid_len/2+1, grid_len+2:scaled_grid_len, :) = &
(2.25_dp, 0.0_dp) * freq_strmfunc(1:grid_len/2+1, grid_len/2+2:grid_len, :)
scaled_freq_strmfunc(grid_len+2:scaled_grid_len, 1:grid_len/2+1, :) = &
(2.25_dp, 0.0_dp) * freq_strmfunc(grid_len/2+2:grid_len, 1:grid_len/2+1,:)
scaled_freq_strmfunc(grid_len+2:scaled_grid_len, grid_len+2:scaled_grid_len, :) = &
(2.25_dp, 0.0_dp) * freq_strmfunc(grid_len/2+2:grid_len, grid_len/2+2:grid_len, :)

! To avoid convolution, we will calculate the Jacobian in physical space, and
! then transform it back to freq space. We want
! J(streamfunction, potential vorticity).
  ! Calculate x-derivative of the potential vorticity.
ALLOCATE(scaled_strmfunc_x_deriv(scaled_grid_len, scaled_grid_len, 2))
scaled_strmfunc_x_deriv = scaled_spectral_x_deriv * scaled_freq_strmfunc
CALL CFFT2DB(scaled_grid_len, scaled_grid_len, scaled_strmfunc_x_deriv(:,:,1))
CALL CFFT2DB(scaled_grid_len, scaled_grid_len, scaled_strmfunc_x_deriv(:,:,2))
  ! Take only the real part to get rid of machine-epsilon errors from the
  ! inverse FFT.
scaled_strmfunc_x_deriv = REAL(scaled_strmfunc_x_deriv)

  ! Calculate y-derivative of the streamfunction.
ALLOCATE(scaled_strmfunc_y_deriv(scaled_grid_len, scaled_grid_len, 2))
scaled_strmfunc_y_deriv = scaled_spectral_y_deriv * scaled_freq_strmfunc
CALL CFFT2DB(scaled_grid_len, scaled_grid_len, scaled_strmfunc_y_deriv(:,:,1))
CALL CFFT2DB(scaled_grid_len, scaled_grid_len, scaled_strmfunc_y_deriv(:,:,2))
  ! Take only the real part to get rid of machine-epsilon errors from the
  ! inverse FFT.
scaled_strmfunc_y_deriv = REAL(scaled_strmfunc_y_deriv)

  ! Calculate x-derivative of the streamfunction.
ALLOCATE(scaled_pot_vort_x_deriv(scaled_grid_len, scaled_grid_len, 2))
scaled_pot_vort_x_deriv = scaled_spectral_x_deriv * scaled_freq_pot_vort
CALL CFFT2DB(scaled_grid_len, scaled_grid_len, scaled_pot_vort_x_deriv(:,:,1))
CALL CFFT2DB(scaled_grid_len, scaled_grid_len, scaled_pot_vort_x_deriv(:,:,2))
  ! Take only the real part to get rid of machine-epsilon errors from the
  ! inverse FFT.
scaled_pot_vort_x_deriv = REAL(scaled_pot_vort_x_deriv)

  ! Calculate y-derivative of the potential vorticity.
ALLOCATE(scaled_pot_vort_y_deriv(scaled_grid_len, scaled_grid_len, 2))
scaled_pot_vort_y_deriv = scaled_spectral_y_deriv * scaled_freq_pot_vort
CALL CFFT2DB(scaled_grid_len, scaled_grid_len, scaled_pot_vort_y_deriv(:,:,1))
CALL CFFT2DB(scaled_grid_len, scaled_grid_len, scaled_pot_vort_y_deriv(:,:,2))
  ! Take only the real part to get rid of machine-epsilon errors from the
  ! inverse FFT.
scaled_pot_vort_y_deriv = REAL(scaled_pot_vort_y_deriv)

  ! Calculate the actual Jacobian.
ALLOCATE(scaled_freq_jacobian(scaled_grid_len, scaled_grid_len, 2))
scaled_freq_jacobian = scaled_strmfunc_x_deriv &
* scaled_pot_vort_y_deriv - scaled_strmfunc_y_deriv &
* scaled_pot_vort_x_deriv
CALL CFFT2DF(scaled_grid_len, scaled_grid_len, scaled_freq_jacobian(:,:,1))
CALL CFFT2DF(scaled_grid_len, scaled_grid_len, scaled_freq_jacobian(:,:,2))
! The larger grid size means the entries of the scaled Jacobian are 9/4 times
! larger than they should be, so we must correct that.
scaled_freq_jacobian = CMPLX((4.0_dp/9.0_dp), 0.0_dp, dp) * scaled_freq_jacobian


! Reduce the Jacobian to the original grid size.
freq_jacobian(1:grid_len/2+1, 1:grid_len/2+1, :) = &
scaled_freq_jacobian(1:grid_len/2+1, 1:grid_len/2+1, :)
freq_jacobian(1:grid_len/2+1, grid_len/2+2:grid_len, :) = &
scaled_freq_jacobian(1:grid_len/2+1, grid_len+2:scaled_grid_len, :)
freq_jacobian(grid_len/2+2:grid_len, 1:grid_len/2+1, :) = &
scaled_freq_jacobian(grid_len+2:scaled_grid_len, 1:grid_len/2+1, :)
freq_jacobian(grid_len/2+2:grid_len, grid_len/2+2:grid_len, :) = &
scaled_freq_jacobian(grid_len+2:scaled_grid_len, grid_len+2:scaled_grid_len, :)

! Deallocate all scaled arrays.
DEALLOCATE(scaled_freq_pot_vort)
DEALLOCATE(scaled_freq_strmfunc)
DEALLOCATE(scaled_strmfunc_x_deriv)
DEALLOCATE(scaled_strmfunc_y_deriv)
DEALLOCATE(scaled_pot_vort_x_deriv)
DEALLOCATE(scaled_pot_vort_y_deriv)
DEALLOCATE(scaled_freq_jacobian)


! Add in the Jacobian to the output matrix.
jacobian_ekman_shear_grid = jacobian_ekman_shear_grid - freq_jacobian

! Thing are still a little off by about 10^(-3)... Somewhere after
! jacobian_ekman_shear_grid first assignment.
!PRINT *, "jacobian_ekman_shear_grid 2 = "
!DO i = 1, 16
!  WRITE(*,"(2D32.16)") jacobian_ekman_shear_grid(3,i,2)
!END DO
!ERROR STOP "Printed requested information."

END SUBROUTINE JACOBIAN_EKMAN_SHEAR

END MODULE
