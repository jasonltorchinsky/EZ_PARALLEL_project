!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : UNIT_TEST
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : April 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
! DESCRIPTION:
!> \brief The unit test for the SPECTRAL_DIM1_DERIVATIVE_EZP subroutine.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM SPECTRAL_DIM1_DERIVATIVE_UNIT_TEST

  USE MPI
  USE EZ_PARALLEL
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  ! Defines standard integer-, real-precision types.
  INCLUDE 'integer_types.h'
  INCLUDE 'real_types.h'

  REAL(dp) :: start_time !< Start time of the unit test.
  REAL(dp) :: end_time !< End time of the unit test.

  CALL CPU_TIME(start_time)
  CALL SPECTRAL_DIM1_DERIVATIVE_TEST
  CALL CPU_TIME(end_time)

  PRINT *, "Execution time: ", end_time - start_time, "."
  PRINT *, "SPECTRAL_DIM1_DERIVATIVE_EZP unit test complete. ", &
       "Normal termination."

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Unit test for the <tt>EZ_PARALLEL</tt> first-dimension wavenumber array
  !! generator subroutine.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE SPECTRAL_DIM1_DERIVATIVE_TEST

    USE MPI
    USE EZ_PARALLEL
    USE EZ_PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER(qb) :: decomp_count_L
    INTEGER(qb) :: row_count_dble
    INTEGER(qb) :: col_count_dble
    INTEGER(qb) :: overlap_dble
    INTEGER(qb) :: decomp_id_dble
    INTEGER(qb) :: row_count_dcmplx
    INTEGER(qb) :: col_count_dcmplx
    INTEGER(qb) :: overlap_dcmplx
    INTEGER(qb) :: decomp_id_dcmplx
    INTEGER(qb) :: order_dble
    INTEGER(qb) :: order_dcmplx
    INTEGER(qb) :: ierror
    REAL(dp), ALLOCATABLE :: grid_dble(:,:)
    COMPLEX(dp), ALLOCATABLE :: grid_dcmplx(:,:)
    COMPLEX(dp), ALLOCATABLE :: spec_dim1_deriv_dble(:,:)
    COMPLEX(dp), ALLOCATABLE :: spec_dim1_deriv_dcmplx(:,:)

    NAMELIST /test_params/ decomp_count_L, row_count_dble, col_count_dble, &
         overlap_dble, decomp_id_dble, row_count_dcmplx, col_count_dcmplx, &
         overlap_dcmplx, decomp_id_dcmplx, order_dble, order_dcmplx
    OPEN(1000, file = 'NAMELIST')
    READ(1000, nml = test_params)
    CLOSE(1000)

    CALL INIT_EZP(decomp_count_L)
    
    ! Print the input parameters.
    IF (proc_id .EQ. 0_qb) THEN
       PRINT *, "decomp_count_L: ", decomp_count_L
       PRINT *, "row_count_dble: ", row_count_dble, " col_count_dble: ", &
            col_count_dble, " overlap_dble: ", overlap_dble, &
            " decomp_id_dble: ", decomp_id_dble
       PRINT *, "row_count_dcmplx: ", row_count_dcmplx, " col_count_dcmplx: ", &
            col_count_dcmplx, " overlap_dcmplx: ", overlap_dcmplx, &
            " decomp_id_dcmplx: ", decomp_id_dcmplx
       PRINT *, "order_dble: ", order_dble, " order_dcmplx: ", order_dcmplx
    END IF
    
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

    ! Make decomposition decomp_id of DBLE and DCMPLX grids.
    CALL DECOMP_GRID_EZP(row_count_dble, col_count_dble, overlap_dble, &
         decomp_id_dble, grid_dble)
    CALL DECOMP_GRID_EZP(row_count_dcmplx, col_count_dcmplx, overlap_dcmplx, &
         decomp_id_dcmplx, grid_dcmplx)
    
    ! Output decomposition parameters.
    IF (proc_id .EQ. 0) THEN
       PRINT *, "DBLE grid decomposition parameters: "
       PRINT *, "decomp_id: ", grid_decomps(decomp_id_dble)%decomp_id
       PRINT *, "row_count_g: ", grid_decomps(decomp_id_dble)%row_count_g
       PRINT *, "col_count_g: ", grid_decomps(decomp_id_dble)%col_count_g
       PRINT *, "overlap: ", grid_decomps(decomp_id_dble)%overlap
       PRINT *, "row_decomp(:): ", grid_decomps(decomp_id_dble)%row_decomp
       PRINT *, "col_decomp(:): ", grid_decomps(decomp_id_dble)%col_decomp
       PRINT *, "col_decomp_ovlp(:): ", grid_decomps(decomp_id_dble)%col_decomp_ovlp

       PRINT *, "DCMPLX grid decomposition parameters: "
       PRINT *, "decomp_id: ", grid_decomps(decomp_id_dcmplx)%decomp_id
       PRINT *, "row_count_g: ", grid_decomps(decomp_id_dcmplx)%row_count_g
       PRINT *, "col_count_g: ", grid_decomps(decomp_id_dcmplx)%col_count_g
       PRINT *, "overlap: ", grid_decomps(decomp_id_dcmplx)%overlap
       PRINT *, "row_decomp(:): ", grid_decomps(decomp_id_dcmplx)%row_decomp
       PRINT *, "col_decomp(:): ", grid_decomps(decomp_id_dcmplx)%col_decomp
       PRINT *, "col_decomp_ovlp(:): ", grid_decomps(decomp_id_dcmplx)%col_decomp_ovlp
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
    PRINT *, "proc_id: ", proc_id, " new col_count_dble: ", col_count_dble, &
         " new col_count_dcmplx: ", col_count_dcmplx

    ! Calculate the spectral derivative in the first dimension for each grid
    ALLOCATE(spec_dim1_deriv_dble(row_count_dble,col_count_dble))
    spec_dim1_deriv_dble = (0.0, 0.0)
    CALL SPECTRAL_DIM1_DERIVATIVE_EZP(order_dble, decomp_id_dble, spec_dim1_deriv_dble)
    ALLOCATE(spec_dim1_deriv_dcmplx(row_count_dcmplx,col_count_dcmplx))
    CALL SPECTRAL_DIM1_DERIVATIVE_EZP(order_dcmplx, decomp_id_dcmplx, &
         spec_dim1_deriv_dcmplx)

    PRINT *, " "
    PRINT *, "proc_id: ", proc_id, " spec_dim1_deriv_dble: ", &
         spec_dim1_deriv_dble
    PRINT *, "proc_id: ", proc_id, " spec_dim1_deriv_dcmplx: ", &
         spec_dim1_deriv_dcmplx

    ! Deallocate allocated arrays
    DEALLOCATE(spec_dim1_deriv_dble)
    DEALLOCATE(spec_dim1_deriv_dcmplx)
    
    CALL FIN_EZP

    RETURN

  END SUBROUTINE SPECTRAL_DIM1_DERIVATIVE_TEST
  
END PROGRAM SPECTRAL_DIM1_DERIVATIVE_UNIT_TEST

