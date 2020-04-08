!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : UNIT_TEST
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : April 2nd, 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
! DESCRIPTION:
!> \brief The unit test for the DECOMP_GRID_EZP subroutine.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM IDENTIFY_REF_POINT_UNIT_TEST

  USE MPI
  USE EZ_PARALLEL
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  ! Defines standard integer-, real-precision types.
  INCLUDE 'integer_types.h'
  INCLUDE 'real_types.h'

  REAL(dp) :: start_time, & !< Start time of the unit test.
       end_time !< End time of the unit test.

  CALL CPU_TIME(start_time)
  CALL IDENTIFY_REF_POINT_TEST
  CALL CPU_TIME(end_time)

  PRINT *, "Execution time: ", end_time - start_time, "."
  PRINT *, "IDENTIFY_REF_POINT_EZP unit test complete. Normal termination."

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Unit test for the EZ_PARALLEL grid decomposition subroutine.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE IDENTIFY_REF_POINT_TEST

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
    INTEGER(qb) :: ierror
    INTEGER(qb) :: i
    REAL(dp) :: col_spc_dble
    REAL(dp) :: col_ref_dble
    COMPLEX(dp) :: col_spc_dcmplx
    COMPLEX(dp) :: col_ref_dcmplx
    REAL(dp), ALLOCATABLE :: grid_dble(:,:)
    COMPLEX(dp), ALLOCATABLE :: grid_dcmplx(:,:)

    NAMELIST /test_params/ decomp_count_L, row_count_dble, col_count_dble, &
         overlap_dble, decomp_id_dble, row_count_dcmplx, col_count_dcmplx, &
         overlap_dcmplx, decomp_id_dcmplx, col_spc_dble, col_ref_dble, &
         col_spc_dcmplx, col_ref_dcmplx
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
       PRINT *, "col_spc_dble: ", col_spc_dble, " col_ref_dble: ", &
            col_ref_dble
       PRINT *, "col_spc_dcmplx: ", col_spc_dcmplx, " col_ref_dcmplx: ", &
            col_ref_dcmplx
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

    ! Identify the sub-grid reference points for each grid.
    CALL IDENTIFY_REF_POINT_EZP(col_spc_dble, col_ref_dble, decomp_id_dble)
    CALL IDENTIFY_REF_POINT_EZP(col_spc_dcmplx, col_ref_dcmplx, decomp_id_dcmplx)

    PRINT *, "proc_id: ", proc_id, " col_ref_dble: ", col_ref_dble, &
         " col_ref_dcmplx: ", col_ref_dcmplx

    CALL FIN_EZP

    RETURN

  END SUBROUTINE IDENTIFY_REF_POINT_TEST

END PROGRAM IDENTIFY_REF_POINT_UNIT_TEST


