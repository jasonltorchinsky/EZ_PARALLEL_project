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

PROGRAM DECOMP_GRID_UNIT_TEST

  USE EZ_PARALLEL
  
  IMPLICIT NONE

  ! Defines standard integer-, real-precision types.
  INCLUDE 'integer_types.h'
  INCLUDE 'real_types.h'

  REAL(dp) :: start_time, & !< Start time of the unit test.
       end_time !< End time of the unit test.

  CALL CPU_TIME(start_time)
  CALL DECOMP_GRID_TEST
  CALL CPU_TIME(end_time)

  PRINT *, "Execution time: ", end_time - start_time, "."
  PRINT *, "DECOMP_GRID_EZP unit test complete. Normal termination."

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Unit test for the EZ_PARALLEL grid decomposition subroutine.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE DECOMP_GRID_TEST

    USE MPI

    IMPLICIT NONE

    INTEGER(qb) :: decomp_count, &
         row_count, &
         col_count, &
         overlap, &
         decomp_id, &
         ierror, &
         i

    NAMELIST /test_params/ decomp_count, row_count, col_count, overlap, decomp_id
    OPEN(1000, file = 'NAMELIST')
    READ(1000, nml = test_params)
    CLOSE(1000)


    CALL INIT_EZP(decomp_count)

    ! Print the input parameters.
    IF (proc_id .EQ. 0_qb) THEN
       PRINT *, "decomp_count: ", decomp_count
       PRINT *, "row_count: ", row_count, " col_count: ", col_count, &
            " overlap: ", overlap, " decomp_id: ", decomp_id
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

    ! Make decomposition decomp_id of grid.
    CALL DECOMP_GRID_EZP(row_count, col_count, overlap, decomp_id)
    
    ! Output decomposition parameters.
    IF (proc_id .EQ. 0) THEN
       PRINT *, "Grid decomposition parameters: "
       PRINT *, "decomp_id: ", grid_decomps(1)%decomp_id
       PRINT *, "row_count_g: ", grid_decomps(1)%row_count_g
       PRINT *, "col_count_g: ", grid_decomps(1)%col_count_g
       PRINT *, "overlap: ", grid_decomps(1)%overlap
       PRINT *, "row_decomp(:): ", grid_decomps(1)%row_decomp
       PRINT *, "col_decomp(:): ", grid_decomps(1)%col_decomp
       PRINT *, "col_decomp_ovlp(:): ", grid_decomps(1)%col_decomp_ovlp
       PRINT *, "New col_count: ", col_count
    END IF

    CALL FIN_EZP

    RETURN

  END SUBROUTINE DECOMP_GRID_TEST

END PROGRAM DECOMP_GRID_UNIT_TEST

