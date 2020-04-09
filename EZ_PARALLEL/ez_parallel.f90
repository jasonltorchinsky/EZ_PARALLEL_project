!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : MAIN
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : April 2nd, 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
!
!> \brief The EZ_PARALLEL module.
!> Contains EZ_PARALLEL subroutines and their interfaces.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE EZ_PARALLEL
  
  USE MPI
  USE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: INIT_EZP 
  PUBLIC :: FIN_EZP !< @public Finalization subroutine for the EZ_PARALLEL module.
  PUBLIC :: DECOMP_GRID_EZP !< @public Grid decomposition for the EZ_PARALLEL module.
  PUBLIC :: IDENTIFY_REF_POINT_EZP !< @public Sub-grid reference point identification
  !! for the EZ_PARALLEL module.
  PUBLIC :: SHARE_SUBGRID_BOUNDARIES_EZP !< @public Sub-grid boundary communication
  !! subroutine for the EZ_PARALLEL module.


  ! Interface blocks for each of the EZ_PARALLEL subroutines.

  !> \public
  !> \brief The initialization interface for the EZ_PARALLEL module.
  
  !> The INIT_EZP interface initializes the EZ_PARALLEL module by
  !! intializing MPI and allocating memory to the list of grid decompositions.
  !! For greater detail, please see init_ezp.f90.
  
  !> \remark Two grid decompositions differ for the sake of <tt>decomp_count_L</tt> if
  !! they differ in any of:
  !! \li size (the number of rows and/or columns),
  !! \li overlap parameter (the number of columns one sub-grid must borrow from
  !! another to perform a time-step),
  !! \li datatype they contain (<tt>DOUBLE PRECISION</tt> or
  !! <tt>DOUBLE COMPLEX</tt>).
  !!
  !! Hence, if the user will need three grid decompositions, they should use
  !! <tt>CALL INIT_EZP(3)</tt>.
  INTERFACE INIT_EZP
     SUBROUTINE INIT(decomp_count_L)
       INTEGER, INTENT(IN) :: decomp_count_L
     END SUBROUTINE INIT
  END INTERFACE INIT_EZP

  INTERFACE FIN_EZP
     SUBROUTINE FIN
     END SUBROUTINE FIN
  END INTERFACE FIN_EZP

  INTERFACE DECOMP_GRID_EZP
     SUBROUTINE DECOMP_GRID_DBLE(row_count, col_count, overlap, decomp_id, grid)
       INTEGER, INTENT(IN) :: row_count
       INTEGER, INTENT(INOUT) :: col_count
       INTEGER, INTENT(IN) :: overlap
       INTEGER, INTENT(IN) :: decomp_id
       DOUBLE PRECISION, ALLOCATABLE, INTENT(IN) :: grid(:,:)
     END SUBROUTINE DECOMP_GRID_DBLE
     SUBROUTINE DECOMP_GRID_DCMPLX(row_count, col_count, overlap, decomp_id, &
          grid)
       INTEGER, INTENT(IN) :: row_count
       INTEGER, INTENT(INOUT) :: col_count
       INTEGER, INTENT(IN) :: overlap
       INTEGER, INTENT(IN) :: decomp_id
       DOUBLE COMPLEX, ALLOCATABLE, INTENT(IN) :: grid(:,:)
     END SUBROUTINE DECOMP_GRID_DCMPLX
  END INTERFACE DECOMP_GRID_EZP

  INTERFACE IDENTIFY_REF_POINT_EZP
     SUBROUTINE IDENTIFY_REF_POINT_DBLE(col_spc, col_ref, decomp_id)
       INTEGER, INTENT(IN) :: decomp_id
       DOUBLE PRECISION, INTENT(IN) :: col_spc
       DOUBLE PRECISION, INTENT(INOUT) :: col_ref
     END SUBROUTINE IDENTIFY_REF_POINT_DBLE
     SUBROUTINE IDENTIFY_REF_POINT_DCMPLX(col_spc, col_ref, decomp_id)
       INTEGER, INTENT(IN) :: decomp_id
       DOUBLE COMPLEX, INTENT(IN) :: col_spc
       DOUBLE COMPLEX, INTENT(INOUT) :: col_ref
     END SUBROUTINE IDENTIFY_REF_POINT_DCMPLX
  END INTERFACE IDENTIFY_REF_POINT_EZP
  
  INTERFACE SHARE_SUBGRID_BOUNDARIES_EZP
     SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DBLE(decomp_id, sub_grid)
       INTEGER, INTENT(IN) :: decomp_id
       DOUBLE PRECISION, INTENT(INOUT) :: sub_grid(:,:)
     END SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DBLE
     SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DCMPLX(decomp_id, sub_grid)
       INTEGER, INTENT(IN) :: decomp_id
       DOUBLE COMPLEX, INTENT(INOUT) :: sub_grid(:,:)
     END SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DCMPLX
  END INTERFACE SHARE_SUBGRID_BOUNDARIES_EZP
  
  
END MODULE EZ_PARALLEL

INCLUDE 'init_ezp.f90'
INCLUDE 'fin_ezp.f90'
INCLUDE 'decomp_grid_ezp.f90'
INCLUDE 'identify_ref_point_ezp.f90'
INCLUDE 'share_subgrid_boundaries_ezp.f90'
