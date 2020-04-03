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
! DESCRIPTION:
!> \brief The main file for the EZ_PARALLEL module.
!> Contains all derived datatypes and includes the files containing the
!! module subroutines.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE EZ_PARALLEL
  
  USE MPI
  
  IMPLICIT NONE

  PRIVATE

  !> The GRID_DECOMPOSITION derived datatype. It contains all information about
  !! the grid and the sub-grid decomposition of the grid.
  TYPE, PUBLIC :: GRID_DECOMPOSITION
     INTEGER :: decomp_id !< ID number of the grid decomposition, between 1 and
     !! the total number of "unique" grid decompositions.
     INTEGER :: row_count_g !< Number of rows in the grid.
     INTEGER :: col_count_g !< Number of columns in the grid.
     INTEGER :: overlap !< Number of extra columns needed by each sub-grid to
     !! successfully step forward in time.
     INTEGER, ALLOCATABLE :: row_decomp(:) !< Number of rows in each sub-grid
     !! of the horizontal-slab decomposition, excluding overlap.
     INTEGER, ALLOCATABLE :: col_decomp(:) !< Number of columns in each
     !! sub-grid of the vertical-slab decomposition, excluding overlap.
     INTEGER, ALLOCATABLE :: col_decomp_ovlp(:) !< Number of columns in each
     !! sub-grid of the vertical slab decomposition, including overlap.
  END TYPE GRID_DECOMPOSITION

  INTEGER, PUBLIC :: decomp_count !< Number of "unique" grid decompositions.
  INTEGER, PUBLIC :: proc_id !< The local processor ID.
  INTEGER, PUBLIC :: proc_count !< The total number of processors called by the
  !! user.
  TYPE(GRID_DECOMPOSITION), ALLOCATABLE, PUBLIC :: grid_decomps(:) !< A list of all
  !! grid decompositions.

  PUBLIC :: INIT_EZP !< Initialization subroutine for the EZ_PARALLEL module.
  PUBLIC :: FIN_EZP !< Finalization subroutine for the EZ_PARALLEL module.
  PUBLIC :: DECOMP_GRID_EZP !< Grid decomposition for the EZ_PARALLEL module.
  PUBLIC :: IDENTIFY_REF_POINT_EZP !< Sub-grid reference point identification
  !! for the EZ_PARALLEL module.

CONTAINS

  INCLUDE 'init_ezp.f90'
  INCLUDE 'fin_ezp.f90'
  INCLUDE 'decomp_grid_ezp.f90'
  INCLUDE 'identify_ref_point_ezp.f90'
 
END MODULE EZ_PARALLEL
