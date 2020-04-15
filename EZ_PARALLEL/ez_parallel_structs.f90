!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : STRUCTS
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
! REVISION         : ALPHA 1.01
!
!> @author
!> Jason Turner
! DESCRIPTION:
!> \brief The EZ_PARALLEL structures module.
!> Contains all EZ_PARALLEL derived datatypes and global variables.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE EZ_PARALLEL_STRUCTS
  
  IMPLICIT NONE

  PRIVATE

  !> @class GRID_DECOMPOSITTION.
  !! The GRID_DECOMPOSITION derived datatype contains all information about
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
     INTEGER :: SEND_BOUNDARIES(2) !< MPI derived datatype for sending
     !! sub-grid boundaries to neightboring sub-grids (1 = left, 2 = right).
     INTEGER :: RECV_BOUNDARIES(2) !< MPI derived datatype for recieving
     !! sub-grid boundaries from neightboring sub-grids (1 = left, 2 = right).
  END TYPE GRID_DECOMPOSITION

  ! Global variable declarations
  INTEGER, PUBLIC :: decomp_count !< Number of "unique" grid decompositions.
  INTEGER, PUBLIC :: proc_id !< The local processor ID.
  INTEGER, PUBLIC :: proc_count !< The total number of processors called by the
  !! user.
  TYPE(GRID_DECOMPOSITION), ALLOCATABLE, PUBLIC :: grid_decomps(:) !< A list of all
  !! grid decompositions.


  
END MODULE EZ_PARALLEL_STRUCTS
