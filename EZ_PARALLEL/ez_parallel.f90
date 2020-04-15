!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     The GENESIS Project of EZ_PARALLEL, to create the initial version.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : EZ_PARALLEL
! PROJECT          : GENESIS
! MODULE           : MAIN
! URL              : https://github.com/jasonlturner/EZ_PARALLEL_project
! AFFILIATION      : University of Wisconsin-Madison
! DATE             : Spring 2020
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
  PUBLIC :: FIN_EZP 
  PUBLIC :: DECOMP_GRID_EZP 
  PUBLIC :: IDENTIFY_REF_POINT_EZP 
  PUBLIC :: SHARE_SUBGRID_BOUNDARIES_EZP
  PUBLIC :: SPECTRAL_DIM1_DERIVATIVE_EZP


  ! Interface blocks for each of the EZ_PARALLEL subroutines.

  !> \public
  !> @brief The initialization interface for the <tt>EZ_PARALLEL</tt> module.
  !
  !> The <tt>INIT_EZP</tt> interface initializes the <tt>EZ_PARALLEL</tt> module
  !! by intializing <tt>MPI</tt>, allocating memory to the list of grid
  !! decompositions, and storing variables frequently used throughout the module
  !! subroutines.
  !> For greater detail on the initialization subroutine, please see
  !! init_ezp.f90.
  !> For greater detail on the variables frequently used throughout the module
  !! subroutines (which are accessible for read/write by the user), please see
  !! ez_parallel_structs.f90.
  !
  !> \remark Two grid decompositions differ for the sake of
  !! <tt>decomp_count_L</tt> if they differ in any of:
  !! \li size (the number of rows and/or columns),
  !! \li overlap parameter (the number of columns one sub-grid must borrow from
  !! another to perform a time-step),
  !! \li datatype they contain (<tt>DOUBLE PRECISION</tt> or
  !! <tt>DOUBLE COMPLEX</tt>).
  !!
  !! Hence, if the user will need three grid decompositions, they should use
  !! <tt>CALL INIT_EZP(3)</tt>.
  INTERFACE INIT_EZP
     !> \param[in] decomp_count_L Number of unique grid decompositions,
     !! "unique" being in terms of size and overlap (local to subroutine).
     SUBROUTINE INIT(decomp_count_L)
       INTEGER, INTENT(IN) :: decomp_count_L
     END SUBROUTINE INIT
  END INTERFACE INIT_EZP


  !> \public
  !> @brief The finalization interface for the <tt>EZ_PARALLEL</tt> module.
  !
  !> The <tt>FIN_EZP</tt> interface finalizes the EZ_PARALLEL module by
  !! finalizing <tt>MPI</tt> and deallocating the memory that held the array of
  !! grid decompositions.
  !> For greater detail on the finalization subroutine, please see fin_ezp.f90.
  INTERFACE FIN_EZP
     SUBROUTINE FIN
     END SUBROUTINE FIN
  END INTERFACE FIN_EZP

  !> \public
  !> @brief The grid decomposition interface for the <tt>EZ_PARALLEL</tt>
  !! module.
  !
  !> Given the column count of the grid, returns the column count of the
  !! sub-grid, including overlap (for use of the vertical slab decomposition).
  !! Also stores the data in the grid_decomps list for use in other subroutines.
  !! For greater detail on the grid decomposition subroutine, please see
  !! decomp_grid_ezp.f90. For greater detail on the <tt>grid_decomps</tt> list,
  !! please see ez_parallel_structs.f90.
  !
  !> \remark This subroutine assumes the code sets up the grid using an
  !! <tt>ALLOCATE</tt> call, e.g., <tt>ALLOCATE(grid(row_count,col_count))</tt>.
  !! Thus, to decompose the grid, this interface will simply change the value
  !! of <tt>col_count</tt> to that of the processor's sub-grid. No other edits
  !! need to be made to the original serial code in this regard.
  !> \remark For example, consider a grid with 8 columns to be split across 3
  !! processors, and say the time-stepping scheme requires an overlap parameter
  !! of 2. Then the sub-grid interiors will have 3, 3, and 2 columns,
  !! respectively. Including their boundaries (the overlapping regions), the
  !! sub-grids will have 5, 7, and 4 columns, respectively.
  INTERFACE DECOMP_GRID_EZP
     !> \param[in] row_count Number of rows in the grid.
     !> \param[inout] col_count Number of columns in the grid, changes to number
     !! of columns in the sub-grid including overlap, so that the user's
     !! time-stepping code may remain unchanged.
     !> \param[in] overlap Number of extra columns needed by each sub-grid to
     !! successfully step forward in time.
     !> \param[in] decomp_id ID number of the grid decomposition, between 1 and
     !! the number of "unique" grid decompositions.
     !> \param[in] grid The grid to be allocated and decomposed, for use in
     !! identifying what datatype it is. NOTE: It must not yet be allocated.
     SUBROUTINE DECOMP_GRID_DBLE(row_count, col_count, overlap, decomp_id, grid)
       INTEGER, INTENT(IN) :: row_count
       INTEGER, INTENT(INOUT) :: col_count
       INTEGER, INTENT(IN) :: overlap
       INTEGER, INTENT(IN) :: decomp_id
       DOUBLE PRECISION, ALLOCATABLE, INTENT(IN) :: grid(:,:)
     END SUBROUTINE DECOMP_GRID_DBLE
     !> \param[in] row_count Number of rows in the grid.
     !> \param[inout] col_count Number of columns in the grid, changes to number
     !! of columns in the sub-grid including overlap, so that the user's
     !! time-stepping code may remain unchanged.
     !> \param[in] overlap Number of extra columns needed by each sub-grid to
     !! successfully step forward in time.
     !> \param[in] decomp_id ID number of the grid decomposition, between 1 and
     !! the number of "unique" grid decompositions.
     !> \param[in] grid The grid to be allocated and decomposed, for use in
     !! identifying what datatype it is. NOTE: It must not yet be allocated.
     SUBROUTINE DECOMP_GRID_DCMPLX(row_count, col_count, overlap, decomp_id, &
          grid)
       INTEGER, INTENT(IN) :: row_count
       INTEGER, INTENT(INOUT) :: col_count
       INTEGER, INTENT(IN) :: overlap
       INTEGER, INTENT(IN) :: decomp_id
       DOUBLE COMPLEX, ALLOCATABLE, INTENT(IN) :: grid(:,:)
     END SUBROUTINE DECOMP_GRID_DCMPLX
  END INTERFACE DECOMP_GRID_EZP

  !> \public
  !! @brief The reference point identification interface for the
  !! <tt>EZ_PARALLEL</tt> module.
  !
  !> Assuming a linear spacing of grid points, adjusts the physical position of
  !! the reference point in the physical direction corresponding to the second
  !! dimension of the grid.
  !> For greater detail on the reference point
  !! identification subroutine, please see identify_ref_point_ezp.f90.
  !
  !> \remark The reference point of the grid is the physical position of the
  !! first entry of the grid, i.e., the point in the top-left corner. This
  !! reference point is used to ensure the proper initialization of the
  !! initial condition for each sub-grid.
  !> \remark Continuing the example from the <tt>decomp_grid</tt> interface,
  !! say that the physical spacing between columns of the grid is 0.25, and
  !! that the original reference point is at a physical position of 0.1 in the
  !! corresponding direction. Then the reference point of the sub-grids will
  !! be at 0.1, 0.35, and 1.1, respectively.
  INTERFACE IDENTIFY_REF_POINT_EZP
     !> \param[in] col_spc Physical spacing of columns in the grid, e.g., if the
     !! second dimension of the grid represents the z-direction, then col_spc
     !! corresponds to the spacing of grid points in the z-direction.
     !> \param[inout] col_ref Physical position of the reference point in the
     !! second dimension of the grid, e.g., the z-coordinate.
     !> \param[in] decomp_id ID number of the grid decomposition, between 1 and
     !! the number of "unique" grid decompositions.
     SUBROUTINE IDENTIFY_REF_POINT_DBLE(col_spc, col_ref, decomp_id)
       INTEGER, INTENT(IN) :: decomp_id
       DOUBLE PRECISION, INTENT(IN) :: col_spc
       DOUBLE PRECISION, INTENT(INOUT) :: col_ref
     END SUBROUTINE IDENTIFY_REF_POINT_DBLE
     !> \param[in] col_spc Physical spacing of columns in the grid, e.g., if the
     !! second dimension of the grid represents the z-direction, then col_spc
     !! corresponds to the spacing of grid points in the z-direction.
     !> \param[inout] col_ref Physical position of the reference point in the
     !! second dimension of the grid, e.g., the z-coordinate.
     !> \param[in] decomp_id ID number of the grid decomposition, between 1 and
     !! the number of "unique" grid decompositions.
     SUBROUTINE IDENTIFY_REF_POINT_DCMPLX(col_spc, col_ref, decomp_id)
       INTEGER, INTENT(IN) :: decomp_id
       DOUBLE COMPLEX, INTENT(IN) :: col_spc
       DOUBLE COMPLEX, INTENT(INOUT) :: col_ref
     END SUBROUTINE IDENTIFY_REF_POINT_DCMPLX
  END INTERFACE IDENTIFY_REF_POINT_EZP

  !> \public
  !! @brief The sub-grid boundary communication interface for the
  !! <tt>EZ_PARALLEL</tt> module.
  !
  !> For greater detail on the sub-grid boundary communication subroutine,
  !! please see share_subgrid_boundaries_ezp.f90.
  !
  !> Communicates the sub-grid boundary to neighboring sub-grids.
  INTERFACE SHARE_SUBGRID_BOUNDARIES_EZP
     !> \param[in] decomp_id ID number of the grid decomposition, between 1 and
     !! the number of "unique" grid decompositions.
     !> \param[in] sub_grid The sub-grid belonging to the processor. This is
     !! stored under the original variable for the grid, e.g., if the serial
     !! code uses heat_grid as the grid, then the heat_grid should be passed to
     !! this subroutine as the sub_grid argument.
     SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DBLE(decomp_id, sub_grid)
       INTEGER, INTENT(IN) :: decomp_id
       DOUBLE PRECISION, INTENT(INOUT) :: sub_grid(:,:)
     END SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DBLE
     SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DCMPLX(decomp_id, sub_grid)
       INTEGER, INTENT(IN) :: decomp_id
       DOUBLE COMPLEX, INTENT(INOUT) :: sub_grid(:,:)
     END SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DCMPLX
  END INTERFACE SHARE_SUBGRID_BOUNDARIES_EZP

  !> \public
  !! @brief The first-dimension wavenumber array generator interface for the
  !! <tt>EZ_PARALLEL</tt> module.
  !
  !> Generates an array containing the wavenumbers along the first dimension
  !! corresponding to the local processor sub-grid in spectral space, for use in
  !! finding the derivative along the first dimension in spectral space.
  !> For greater detail on the first-dimension wavenumber array generator
  !! subroutine, please see spectral_dim1_derivative_ezp.f90.
  !> \remark We use the convention of wavenumbers going from 0, ...,
  !! <tt>row_count</tt>/2, -<tt>row_count</tt>/2+1, ..., -1.
  INTERFACE SPECTRAL_DIM1_DERIVATIVE_EZP
     !> \param[in] order The order of derivative desired.
     !> \param[in] decomp_id ID number of the grid decomposition, between 1 and the
     !! number of "unique" grid decompositions.
     !> \param[inout] spec_dim1_deriv The array to hold the wavenumbers along the
     !! first dimension.
     SUBROUTINE SPECTRAL_DIM1_DERIVATIVE(order, decomp_id, spec_dim1_deriv)
       INTEGER, INTENT(IN) :: order
       INTEGER, INTENT(IN) :: decomp_id
       DOUBLE COMPLEX, INTENT(INOUT) :: spec_dim1_deriv(:,:)
     END SUBROUTINE SPECTRAL_DIM1_DERIVATIVE
  END INTERFACE SPECTRAL_DIM1_DERIVATIVE_EZP
  
  
  
END MODULE EZ_PARALLEL

INCLUDE 'init_ezp.f90'
INCLUDE 'fin_ezp.f90'
INCLUDE 'decomp_grid_ezp.f90'
INCLUDE 'identify_ref_point_ezp.f90'
INCLUDE 'share_subgrid_boundaries_ezp.f90'
INCLUDE 'spectral_dim1_derivative_ezp.f90'
