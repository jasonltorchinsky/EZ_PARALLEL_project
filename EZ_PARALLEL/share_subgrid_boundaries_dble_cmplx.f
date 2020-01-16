!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SHARE_SUBGRID_BOUNDARIES_DBLE_CMPLX. Communicates the sub-grid boundary to
! neighboring sub-grids, for DOUBLE COMPLEX sub-grids.
!
! ARGUMENTS: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - sub_grid: The sub-grid of the processor, or the grid if only one processor
! is in use (DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len)).
!
! Written By: Jason Turner
! Last Updated: January 15, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DBLE_CMPLX_EZP(dim1_len, dim2_len, &
& overlap, sub_grid)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& overlap
DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len) :: sub_grid

! Since DOUBLE COMPLEX arrays are stored in memory as a DOUBLE PRECISION with
! twice the dim1 length, we can just use the SHARE_SUBGRID_BOUNDARIES_DBLE_REAL
! subroutine.
CALL SHARE_SUBGRID_BOUNDARIES_DBLE_REAL_EZP(2*dim1_len, dim2_len, overlap, &
& sub_grid)

END SUBROUTINE SHARE_SUBGRID_BOUNDARIES_DBLE_CMPLX_EZP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
