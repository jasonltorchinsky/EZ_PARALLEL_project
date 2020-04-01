!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> The EZ_PARALLEL finalization subroutine.
!> Deallocates memory for the grid decomposition list, and finalzies MPI.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE FIN_EZP
  
  USE MPI
  
  IMPLICIT NONE

  INTEGER :: ierror

  DEALLOCATE(grid_decomps)
  CALL MPI_FINALIZE(ierror)
 
END SUBROUTINE FIN_EZP



