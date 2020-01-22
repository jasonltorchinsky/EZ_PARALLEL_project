!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ZERO_PADDING_DBLE_REAL. Scales the input matrix by padding zeros (the
! 3/2-rule) for de-aliasing FFTs.
!
! ARGUMENTS: - dim1_len, dim2_len: The shape of the matrix (INTEGER).
! - matrix: The matrix to be scaled (DOUBLE PRECISION,
! DIMENSION(dim1_len, dim2_len)).
! - scaled_matrix: The resulting scaled matrix (DOUBLE PRECISION,
! DIMENSION(:,:), ALLOCATABLE).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - proc_id: The processor ID (INTEGER).
! - proc_count: The total number of processors (INTEGER).
!
! NOTES: - When discussing the wavenumbers, we use the convention (for dim1)
! [0, ..., dim1_len/2, -dim1_len/1 + 1, ... -1], with FLOOR(dim1_len/2), as
! is assumed in Fortran integer division.
!
! Written By: Jason Turner
! Last Updated: January 22, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE ZERO_PADDING_DBLE_REAL_EZP(dim1_len, dim2_len, matrix, &
& scaled_matrix, overlap)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& overlap, &
& proc_id, &
& proc_count, &
& ierror
DOUBLE PRECISION, DIMENSION(dim1_len, dim2_len) :: matrix
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: scaled_matrix

CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, proc_count, ierror)

! Check for errors in user input.
CALL ERROR_HANDLING(proc_id)

! Run in serial if there is only one processor.
IF (proc_count .EQ. 1) THEN
  CALL SERIAL(dim1_len, dim2_len, matrix, scaled_matrix)
ELSE
  CALL PARALLEL(dim1_len, dim2_len, matrix, scaled_matrix, overlap, &
  & proc_id, proc_count)
END IF


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SERIAL. Handles the execution if there is only one processor.
!
! STRUCTURE: 1) Calculates the size of the scaled matrix, and the index in the
! scaled matrix of the largest negative wavenumber in both dim1 and dim2.
! 2) Fills in the non-zero entries of the scaled matrix (which correspond to the
! wavenumbers of the original matrix) with the entries of the original matrix.
!
! VARIABLES: - dim1_len, dim2_len: The shape of the matrix (INTEGER).
! - matrix: The matrix to be scaled (DOUBLE PRECISION,
! DIMENSION(dim1_len, dim2_len)).
! - scaled_matrix: The resulting scaled matrix (DOUBLE PRECISION,
! DIMENSION(3*dim1_len/2, 3*dim2_len/2)). Note, we use an assumed-shape array
! here.
! - scaled_dim1_len, scaled_dim2_len: The dimensions of scaled_matrix (INTEGER).
! - scaled_dim1_neg_wavenum_low_index, scaled_dim2_neg_wavenum_low_index: The
! lower index of the largest negative wavenumber in scaled_matrix (INTEGER).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE SERIAL(dim1_len, dim2_len, matrix, scaled_matrix)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len,  &
& scaled_dim1_len, &
& scaled_dim2_len, &
& scaled_dim1_neg_wavenum_low_index, &
& scaled_dim2_neg_wavenum_low_index
DOUBLE PRECISION, DIMENSION(dim1_len, dim2_len) :: matrix
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: scaled_matrix

scaled_dim1_len = (3 * dim1_len)/2
scaled_dim2_len = (3 * dim2_len)/2

PRINT *, 'FLAG 0'

ALLOCATE(scaled_matrix(scaled_dim1_len,scaled_dim2_len))
scaled_matrix = 0.0
IF (MOD(dim1_len, 2) .EQ. 0) THEN
  scaled_dim1_neg_wavenum_low_index = dim1_len + 2
ELSE IF (MOD(dim1_len, 2) .EQ. 1) THEN
  scaled_dim1_neg_wavenum_low_index = dim1_len + 1
END IF

IF (MOD(dim2_len, 2) .EQ. 0) THEN
  scaled_dim2_neg_wavenum_low_index = dim2_len + 2
ELSE IF (MOD(dim2_len, 2) .EQ. 1) THEN
  scaled_dim2_neg_wavenum_low_index = dim2_len + 1
END IF

PRINT *, 'FLAG 1'

! This matrix is scale by 3/2 in each dimension, with the indices corresponding
! to the largest wavenumbers (in magnitude) zeroed out.
scaled_matrix(1:dim1_len/2+1, 1:dim2_len/2+1) = &
& matrix(1:dim1_len/2+1, 1:dim2_len/2+1)
PRINT *, 'FLAG 2'
scaled_matrix(1:dim1_len/2+1, &
  & scaled_dim2_neg_wavenum_low_index:scaled_dim2_len) = &
& matrix(1:dim1_len/2+1, dim2_len/2+2:dim2_len)
PRINT *, 'FLAG 3'
scaled_matrix(scaled_dim1_neg_wavenum_low_index:scaled_dim1_len, &
  & 1:dim2_len/2+1) = &
& matrix(dim1_len/2+2:dim1_len, 1:dim2_len/2+1)
PRINT *, 'FLAG 4'
scaled_matrix(scaled_dim1_neg_wavenum_low_index:scaled_dim1_len, &
  & scaled_dim2_neg_wavenum_low_index:scaled_dim2_len) = &
& matrix(dim1_len/2+2:dim1_len, dim2_len/2+2:dim2_len)
PRINT *, 'FLAG 5'

END SUBROUTINE SERIAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PARALLEL. Handles the execution if there is more than one processor.
!
! STRUCTURE: 1) Calculate the size of the total grid in order to obtain the
! proper size for the scaled sub-grid. For ease of calculation, we will obtain
! the total scaled grid by obtaining the interiors of each subgrid, and then
! share their overlapping regions.
!
! VARIABLES: - variable1: variable1 description.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE PARALLEL(dim1_len, dim2_len, matrix, scaled_matrix, overlap, &
& proc_id, proc_count)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& overlap, &
& proc_id, &
& proc_count, &
& dim1_len_total, &
& dim2_len_total, &
& dim1_len_interior, &
& dim2_len_interior, &
& dim2_len_other, &
& dim2_global_index_low, &
& dim2_global_index_high, &
& ierror, &
& i, &
& j
INTEGER, DIMENSION(proc_count) :: dim2_len_interior_list
DOUBLE PRECISION, DIMENSION(dim1_len, dim2_len) :: matrix
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: scaled_matrix, &
& matrix_interior, &
& scaled_matrix_interior

! Obtain the total size of the matrix.
dim1_len_interior = dim1_len
dim1_len_total = dim1_len
! Middle processor.
IF (MOD(proc_id, proc_count-1) .NE. 0) THEN
  dim2_len_interior = dim2_len - 2 * overlap
! Processor 0 or proc_count-1.
ELSE
  dim2_len_interior = dim2_len - overlap
END IF

CALL MPI_ALLREDUCE(dim2_len_interior, dim2_len_total, 1, MPI_INTEGER, MPI_SUM, &
& MPI_COMM_WORLD, ierror)

! Obtain the interior bands of the sub-grids.
ALLOCATE(matrix_interior(dim1_len_interior, dim2_len_interior))
! Processor 0, doesn't need last parts of dim2.
IF (proc_id .EQ. 0) THEN
  matrix_interior = matrix(:,1:dim2_len-overlap)
! Last processor, doesn't need first part of dim2.
ELSE IF (proc_id .EQ. proc_count-1) THEN
  matrix_interior = matrix(:,overlap+1:dim2_len)
! Middle processor, doesn't need first or last part of dim2.
ELSE
  matrix_interior = matrix(:,overlap+1:dim2_len-overlap)
END IF

! To obtain global dim2 coordinates of sub-grid, need list of all dim2_len. We
! recalculate those based on dim2_len_total.

! Divide dim2_len_total of total grid as evenly as possible.
dim2_len_other = dim2_len_total/proc_count
! Now check if sub-grid of a given processor needs to take additional dim2,
! i.e., dim2_len % proc_count != 0.
DO i = 1, proc_count
  IF ((i-1) .LT. MODULO(dim2_len_total, proc_count)) THEN
    dim2_len_interior_list(i) = dim2_len_other + 1
    ! The indexing is weird here since processors start at index 0 but arrays
    ! start at index 1. Thus, entry i of this corresponds to processor (i-1).
  ELSE
    dim2_len_interior_list(i) = dim2_len_other
  END IF
END DO

dim2_global_index_low = SUM(dim2_len_interior_list(1:proc_id)) + 1
dim2_global_index_high = SUM(dim2_len_interior_list(1:proc_id+1))

! Each processor prints out its global dim2 indices.
DO i = 0, proc_count-1
  IF (proc_id .EQ. i) THEN
    PRINT *, 'proc_id: ', proc_id, ' global dim2 indices: ', &
    & dim2_global_index_low, dim2_global_index_high
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  ELSE
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  END IF
END DO
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

! Each processor prints out its matrix_interior.
!DO i = 0, proc_count-1
!  IF (proc_id .EQ. i) THEN
!    PRINT *, 'proc_id: ', proc_id, ' matrix_interior: '
!    DO j = 1, dim2_len_interior
!      PRINT *, matrix_interior(:,j)
!    END DO
!    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
!  ELSE
!    CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
!  END IF
!END DO
!CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)


END SUBROUTINE PARALLEL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ERROR_HANDLING. Checks for input error and other possible complications.
!
! ERRORS TO CHECK: - error1.
!
! VARIABLES: - proc_id: Processor ID (INTEGER).
! - proc_err_flag - Flag is triggered if processor encounters an error
! (LOGICAL).
! - global_err_flag - Flag is triggered if any processor encouters an error
! (LOGICAL).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE ERROR_HANDLING(proc_id)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: proc_id
LOGICAL :: proc_err_flag, &
& global_err_flag

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

proc_err_flag = .FALSE.
global_err_flag = .FALSE.

! error 1.
!IF (error1) THEN
!  proc_err_flag = .TRUE.
!  PRINT *, 'EZ_PARALLEL: Issue with processor ', proc_id, ' in subroutine ', &
!  & 'call SUBROUTINE_NAME. order: ', order, ' is less ', &
!  & 'than 1. Please increase the order.'
!END IF

! STOP all processors if any encounter an error.
CALL MPI_ALLREDUCE(proc_err_flag, global_err_flag, 1, MPI_LOGICAL, MPI_LOR, &
& MPI_COMM_WORLD, ierror)

IF (global_err_flag) THEN
  CALL MPI_FINALIZE(ierror)
  STOP
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

END SUBROUTINE ERROR_HANDLING
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE ZERO_PADDING_DBLE_REAL_EZP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
