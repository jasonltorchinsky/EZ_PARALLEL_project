!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ZERO_PADDING_DBLE_CMPLX. Scales the input matrix by padding zeros (the
! 3/2-rule) for de-aliasing FFTs.
!
! ARGUMENTS: - dim1_len, dim2_len: The shape of the matrix (INTEGER).
! - matrix: The matrix to be scaled (DOUBLE COMPLEX,
! DIMENSION(dim1_len, dim2_len)).
! - scl_dim1_len, slc_dim2_len: The dimensions of scaled_matrix (INTEGER).
! - scaled_matrix: The resulting scaled matrix (DOUBLE COMPLEX,
! DIMENSION(scl_dim1_len, scl_dim2_len)).
! - overlap: The overlap required for the numerical scheme (INTEGER).
!
! NOTES: - When discussing the wavenumbers, we use the convention (for dim1)
! [0, ..., dim1_len/2, -dim1_len/1 + 1, ... -1], with FLOOR(dim1_len/2), as
! is assumed in Fortran integer division.
!
! Written By: Jason Turner
! Last Updated: January 28, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE ZERO_PADDING_DBLE_CMPLX_EZP(dim1_len, dim2_len, matrix, &
& scl_dim1_len, scl_dim2_len, scaled_matrix, overlap)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& overlap, &
& scl_dim1_len, &
& scl_dim2_len, &
& proc_id, &
& proc_count, &
& ierror
DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len) :: matrix
DOUBLE COMPLEX, DIMENSION(scl_dim1_len, scl_dim2_len) :: scaled_matrix

CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, proc_count, ierror)

! Check for errors in user input.
CALL ERROR_HANDLING(proc_id)

! Run in serial if there is only one processor.
IF (proc_count .EQ. 1) THEN
  CALL SERIAL(dim1_len, dim2_len, matrix, scl_dim1_len, scl_dim2_len, &
    & scaled_matrix)
ELSE
  CALL PARALLEL(dim1_len, dim2_len, matrix, scl_dim1_len, scl_dim2_len, &
    & scaled_matrix, overlap, proc_id, proc_count)
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
! - matrix: The matrix to be scaled (DOUBLE COMPLEX,
! DIMENSION(dim1_len, dim2_len)).
! - scl_dim1_len, slc_dim2_len: The dimensions of scaled_matrix (INTEGER).
! - scaled_matrix: The resulting scaled matrix (DOUBLE COMPLEX,
! DIMENSION(scl_dim1_len, scl_dim2_len)).
! - scaled_dim1_neg_wavenum_low_index, scaled_dim2_neg_wavenum_low_index: The
! lower index of the largest negative wavenumber in scaled_matrix (INTEGER).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE SERIAL(dim1_len, dim2_len, matrix, scl_dim1_len, scl_dim2_len, &
  & scaled_matrix)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len,  &
& scl_dim1_len, &
& scl_dim2_len, &
& scaled_dim1_neg_wavenum_low_index, &
& scaled_dim2_neg_wavenum_low_index
DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len) :: matrix
DOUBLE COMPLEX, DIMENSION(scl_dim1_len, scl_dim2_len) :: scaled_matrix

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

! This matrix is scale by 3/2 in each dimension, with the indices corresponding
! to the largest wavenumbers (in magnitude) zeroed out.
scaled_matrix(1:dim1_len/2+1, 1:dim2_len/2+1) = &
& matrix(1:dim1_len/2+1, 1:dim2_len/2+1)
scaled_matrix(1:dim1_len/2+1, &
  & scaled_dim2_neg_wavenum_low_index:scl_dim2_len) = &
& matrix(1:dim1_len/2+1, dim2_len/2+2:dim2_len)
scaled_matrix(scaled_dim1_neg_wavenum_low_index:scl_dim1_len, &
  & 1:dim2_len/2+1) = &
& matrix(dim1_len/2+2:dim1_len, 1:dim2_len/2+1)
scaled_matrix(scaled_dim1_neg_wavenum_low_index:scl_dim1_len, &
  & scaled_dim2_neg_wavenum_low_index:scl_dim2_len) = &
& matrix(dim1_len/2+2:dim1_len, dim2_len/2+2:dim2_len)

END SUBROUTINE SERIAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PARALLEL. Handles the execution if there is more than one processor.
!
! STRUCTURE: 1) Strip the overlap from each sub-grid to make them easier to
! work with. Scale each sub-grid interior along the first dimension.
! 2) Rearrange the partially scaled grid into slices along the first dimension
! (the dim2 slices).
! 3) Scale the dim2 slices along the second dimension.
! 4) Construct the interior of the scaled sub-grids using these scaled dim2
! slices.
! 4) Fill in the overlap region with the SHARE_SUBGRID_BOUNDARIES subroutine.
!
!
! VARIABLES: - dim1_len, dim2_len: The shape of the matrix (INTEGER).
! - matrix: The matrix to be scaled (DOUBLE COMPLEX,
! DIMENSION(dim1_len, dim2_len)).
! - scl_dim1_len, slc_dim2_len: The dimensions of scaled_matrix (INTEGER).
! - scaled_matrix: The resulting scaled matrix (DOUBLE COMPLEX,
! DIMENSION(scl_dim1_len, scl_dim2_len)).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - proc_id: Processor ID (INTEGER).
! - proc_count: Number of processors (INTEGER).
! - dim2_len_interior: dim2_len of sub-grid not counting the overlapping region
! (INTEGER).
! - i: Counting index for DO loops (INTEGER).
! - dim1_len_dim1_slice, dim2_len_dim1_slice, dim1_len_dim2_slice,
! dim2_len_dim2_slice: Dimensions of the dim1 and dim2 slices (INTEGER).
! - scl_dim1_len_dim1_slice, scl_dim2_len_dim1_slice, scl_dim1_len_dim2_slice,
! scl_dim2_len_dim2_slice: The dimensions of the scaled dim1 and dim2 slices
! (INTEGER).
! - scl_dim1_slice_pos_wvnm_low_ind, scl_dim1_slice_pos_wvnm_high_ind,
! - scl_dim1_slice_neg_wvnm_low_ind, scl_dim1_slice_neg_wvnm_high_ind,
! - scl_dim2_slice_pos_wvnm_low_ind, scl_dim2_slice_pos_wvnm_high_ind,
! - scl_dim2_slice_neg_wvnm_low_ind, scl_dim2_slice_neg_wvnm_high_ind: The
! indices for the positive/negative wavenumbers in the scaled dim1 and dim2
! slices (INTEGER).
! - fill_in_proc_trgt: The proc_id of the processor to communicate with during
! the dim2_slice fill-in process (INTEGER).
! - dim1_low_index_dim1_slice: Lesser global dim2 index for the portion of
! dim1_slice to be put into dim2_slice (INTEGER).
! - dim1_high_index_dim1_slice: Greater global dim2 index for the portion of
! dim1_slice to be put into dim2_slice (INTEGER).
! - dim2_low_index_dim2_slice: Lesser global dim2 index for the portion of
! dim2_slice to be filled in by dim1_slice (INTEGER).
! - dim2_high_index_dim2_slice: Greater global dim2 index for the portion of
! dim2_slice to be filled in by dim1_slice (INTEGER).
! - dim1_len_dim2_slice_list, dim2_len_dim1_slice_list: The list of dim1/dim2
! lengths for the dim1/dim2 slices (INTEGER, DIMENSION(proc_count)).
! - dim1_slice, dim2_slice, scl_dim1_slice, scl_dim2_slice: The slices of the
! global matrix and scaled matrix in each dimension (dim1 is a slice whose dim1
! length matches that of the global grid) (DOUBLE COMPLEX, DIMENSION(:,:),
! ALLOCATABLE).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE PARALLEL(dim1_len, dim2_len, matrix, scl_dim1_len, scl_dim2_len, &
  & scaled_matrix, overlap, proc_id, proc_count)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! FIX THIS UP FOR USE WITH INPUTTING THE SCALED MATIRX DIMENSIONS.
USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& scl_dim1_len, &
& scl_dim2_len, &
& overlap, &
& proc_id, &
& proc_count, &
& ierror, &
& dim1_len_dim1_slice, &
& dim2_len_dim1_slice, &
& dim1_len_dim2_slice, &
& dim2_len_dim2_slice, &
& scl_dim1_len_dim1_slice, &
& scl_dim2_len_dim1_slice, &
& scl_dim1_len_dim2_slice, &
& scl_dim2_len_dim2_slice, &
& scl_dim1_slice_pos_wvnm_low_ind, &
& scl_dim1_slice_pos_wvnm_high_ind, &
& scl_dim1_slice_neg_wvnm_low_ind, &
& scl_dim1_slice_neg_wvnm_high_ind, &
& scl_dim2_slice_pos_wvnm_low_ind, &
& scl_dim2_slice_pos_wvnm_high_ind, &
& scl_dim2_slice_neg_wvnm_low_ind, &
& scl_dim2_slice_neg_wvnm_high_ind, &
& fill_in_proc_trgt, &
& i, &
& j, &
& status(MPI_STATUS_SIZE), &
& dim1_low_index_dim1_slice, &
& dim1_high_index_dim1_slice, &
& dim2_low_index_dim2_slice, &
& dim2_high_index_dim2_slice
INTEGER, DIMENSION(proc_count) :: dim1_len_dim2_slice_list, &
& dim2_len_dim1_slice_list
DOUBLE COMPLEX, DIMENSION(dim1_len,dim2_len) :: matrix
DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: dim1_slice, &
& dim2_slice, &
& scl_dim1_slice, &
& scl_dim2_slice
DOUBLE COMPLEX, DIMENSION(scl_dim1_len, scl_dim2_len) :: scaled_matrix

scaled_matrix = 0.0
! Obtain the size of the dim1_slice (excluding overlap).
dim1_len_dim1_slice = dim1_len
! Middle processor.
IF (MOD(proc_id, proc_count-1) .NE. 0) THEN
  dim2_len_dim1_slice = dim2_len - 2 * overlap
! Processor 0 or proc_count-1.
ELSE
  dim2_len_dim1_slice = dim2_len - overlap
END IF
! Obtain the dim1_slice entries.
ALLOCATE(dim1_slice(dim1_len_dim1_slice, dim2_len_dim1_slice))
! Processor 0, doesn't need last parts of dim2.
IF (proc_id .EQ. 0) THEN
  dim1_slice = matrix(:,1:dim2_len-overlap)
! Last processor, doesn't need first part of dim2.
ELSE IF (proc_id .EQ. proc_count-1) THEN
  dim1_slice = matrix(:,overlap+1:dim2_len)
! Middle processor, doesn't need first or last part of dim2.
ELSE
  dim1_slice = matrix(:,overlap+1:dim2_len-overlap)
END IF

! Obtain the size of the scaled dim1 slice.
scl_dim1_len_dim1_slice = scl_dim1_len
scl_dim2_len_dim1_slice = dim2_len_dim1_slice

! Scale the dim1_slice.
ALLOCATE(scl_dim1_slice(scl_dim1_len_dim1_slice, scl_dim2_len_dim1_slice))
scl_dim1_slice = 0.0
! Determine the indices of scl_dim1_slice that dim1_slice entires go into.
scl_dim1_slice_pos_wvnm_low_ind = 1
scl_dim1_slice_pos_wvnm_high_ind = dim1_len_dim1_slice/2 + 1
scl_dim1_slice_neg_wvnm_high_ind = scl_dim1_len_dim1_slice
IF (MOD(dim1_len_dim1_slice, 2) .EQ. 0) THEN
  scl_dim1_slice_neg_wvnm_low_ind = dim1_len_dim1_slice + 2
ELSE IF (MOD(dim1_len_dim1_slice, 2) .EQ. 1) THEN
  scl_dim1_slice_neg_wvnm_low_ind = dim1_len_dim1_slice + 1
END IF

! Fill in the scaled dim1_slice.
scl_dim1_slice( &
  & scl_dim1_slice_pos_wvnm_low_ind:scl_dim1_slice_pos_wvnm_high_ind,:) = &
& dim1_slice(1:dim1_len_dim1_slice/2 + 1,:)
scl_dim1_slice( &
  & scl_dim1_slice_neg_wvnm_low_ind:scl_dim1_slice_neg_wvnm_high_ind,:) = &
& dim1_slice(dim1_len_dim1_slice/2 + 2:dim1_len_dim1_slice,:)

! At this point, we've scaled up in dim1. Now we must do the same for dim2.

! Obtain the size of the dim2 slice.
  ! Divide dim1_len of matrix as evenly as possible.
dim1_len_dim2_slice = scl_dim1_len/proc_count
  ! Now check if dim2_slice of a given processor needs to take additional row,
  ! i.e., dim1_len % proc_count != 0.
DO i = 1, proc_count
  IF ((i-1) .LT. MODULO(scl_dim1_len, proc_count)) THEN
    dim1_len_dim2_slice_list(i) = dim1_len_dim2_slice + 1
    ! The indexing is weird here since processors start at index 0 but arrays
    ! start at index 1. Thus, entry ii of this corresponds to processor (ii-1).
  ELSE
    dim1_len_dim2_slice_list(i) = dim1_len_dim2_slice
  END IF
END DO
dim1_len_dim2_slice = dim1_len_dim2_slice_list(proc_id + 1)
CALL MPI_ALLREDUCE(dim2_len_dim1_slice, dim2_len_dim2_slice, 1, MPI_INTEGER, &
& MPI_SUM, MPI_COMM_WORLD, ierror)
ALLOCATE(dim2_slice(dim1_len_dim2_slice, dim2_len_dim2_slice))

! For later use, obtain a list of all of the dim2_len_dim1_slice's.
  ! Divide dim2_len_dim2_slice of matrix as evenly as possible.
dim2_len_dim1_slice = dim2_len_dim2_slice/proc_count
  ! Now check if dim2_slice of a given processor needs to take additional row,
  ! i.e., dim1_len % proc_count != 0.
DO i = 1, proc_count
  IF ((i-1) .LT. MODULO(dim2_len_dim2_slice, proc_count)) THEN
    dim2_len_dim1_slice_list(i) = dim2_len_dim1_slice + 1
    ! The indexing is weird here since processors start at index 0 but arrays
    ! start at index 1. Thus, entry ii of this corresponds to processor (ii-1).
  ELSE
    dim2_len_dim1_slice_list(i) = dim2_len_dim1_slice
  END IF
END DO
dim2_len_dim1_slice = dim2_len_dim1_slice_list(proc_id + 1)

! Fill in the dim2_slice with sections of the scl_dim1_slice's.
DO i = 1, proc_count
  fill_in_proc_trgt = MOD(proc_count - proc_id + i, proc_count)
  dim1_high_index_dim1_slice = &
  & SUM(dim1_len_dim2_slice_list(1:fill_in_proc_trgt+1))
  dim1_low_index_dim1_slice = dim1_high_index_dim1_slice &
  & - dim1_len_dim2_slice_list(fill_in_proc_trgt+1) + 1
  dim2_high_index_dim2_slice = &
  & SUM(dim2_len_dim1_slice_list(1:fill_in_proc_trgt+1))
  dim2_low_index_dim2_slice = dim2_high_index_dim2_slice &
  & - dim2_len_dim1_slice_list(fill_in_proc_trgt+1) + 1

  ! If fill_in_proc_trgt < proc_id, RECV then SEND.
  IF (fill_in_proc_trgt .LT. proc_id) THEN
    CALL MPI_RECV( &
    & dim2_slice(:,dim2_low_index_dim2_slice:dim2_high_index_dim2_slice), &
    & dim1_len_dim2_slice_list(proc_id+1) &
    & * dim2_len_dim1_slice_list(fill_in_proc_trgt+1), MPI_DOUBLE_COMPLEX, &
    & fill_in_proc_trgt, proc_id, MPI_COMM_WORLD, status, ierror)

    CALL MPI_SEND( &
    & scl_dim1_slice(dim1_low_index_dim1_slice:dim1_high_index_dim1_slice,:), &
    & dim1_len_dim2_slice_list(fill_in_proc_trgt+1) &
    & * dim2_len_dim1_slice_list(proc_id+1), MPI_DOUBLE_COMPLEX, &
    & fill_in_proc_trgt, fill_in_proc_trgt, MPI_COMM_WORLD, ierror)

  ! If fill_in_proc_trgt > proc_id, SEND then RECV.
  ELSE IF (fill_in_proc_trgt .GT. proc_id) THEN
    CALL MPI_SEND( &
    & scl_dim1_slice(dim1_low_index_dim1_slice:dim1_high_index_dim1_slice,:), &
    & dim1_len_dim2_slice_list(fill_in_proc_trgt+1) &
    & * dim2_len_dim1_slice_list(proc_id+1), MPI_DOUBLE_COMPLEX, &
    & fill_in_proc_trgt, fill_in_proc_trgt, MPI_COMM_WORLD, ierror)

    CALL MPI_RECV( &
    & dim2_slice(:,dim2_low_index_dim2_slice:dim2_high_index_dim2_slice), &
    & dim1_len_dim2_slice_list(proc_id+1) &
    & * dim2_len_dim1_slice_list(fill_in_proc_trgt+1), MPI_DOUBLE_COMPLEX, &
    & fill_in_proc_trgt, proc_id, MPI_COMM_WORLD, status, ierror)

  ! If fill_in_proc_trgt = proc_id, then fill in dim2_slice from your
  ! dim1_slice.
  ELSE
    dim2_slice(:,dim2_low_index_dim2_slice:dim2_high_index_dim2_slice) = &
    & scl_dim1_slice(dim1_low_index_dim1_slice:dim1_high_index_dim1_slice,:)
  END IF
END DO

! Scale the dim2_slice.
scl_dim1_len_dim2_slice = dim1_len_dim2_slice
scl_dim2_len_dim2_slice = (3 * dim2_len_dim2_slice) / 2
ALLOCATE(scl_dim2_slice(scl_dim1_len_dim2_slice, scl_dim2_len_dim2_slice))
scl_dim2_slice = 0.0
! Determine the indices of scl_dim1_slice that dim1_slice entires go into.
scl_dim2_slice_pos_wvnm_low_ind = 1
scl_dim2_slice_pos_wvnm_high_ind = dim2_len_dim2_slice/2 + 1
scl_dim2_slice_neg_wvnm_high_ind = scl_dim2_len_dim2_slice
IF (MOD(dim2_len_dim2_slice, 2) .EQ. 0) THEN
  scl_dim2_slice_neg_wvnm_low_ind = dim2_len_dim2_slice + 2
ELSE IF (MOD(dim2_len_dim2_slice, 2) .EQ. 1) THEN
  scl_dim2_slice_neg_wvnm_low_ind = dim2_len_dim2_slice + 1
END IF

! Fill in the scaled dim2_slice.
scl_dim2_slice(:, &
  & scl_dim2_slice_pos_wvnm_low_ind:scl_dim2_slice_pos_wvnm_high_ind) = &
& dim2_slice(:,1:dim2_len_dim2_slice/2 + 1)
scl_dim2_slice(:, &
  & scl_dim2_slice_neg_wvnm_low_ind:scl_dim2_slice_neg_wvnm_high_ind) = &
& dim2_slice(:,dim2_len_dim2_slice/2 + 2:dim2_len_dim2_slice)

! At this point, we have scaled the sub-grids along dimension 2 as well.
! Reallocate the scl_dim1_slice to be scaled in both dimensions.
DEALLOCATE(scl_dim1_slice)
scl_dim1_len_dim1_slice = scl_dim1_len

! Calculate the new dim2_len_dim1_slice_list.
! Divide dim2_len_dim2_slice of matrix as evenly as possible.
scl_dim2_len_dim1_slice = scl_dim2_len_dim2_slice/proc_count
  ! Now check if dim2_slice of a given processor needs to take additional row,
  ! i.e., dim1_len % proc_count != 0.
DO i = 1, proc_count
  IF ((i-1) .LT. MODULO(scl_dim2_len_dim2_slice, proc_count)) THEN
    dim2_len_dim1_slice_list(i) = scl_dim2_len_dim1_slice + 1
    ! The indexing is weird here since processors start at index 0 but arrays
    ! start at index 1. Thus, entry ii of this corresponds to processor (ii-1).
  ELSE
    dim2_len_dim1_slice_list(i) = scl_dim2_len_dim1_slice
  END IF
END DO
scl_dim2_len_dim1_slice = dim2_len_dim1_slice_list(proc_id + 1)
ALLOCATE(scl_dim1_slice(scl_dim1_len_dim1_slice, scl_dim2_len_dim1_slice))

! Fill in the scl_dim1_slice with sections of the scl_dim2_slice's.
DO i = 1, proc_count
  fill_in_proc_trgt = MOD(proc_count - proc_id + i, proc_count)
  dim1_high_index_dim1_slice = &
  & SUM(dim1_len_dim2_slice_list(1:fill_in_proc_trgt+1))
  dim1_low_index_dim1_slice = dim1_high_index_dim1_slice &
  & - dim1_len_dim2_slice_list(fill_in_proc_trgt+1) + 1
  dim2_high_index_dim2_slice = &
  & SUM(dim2_len_dim1_slice_list(1:fill_in_proc_trgt+1))
  dim2_low_index_dim2_slice = dim2_high_index_dim2_slice &
  & - dim2_len_dim1_slice_list(fill_in_proc_trgt+1) + 1

  ! If fill_in_proc_trgt < proc_id, RECV then SEND.
  IF (fill_in_proc_trgt .LT. proc_id) THEN
    CALL MPI_RECV( &
    & scl_dim1_slice(dim1_low_index_dim1_slice:dim1_high_index_dim1_slice,:), &
    & dim1_len_dim2_slice_list(fill_in_proc_trgt+1) &
    & * dim2_len_dim1_slice_list(proc_id+1), MPI_DOUBLE_COMPLEX, &
    & fill_in_proc_trgt, proc_id, MPI_COMM_WORLD, status, ierror)

    CALL MPI_SEND( &
    & scl_dim2_slice(:,dim2_low_index_dim2_slice:dim2_high_index_dim2_slice), &
    & dim1_len_dim2_slice_list(proc_id+1) &
    & * dim2_len_dim1_slice_list(fill_in_proc_trgt+1), MPI_DOUBLE_COMPLEX, &
    & fill_in_proc_trgt, fill_in_proc_trgt, MPI_COMM_WORLD, ierror)

  ! If fill_in_proc_trgt > proc_id, SEND then RECV.
  ELSE IF (fill_in_proc_trgt .GT. proc_id) THEN
    CALL MPI_SEND( &
    & scl_dim2_slice(:,dim2_low_index_dim2_slice:dim2_high_index_dim2_slice), &
    & dim1_len_dim2_slice_list(proc_id+1) &
    & * dim2_len_dim1_slice_list(fill_in_proc_trgt+1), MPI_DOUBLE_COMPLEX, &
    & fill_in_proc_trgt, fill_in_proc_trgt, MPI_COMM_WORLD, ierror)

    CALL MPI_RECV( &
    & scl_dim1_slice(dim1_low_index_dim1_slice:dim1_high_index_dim1_slice,:), &
    & dim1_len_dim2_slice_list(fill_in_proc_trgt+1) &
    & * dim2_len_dim1_slice_list(proc_id+1), MPI_DOUBLE_COMPLEX, &
    & fill_in_proc_trgt, proc_id, MPI_COMM_WORLD, status, ierror)

  ! If fill_in_proc_trgt = proc_id, then fill in dim2_slice from your
  ! dim1_slice.
  ELSE
    scl_dim1_slice(dim1_low_index_dim1_slice:dim1_high_index_dim1_slice,:) = &
    & scl_dim2_slice(:,dim2_low_index_dim2_slice:dim2_high_index_dim2_slice)
  END IF
END DO

! Fill in interior of scaled_matrix.
IF (proc_id .EQ. 0) THEN
  ! First processor, fill in all except last part along dimension 2 of matrix.
  scaled_matrix(:, 1:scl_dim2_len-overlap) = scl_dim1_slice
ELSE IF (proc_id .EQ. proc_count-1) THEN
  ! Last processor, fill in all but first part along dimension 2 of matrix.
  scaled_matrix(:, 1+overlap:scl_dim2_len) = scl_dim1_slice
ELSE
  ! Middle processor, fill in all but first and last part along dimension 2 of
  ! matrix.
  scaled_matrix(:, 1+overlap:scl_dim2_len-overlap) = scl_dim1_slice
END IF


CALL SHARE_SUBGRID_BOUNDARIES_DBLE_CMPLX_EZP(scl_dim1_len, scl_dim2_len, &
  & overlap, scaled_matrix)

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

END SUBROUTINE ZERO_PADDING_DBLE_CMPLX_EZP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
