!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! CFFT2DF. Returns the Fourier Transform (not normalized) of a 2D DOUBLE COMPLEX
! matrix.
!
! ARGUMENTS: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - matrix: The 2D DOUBLE COMPLEX array which will be Fourier transformed
! (DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len)).
!
! NOTES: - This subroutine is based on the FFT2DC subroutine included in the
! UCLA LES t12 code, and the name Bjorn Stevens is included as the author of
! the documentation.
! - This subroutine is built on the FFTPACK Fortran library, developed
! by Paul N. Schwartztrauber (SW). The DOUBLE PRECISION version, used here,
! was converted from the original by Hugh C. Pumphrey, and is available at
!      https://www.netlib.org/fftpack/dp.tgz
!
! Written By: Jason Turner
! Last Updated: January 23, 2020
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE CFFT2DF_EZP(dim1_len, dim2_len, overlap, matrix)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& overlap, &
& proc_id, &
& proc_count, &
& ierror
DOUBLE COMPLEX :: matrix(dim1_len * dim2_len)

CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierror)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, proc_count, ierror)

! Check for errors in user input.
CALL ERROR_HANDLING(proc_id)

! Run in serial if there is only one processor.
IF (proc_count .EQ. 1) THEN
  CALL SERIAL(dim1_len, dim2_len, matrix)
ELSE
  CALL PARALLEL(dim1_len, dim2_len, overlap, matrix, proc_id, proc_count)
END IF

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SERIAL. Handles the execution if there is only one processor.
!
! STRUCTURE: 1) Fills WSAVE1, and transforms the matrix along dimension 1
! (dim1).
! 2) Fills WSAVE2, transposes the matrix, and transforms the transpose along
! dim1 (which is dimension 2 of the matrix). Returns the transpose of the
! tranformed transpose matrix.
!
! VARIABLES: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - matrix: The 2D DOUBLE COMPLEX array which will be Fourier transformed
! (DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len)).
! - WSAVE1: Work array to store prime factorization of dim1_len and tabulate
! trig functions for usage in SW's FFT routines (DOUBLE PRECISION,
! DIMENSION(4*dim1_len+15)).
! - WSAVE2: Work array to store prime factorization of dim2_len and tabulate
! trig functions for usage in SW's FFT routines (DOUBLE PRECISION,
! DIMENSION(4*dim2_len+15)).
! - matrix_tr: The 2D DOUBLE COMPLEX array which will store the transpose of
! matrix during the Fourier transform (DOUBLE COMPLEX,
! DIMENSION(dim2_len, dim1_len)).
! - ii: Counting index used in DO loops (INTEGER).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE SERIAL(dim1_len, dim2_len, matrix)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& ii
DOUBLE PRECISION, DIMENSION(4*dim1_len+15) :: WSAVE1
DOUBLE PRECISION, DIMENSION(4*dim2_len+15) :: WSAVE2
DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len) :: matrix
DOUBLE COMPLEX, DIMENSION(dim2_len, dim1_len) :: matrix_tr

! SW complex forward 1D FFT (FFFT) initialization, fills WSAVE1.
WSAVE1 = 0.0
CALL ZFFTI(dim1_len, WSAVE1)

! Perform 1D FFFT on the dim1 of matrix.
DO ii = 1, dim2_len
  CALL ZFFTF(dim1_len, matrix(:,ii), WSAVE1) ! SW FFFT routine.
END DO
matrix_tr = TRANSPOSE(matrix)

WSAVE2 = 0.0
CALL ZFFTI(dim2_len, WSAVE2)
! SW 1D (FFFT) initialization, fills WSAVE2.

! Perform 1D FFFT on the dim1 of matrix_transpose.
DO ii = 1, dim1_len
  CALL ZFFTF(dim2_len, matrix_tr(:,ii), WSAVE2)
END DO
matrix = TRANSPOSE(matrix_tr)

END SUBROUTINE SERIAL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PARALLEL. Handles the execution if there is more than one processor.
!
! STRUCTURE: 1) Trim the overlap from the input matrix. FFT input matrix along
! dim1.
! 2) Get the sub-grids that the processor would have if the grid was split along
! dimension 2 (dim2). FFT those sub-grids along dim2.
! 3) Send dim2 sub-grid sections to their appropriate places in the
! dim1 sub-grids
!
! VARIABLES: - dim1_len, dim2_len: The shape of the sub-grid (INTEGER).
! - overlap: The overlap required for the numerical scheme (INTEGER).
! - matrix: The 2D DOUBLE COMPLEX array which will be Fourier transformed
! (DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len)).
! - proc_id: Processor ID (INTEGER).
! - proc_count: Total number of processors (INTEGER).
! - WSAVE1: Work array to store prime factorization of dim1_len and tabulate
! trig functions for usage in SW's FFT routines (DOUBLE PRECISION,
! DIMENSION(4*dim1_len+15)).
! - WSAVE2: Work array to store prime factorization of dim2_len and tabulate
! trig functions for usage in SW's FFT routines (DOUBLE PRECISION,
! DIMENSION(:), ALLOCATABLE).
! - dim1_slice: The sub-grid (along dim1) that does not include its overlap with
! neighboring subgrids (DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE).
! - dim2_slice: The sub-grid (along dim2) that does not include its overlap with
! neighboring subgrids (DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE).
! - ii: Counting index used in DO loops (INTEGER).
! - dim2_len_dim1_slice: dim2 length of dim1_slice, later recalculated for use
! in communication (INTEGER).
! - dim1_len_dim2_slice: dim1 length of dim2_slice, later recalculated for use
! in communication (INTEGER).
! - dim1_len_dim2_slice_list: List of all dim1_len for all dim2_slice's
! (INTEGER, DIMENSION(proc_count)).
! - dim2_len_dim2_slice: dim2 length of dim2_slice, later recalculated for use
! in communication (INTEGER).
! - fill_in_proc_trgt: The proc_id of the processor to communicate with during
! the dim2_slice fill-in process (INTEGER).
! - dim2_len_total: dim2 length of the grid (INTEGER).
! - dim2_len_dim1_slice_list: List of all dim2_len for all dim1_slice's
! (INTEGER, DIMENSION(proc_count)).
! - dim1_low_index_dim1_slice: Lesser global dim2 index for the portion of
! dim1_slice to be put into dim2_slice (INTEGER).
! - dim1_high_index_dim1_slice: Greater global dim2 index for the portion of
! dim1_slice to be put into dim2_slice (INTEGER).
! - dim2_low_index_dim2_slice: Lesser global dim2 index for the portion of
! dim2_slice to be filled in by dim1_slice (INTEGER).
! - dim2_high_index_dim2_slice: Greater global dim2 index for the portion of
! dim2_slice to be filled in by dim1_slice (INTEGER).
! - dim2_slice_tr: The 2D DOUBLE COMPLEX array which will store the transpose of
! dim2_slice during the Fourier transform (DOUBLE COMPLEX, DIMENSION(:,:),
! ALLOCATABLE).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE PARALLEL(dim1_len, dim2_len, overlap, matrix, proc_id, proc_count)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE MPI

IMPLICIT NONE

INTEGER :: dim1_len, &
& dim2_len, &
& overlap, &
& proc_id, &
& proc_count, &
& ierror, &
& status(MPI_STATUS_SIZE), &
& ii, &
& dim2_len_dim1_slice, &
& dim1_len_dim2_slice, &
& dim2_len_dim2_slice, &
& fill_in_proc_trgt, &
& dim2_len_total, &
& dim1_low_index_dim1_slice, &
& dim1_high_index_dim1_slice, &
& dim2_low_index_dim2_slice, &
& dim2_high_index_dim2_slice
INTEGER, DIMENSION(proc_count) :: dim1_len_dim2_slice_list, &
& dim2_len_dim1_slice_list
DOUBLE PRECISION, DIMENSION(4*dim1_len+15) :: WSAVE1
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WSAVE2
DOUBLE COMPLEX, DIMENSION(dim1_len, dim2_len) :: matrix
DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE :: dim1_slice, &
& dim2_slice, &
& dim2_slice_tr


! Trim the overlap from the input matrix.
IF (proc_id .EQ. 0) THEN
  ! First processor, trim last part along dimension 2 of matrix.
  dim2_len_dim1_slice = dim2_len - overlap
  ALLOCATE(dim1_slice(dim1_len, dim2_len_dim1_slice))
  dim1_slice = matrix(:, 1:dim2_len-overlap)
ELSE IF (proc_id .EQ. proc_count-1) THEN
  ! Last processor, trim first part along dimension 2 of matrix.
  dim2_len_dim1_slice = dim2_len - overlap
  ALLOCATE(dim1_slice(dim1_len, dim2_len_dim1_slice))
  dim1_slice = matrix(:, 1+overlap:dim2_len)
ELSE
  ! Middle processor, trim first and last part along dimension 2 of matrix.
  dim2_len_dim1_slice = dim2_len - 2 * overlap
  ALLOCATE(dim1_slice(dim1_len, dim2_len_dim1_slice))
  dim1_slice = matrix(:, 1+overlap:dim2_len-overlap)
END IF

! Perfom 1D FFT along first dimension of dim1_slice.
  ! SW complex forward 1D FFT (FFFT) initialization, fills WSAVE1.
WSAVE1 = 0.
CALL ZFFTI(dim1_len, WSAVE1)

DO ii = 1, dim2_len_dim1_slice
  CALL ZFFTF(dim1_len, dim1_slice(:,ii), WSAVE1) ! SW FFFT routine.
END DO

! ~~~ Up to this point, we have FFT'd the grid along dimension 1.

! Now need to FFT along dim2, so we need a dim2_slice. First, we must
! calculate the size of dim2_slice for all processors.
! Calculate dim1_len for each dim2_slice.
  ! Divide dim1_len of total grid as evenly as possible.
dim1_len_dim2_slice = dim1_len/proc_count
  ! Now check if dim2_slice of a given processor needs to take additional row,
  ! i.e., dim1_len % proc_count != 0.
DO ii = 1, proc_count
  IF ((ii-1) .LT. MODULO(dim1_len, proc_count)) THEN
    dim1_len_dim2_slice_list(ii) = dim1_len_dim2_slice + 1
    ! The indexing is weird here since processors start at index 0 but arrays
    ! start at index 1. Thus, entry ii of this corresponds to processor (ii-1).
  ELSE
    dim1_len_dim2_slice_list(ii) = dim1_len_dim2_slice
  END IF
END DO


! Calculate dim2_len for all dim2_slice's by adding up dim2 for each processor's
! dim1_slice.
CALL MPI_ALLREDUCE(dim2_len_dim1_slice, dim2_len_dim2_slice, 1, MPI_INTEGER, &
& MPI_SUM, MPI_COMM_WORLD, ierror)


! Allocate memory to dim2_slice and fill it in.
ALLOCATE(dim2_slice(dim1_len_dim2_slice_list(proc_id + 1), dim2_len_dim2_slice))
dim2_slice = (0.0, 0.0)

! To fill in dim2_slice, we need the dim2_len of all dim1_slice's.
CALL MPI_ALLREDUCE(dim2_len_dim1_slice, dim2_len_total, 1, MPI_INTEGER, &
& MPI_SUM, MPI_COMM_WORLD, ierror)
! Calculate dim2_len for each dim1_slice, to avoid communication between
! processors.
  ! Divide dim2_len of total grid as evenly as possible.
dim2_len_dim1_slice = dim2_len_total/proc_count
  ! Now check if dim1_slice of a given processor needs to take additional row,
  ! i.e., dim2_len % proc_count != 0.
DO ii = 1, proc_count
  IF ((ii-1) .LT. MODULO(dim2_len_total, proc_count)) THEN
    dim2_len_dim1_slice_list(ii) = dim2_len_dim1_slice + 1
    ! The indexing is weird here since processors start at index 0 but arrays
    ! start at index 1. Thus, entry ii of this corresponds to processor (ii-1).
  ELSE
    dim2_len_dim1_slice_list(ii) = dim2_len_dim1_slice
  END IF
END DO

! Fill in dim2_slice using sections of the dim1_slice's.
DO ii = 1, proc_count
  fill_in_proc_trgt = MOD(proc_count - proc_id + ii, proc_count)
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
    & dim1_len_dim2_slice_list(proc_id+1)&
    & * dim2_len_dim1_slice_list(fill_in_proc_trgt+1), MPI_DOUBLE_COMPLEX, &
    & fill_in_proc_trgt, proc_id, MPI_COMM_WORLD, status, ierror)

    CALL MPI_SEND( &
    & dim1_slice(dim1_low_index_dim1_slice:dim1_high_index_dim1_slice,:), &
    & dim1_len_dim2_slice_list(fill_in_proc_trgt+1) &
    & * dim2_len_dim1_slice_list(proc_id+1), MPI_DOUBLE_COMPLEX, &
    & fill_in_proc_trgt, fill_in_proc_trgt, MPI_COMM_WORLD, ierror)

  ! If fill_in_proc_trgt > proc_id, SEND then RECV.
  ELSE IF (fill_in_proc_trgt .GT. proc_id) THEN
    CALL MPI_SEND( &
    & dim1_slice(dim1_low_index_dim1_slice:dim1_high_index_dim1_slice,:), &
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
    & dim1_slice(dim1_low_index_dim1_slice:dim1_high_index_dim1_slice,:)
  END IF
END DO

! ~~~~~~~~~~~~~~~ Up to this point, each processor has a filled-in dim2_slice,
! and we are ready to FFT along dim2.
ALLOCATE( &
& dim2_slice_tr(dim2_len_dim2_slice, dim1_len_dim2_slice_list(proc_id + 1)))
dim2_slice_tr = (0.0, 0.0)
dim2_slice_tr = TRANSPOSE(dim2_slice)
ALLOCATE(WSAVE2(4*dim2_len_dim2_slice+15))

! Perfom 1D FFT along first dimension of dim2_slice_tr (the columns of the work
! grid).
  ! SW complex forward 1D FFT (FFFT) initialization, fills WSAVE1.
WSAVE2 = 0.
CALL ZFFTI(dim2_len_dim2_slice, WSAVE2)

DO ii = 1, dim1_len_dim2_slice_list(proc_id + 1)
  ! SW FFFT routine.
  CALL ZFFTF(dim2_len_dim2_slice, dim2_slice_tr(:,ii), WSAVE2)
END DO

dim2_slice = TRANSPOSE(dim2_slice_tr)

! ~~~ Up to this point, we have FFT'd along dim2, and now we must reconstruct
! the dim1_slice's.

! Fill in dim1_slice using sections of the dim2_slice's.
DO ii = 1, proc_count
  fill_in_proc_trgt = MOD(proc_count - proc_id + ii, proc_count)
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
    & dim1_slice(dim1_low_index_dim1_slice:dim1_high_index_dim1_slice,:), &
    & dim1_len_dim2_slice_list(fill_in_proc_trgt+1) &
    & * dim2_len_dim1_slice_list(proc_id+1), MPI_DOUBLE_COMPLEX, &
    & fill_in_proc_trgt, proc_id, MPI_COMM_WORLD, status, ierror)

    CALL MPI_SEND( &
    & dim2_slice(:,dim2_low_index_dim2_slice:dim2_high_index_dim2_slice), &
    & dim1_len_dim2_slice_list(proc_id+1) &
    & * dim2_len_dim1_slice_list(fill_in_proc_trgt+1), MPI_DOUBLE_COMPLEX, &
    & fill_in_proc_trgt, fill_in_proc_trgt, MPI_COMM_WORLD, ierror)

  ! If fill_in_proc_trgt > proc_id, SEND then RECV.
  ELSE IF (fill_in_proc_trgt .GT. proc_id) THEN
    CALL MPI_SEND( &
    & dim2_slice(:,dim2_low_index_dim2_slice:dim2_high_index_dim2_slice), &
    & dim1_len_dim2_slice_list(proc_id+1) &
    & * dim2_len_dim1_slice_list(fill_in_proc_trgt+1), MPI_DOUBLE_COMPLEX, &
    & fill_in_proc_trgt, fill_in_proc_trgt, MPI_COMM_WORLD, ierror)

    CALL MPI_RECV( &
    & dim1_slice(dim1_low_index_dim1_slice:dim1_high_index_dim1_slice,:), &
    & dim1_len_dim2_slice_list(fill_in_proc_trgt+1) &
    & * dim2_len_dim1_slice_list(proc_id+1), MPI_DOUBLE_COMPLEX, &
    & fill_in_proc_trgt, proc_id, MPI_COMM_WORLD, status, ierror)

  ! If fill_in_proc_trgt = proc_id, then fill in dim2_slice from your
  ! dim1_slice.
  ELSE
    dim1_slice(dim1_low_index_dim1_slice:dim1_high_index_dim1_slice,:) = &
    & dim2_slice(:,dim2_low_index_dim2_slice:dim2_high_index_dim2_slice)
  END IF
END DO

! Fill in interior of input matrix.
IF (proc_id .EQ. 0) THEN
  ! First processor, fill in all except last part along dimension 2 of matrix.
  matrix(:, 1:dim2_len-overlap) = dim1_slice
ELSE IF (proc_id .EQ. proc_count-1) THEN
  ! Last processor, fill in all but first part along dimension 2 of matrix.
  matrix(:, 1+overlap:dim2_len) = dim1_slice
ELSE
  ! Middle processor, fill in all but first and last part along dimension 2 of
  ! matrix.
  matrix(:, 1+overlap:dim2_len-overlap) = dim1_slice
END IF

CALL SHARE_SUBGRID_BOUNDARIES_DBLE_CMPLX_EZP(dim1_len, dim2_len, overlap, &
& matrix)

END SUBROUTINE PARALLEL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ERROR_HANDLING. Checks for input error and other possible complications.
!
! ERRORS TO CHECK: N/A. (Assume that users will not enter negative matrix
! dimension lengths.)
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

END SUBROUTINE CFFT2DF_EZP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
