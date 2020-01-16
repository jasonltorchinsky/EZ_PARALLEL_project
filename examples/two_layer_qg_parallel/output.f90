!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Write output to .csv file.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE OUTPUT

IMPLICIT NONE

PRIVATE

! Defines standard integer-, real-precision types.
INCLUDE 'integer_types.h'
INCLUDE 'real_types.h'

PUBLIC :: WRITE_OUTPUT

CONTAINS


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE WRITE_OUTPUT(timestep, time, dt)
!  Writes the given array to a .csv file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

USE INITIALIZE

! ADDED TO PARALLEL
USE MPI

IMPLICIT NONE

INTEGER(qb) :: timestep, i, j
INTEGER(qb) :: proc_id ! ADDED TO PARALLEL
REAL(dp) :: time, dt
!~~ CHANGED FOR PARALLEL
CHARACTER(LEN=38) :: layer1_file_name, layer2_file_name
!CHARACTER(LEN=35) :: baroclin_file_name, barotrop_file_name
CHARACTER(LEN=43) :: timestep_info_file_name
!~~

! ADDED TO PARALLEL
CALL GET_ID(proc_id)

WRITE(layer1_file_name,"(A,I0.8,A,I0.4,A)") &
"./output_data/layer1_", timestep, "_", proc_id, ".csv"
WRITE(layer2_file_name,"(A,I0.8,A,I0.4,A)") &
"./output_data/layer2_", timestep, "_", proc_id, ".csv"
!WRITE(baroclin_file_name,"(A,I0.8,A)") "./output_data/baroclin_", timestep, ".csv"
!WRITE(barotrop_file_name,"(A,I0.8,A)") "./output_data/barotrop_", timestep, ".csv"
WRITE(timestep_info_file_name,"(A,I0.8,A,I0.4,A)") &
"./output_data/out_", timestep, "_", proc_id, "_info.txt"

OPEN(1001, file = layer1_file_name, form = "formatted")
OPEN(1002, file = layer2_file_name, form = "formatted")
!OPEN(1003, file = baroclin_file_name, form = "formatted")
!OPEN(1004, file = barotrop_file_name, form = "formatted")

DO j = 1 , y_len
  DO i = 1, x_len
    WRITE(1001, "(E32.16)", ADVANCE = 'NO') &
    REAL(physical_pot_vorticity_grid(i,j,1), dp)
    WRITE(1001, "(A)", ADVANCE = 'NO') ','

    WRITE(1002, "(E32.16)", ADVANCE = 'NO') &
    REAL(physical_pot_vorticity_grid(i,j,2), dp)
    WRITE(1002, "(A)", ADVANCE = 'NO') ','

    !WRITE(1003, "(E32.16)", ADVANCE = 'NO') &
    !REAL((0.5_dp, 0.0_dp) * (physical_pot_vorticity_grid(i,j,1) &
    !+ physical_pot_vorticity_grid(i,j,2)), dp)
    !WRITE(1003, "(A)", ADVANCE = 'NO') ','

    !WRITE(1004, "(E32.16)", ADVANCE = 'NO') &
    !REAL((0.5_dp, 0.0_dp) * (physical_pot_vorticity_grid(i,j,1) &
    !- physical_pot_vorticity_grid(i,j,2)), dp)
    !WRITE(1004, "(A)", ADVANCE = 'NO') ','
  END DO
  WRITE(1001, "(1x)")
  WRITE(1002, "(1x)")
  !WRITE(1003, "(1x)")
  !WRITE(1004, "(1x)")
END DO

CLOSE(1001)
CLOSE(1002)
!CLOSE(1003)
!CLOSE(1004)

!~~ CHANGED FOR PARALLEL
IF (proc_id .EQ. 0) THEN
  OPEN(1005, file = timestep_info_file_name, form = "formatted")
  WRITE(1005,"(A,I0.8,1x)") "step = ", timestep
  WRITE(1005,"(A,E32.16,1x)") "time = ", time
  WRITE(1005,"(A,E32.16,1x)") "dt = ", dt
  CLOSE(1005)
END IF
!~~

PRINT *, "Processor ", proc_id, " wrote outputs for step ", timestep, "."

END SUBROUTINE WRITE_OUTPUT

END MODULE
