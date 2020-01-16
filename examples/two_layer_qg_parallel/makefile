#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    This is a makefile that will give instructions on how to link the other
#    files in this folder.
#    Please see http://genius2k.is-programmer.com/posts/40301.html for a 
#    tutorial on makefiles. In short, each stanza is of the form
#
#    TARGET: PREREQUISITES
#    	COMMAND
#
#    The -c commands compile without linking every source file. The -o commands
#    link all of the object files.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Objects, modules, and compiler.

# Options for Linux, Mac with mpifort installed.
# Compiler.
#FC = mpifort
# Compiler flags.
#FCFLAGS =
# Include file search path, e.g., -I/usr/include/
#INPTHS = -I./EZ_PARALLEL
# Library search path, e.g., -L./mpi_int_comm.
#INLIBS = -L./EZ_PARALLEL
# Link with specific libraries (.a or .so files), e.g. -lgfortran, -lnetcdff.
#LKLIBS = -lEZ_PARALLEL

# Options for Windows with msmpi installed.
# Compiler.
FC = gfortran
# Compiler flags.
FCFLAGS = -O1
# Include file search path, e.g., -I\usr\include\.
INCLDPTHS = .\EZ_PARALLEL C:\lib\mpi\include
INPTHS = $(addprefix -I,$(INCLDPTHS))
# Library search path, e.g., -L.\mpi_int_comm.
INCLDLIBS = .\EZ_PARALLEL C:\lib\mpi\lib 
INLIBS = $(addprefix -L,$(INCLDLIBS))
# Link with specific libraries (.a or .so files), e.g. -lgfortran, -lnetcdff.
LIBNAMES = EZ_PARALLEL msmpi
LKLIBS = $(addprefix -l,$(LIBNAMES))

TRGT = parallel_two_layer_qg_eqn_solver.exe        # Executable name.
SRCS = $(wildcard *.f90)                        # Source files.
OBJS = $(patsubst %.f90, %.o, $(SRCS))           # Object files.
MODS = $(patsubst %.f90, %.mod, $(SRCS))         # Module files.

# Makefile.

# Makes the executable.
$(TRGT): $(OBJS)                # Make executable.
	$(FC) $(FCFLAGS) $(OBJS) -o $(TRGT) $(INPTHS) $(INLIBS) $(LKLIBS)

# Make main.o.
main.o: main.f90 time_stepper.o initialize.o
	$(FC) $(FCFLAGS) -c $< $(INPTHS) $(INLIBS) $(LKLIBS)

# Make time_stepper.o.
time_stepper.o: time_stepper.f90 output.o initialize.o jacobian_ekman_shear_solve.o
	$(FC) $(FCFLAGS) -c $< $(INPTHS) $(INLIBS) $(LKLIBS)

# Make output.o.
output.o: output.f90 initialize.o
	$(FC) $(FCFLAGS) -c $< $(INPTHS) $(INLIBS) $(LKLIBS)

# Make all files without dependencies.
%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(INPTHS) $(INLIBS) $(LKLIBS)


.PHONY : clean
clean:
# Linux
#	rm -rf *.o *.mod *.exe
# Windows
	del *.mod *.o *.exe

.PHONY : tidy
tidy:
# Linux
#	rm -rf *.o *.mod
# Windows
	del *.mod *.o