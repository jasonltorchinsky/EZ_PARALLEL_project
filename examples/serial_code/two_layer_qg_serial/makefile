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

# Compiler.
FC = gfortran
# Compiler flags.
FCFLAGS = -O1
# Include file search path, e.g., -I/usr/include/
INPTHS = -I./EZDFFTPACK
# Library search path, e.g., -L./mpi_int_comm.
INLIBS = -L./EZDFFTPACK
# Link with specific libraries (.a or .so files), e.g. -lgfortran, -lnetcdff.
LKLIBS = -lEZDFFTPACK


TRGT = serial_two_layer_qg_eqn_solver.exe        # Executable name.
SRCS = $(wildcard *.f90)                        # Source files.
OBJS = $(patsubst %.f90, %.o, $(SRCS))           # Object files.
MODS = $(patsubst %.f90, %.mod, $(SRCS))         # Module files.

# Makefile.

# Makes the executable.
$(TRGT): $(OBJS)                # Make executable.
	$(FC) $(FCFLAGS) $(OBJS) -o $(TRGT) $(INLIBS) $(LKLIBS)

# Make main.o.
main.o: main.f90 time_stepper.o initialize.o
	$(FC) $(FCFLAGS) -c $< $(INLIBS) $(LKLIBS)

# Make time_stepper.o.
time_stepper.o: time_stepper.f90 output.o initialize.o jacobian_ekman_shear_solve.o
	$(FC) $(FCFLAGS) -c $<

# Make output.o.
output.o: output.f90 initialize.o
	$(FC) $(FCFLAGS) -c $<

# Make all files without dependencies.
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<


.PHONY : clean
clean:
#	rm -rf *.o *.mod
	del *.o *.mod *.exe