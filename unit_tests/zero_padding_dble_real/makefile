#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Path to EZ_PARALLEL library.
EZPPTH = ../../EZ_PARALLEL

# Command-Line System Options: Linux, Windows.
CMD = Linux
# FC options: mpifort, gfortran (with msmpi).
FC = mpifort

# Executable name.
TRGT = zero_padding_dble_real_unit_test.exe
# .f files.
SRCS = $(wildcard *.f)
# .o files.
OBJS = $(patsubst %.f, %.o, $(SRCS))
# .mod files.
MODS = $(patsubst %.f, %.mod, $(SRCS))

# Set the compiler flags, dependent on the compiler.
ifeq ($(FC), gfortran)
  # Path to msmpi include.
  MPIPTHINC = C:/lib/mpi/include
  # Path to msmpi lib.
  MPIPTHLIB = C:/lib/mpi/lib 
  # Compiler flags.
  FCFLAGS = -ffree-form
  # Include file search path, e.g., C:/lib/mpi/include.
  INCPTHS = $(EZPPTH) $(MPIPTHINC)
  # Library search path, e.g., C:/lib/mpi/lib.
  INCLIBS = $(EZPPTH) $(MPIPTHLIB) 
  # Link with specific libraries (.a or .so files), e.g. EZ_PARALLEL.
  LIBNAMES = EZ_PARALLEL msmpi
  
else ifeq ($(FC), mpifort)
  # Compiler flags.
  FCFLAGS = -ffree-form
  # Include file search paths, e.g., /usr/include/.
  INCPTHS = $(EZPPTH)
  # Library search path, e.g., ../../EZ_PARALLEL.
  INCLIBS = $(EZPPTH) 
  # Link with specific libraries (.a or .so files), e.g. EZ_PARALLEL.
  LIBNAMES = EZ_PARALLEL
endif

# Add appropriate prefixes to include file search paths, library search paths,
# and include libraries.
INCFLAGS = $(addprefix -I,$(INCPTHS))
INCLIBFLAGS = $(addprefix -L,$(INCLIBS))
LKLIBFLAGS = $(addprefix -l,$(LIBNAMES))

# Makes the executable.
$(TRGT): $(OBJS)                # Make executable.
	$(FC) $< -o $(TRGT) $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)

# Compile all .f files. 
%.o: %.f
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)


.PHONY : clean tidy
# Clean deletes everything besides source code.
clean:
ifeq ($(CMD), Windows)
	del *.o *.mod *.exe
else ifeq ($(CMD), Linux)
	rm -rf *.o *.mod *.exe
endif

# Tidy deletes all .o and .mod files, keeping the executable.
tidy:
ifeq ($(CMD), Windows)
	del *.o *.mod
else ifeq ($(CMD), Linux)
	rm -rf *.o *.mod
endif
