#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is the makefile for compiling the EZ_PARALLEL library. It contains two 
# sets of options: Command-Line System options and Fortran Compiler (FC) 
# options. The Command-Line System options determine the 'make clean' and 
# 'make tidy' commands that will be executed, (e.g., 'rm' for Linux and 'del'
# for Windows. The FC options determine which compiler will be used, and
# it is assumed that msmpi will be used with gfortran. File paths used in 
# this makefile are either standard installation paths.
#
# After running 'make', the archive file libEZ_PARALLEL.a produced can be moved
# to other directories to be used with other programs.
#
# Created by - Jason Turner
# Last Edited - January 24th, 2020
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LIB = EZ_PARALLEL

# Command-Line System options: Linux, Windows.
CMD = Linux
# FC options: mpifort, gfortran (with msmpi).
FC = mpifort

# Set the compiler flags, dependent on the compiler.
ifeq ($(FC), gfortran)
  # Path to msmpi include.
  MPIPTHINC = C:/lib/mpi/include
  # Path to msmpi lib.
  MPIPTHLIB = C:/lib/mpi/lib 
  # Compiler flags.
  FCFLAGS = -ffree-form
  # Include file search path, e.g., C:/lib/mpi/include.
  INCPTHS = $(MPIPTHINC)
  # Library search path, e.g., C:/lib/mpi/lib.
  INCLIBS = $(MPIPTHLIB) 
  # Link with specific libraries (.a or .so files), e.g. EZ_PARALLEL.
  LIBNAMES = msmpi
  
else ifeq ($(FC), mpifort)
  # Compiler flags.
  FCFLAGS = -ffree-form
  # Include file search paths, e.g., /usr/include/.
  INCPTHS = 
  # Library search path, e.g., ../../EZ_PARALLEL.
  INCLIBS =  
  # Link with specific libraries (.a or .so files), e.g. EZ_PARALLEL.
  LIBNAMES = 
endif

# Add appropriate prefixes to include file search paths, library search paths,
# and include libraries.
INCFLAGS = $(addprefix -I,$(INCPTHS))
INCLIBFLAGS = $(addprefix -L,$(INCLIBS))
LKLIBFLAGS = $(addprefix -l,$(LIBNAMES))

# EZ_PARALLEL source file names.
EZPARSRCS = $(wildcard *.f)
# EZ_PARALLEL object file names.
EZPAROBJS = $(patsubst %.f, %.o, $(EZPARSRCS))

# DFFTPACK source file names.
DFFTPACKSRCS = $(wildcard ./DFFTPACK/*.f)
# EZ_PARALLEL object file names.
DFFTPACKOBJS = $(patsubst %.f, %.o, $(DFFTPACKSRCS))

# The build directory to hold all .o files.
OBJS = $(EZPAROBJS) $(DFFTPACKOBJS)

# Compilation rules.
# Archive all .o files to library.
lib$(LIB).a: $(OBJS)
	ar crv lib$(LIB).a $(OBJS)

# Compile all DFFTPACK source code.
%.o: DFFTPACK/%.f
	$(FC) -c $< 

# Compile all EZ_PARALLEL source code
%.o: ./%.f
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)
	
.PHONY : clean tidy
# Clean deletes everything besides source code.
clean:
ifeq ($(CMD), Windows)
	del *.o *.a
else ifeq ($(CMD), Linux)
	rm -rf *.o *.a
endif

# Tidy deletes all .o files, keeping the .a file.
tidy:
ifeq ($(CMD), Windows)
	del *.o
else ifeq ($(CMD), Linux)
	rm -rf *.o
endif
	