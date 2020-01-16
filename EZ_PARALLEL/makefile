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
# Last Edited - January 10th, 2020
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LIB = EZ_PARALLEL

# Command-Line System options: Linux, Windows.
CMD = Windows
# FC options: mpifort, gfortran (with msmpi).
FC = gfortran

# Set the compiler flags, dependent on the compiler.
ifeq ($(FC), gfortran)
  # Path to msmpi include.
  MPIPTHINC = C:\lib\mpi\include
  # Path to msmpi lib.
  MPIPTHLIB = C:\lib\mpi\lib 
  # Compiler flags.
  FCFLAGS = -ffree-form
  # Include file search path, e.g., C:\lib\mpi\include.
  INCPTHS = $(MPIPTHINC)
  # Library search path, e.g., C:\lib\mpi\lib.
  INCLIBS = $(MPIPTHLIB) 
  # Link with specific libraries (.a or .so files), e.g. EZ_PARALLEL.
  LIBNAMES = msmpi
  
else ifeq ($(FC), mpifort)
  # Compiler flags.
  FCFLAGS = -ffree-form
  # Include file search paths, e.g., \usr\include\.
  INCPTHS = 
  # Library search path, e.g., ..\..\EZ_PARALLEL.
  INCLIBS =  
  # Link with specific libraries (.a or .so files), e.g. EZ_PARALLEL.
  LIBNAMES = 
endif

# Add appropriate prefixes to include file search paths, library search paths,
# and include libraries.
INCFLAGS = $(addprefix -I,$(INCPTHS))
INCLIBFLAGS = $(addprefix -L,$(INCLIBS))
LKLIBFLAGS = $(addprefix -l,$(LIBNAMES))

# EZ_PARALLEL file names.
EZPARNMS = \
init_mpi\
fin_mpi\
decomp_grid\
identify_ref_point\
share_subgrid_boundaries_dble_real\
share_subgrid_boundaries_dble_cmplx\
get_id\
cfft2df\
cfft2db\
rfft2df\
rfft2db\
maxval\
minval\
spectral_dim1_derivative\
spectral_dim2_derivative

# EZ_PARALLEL source file names.
EZPARSRCS = $(addsuffix .f,$(basename $(EZPARNMS)))
# EZ_PARALLEL object file names.
EZPAROBJS = $(addsuffix .o,$(basename $(EZPARNMS)))

# dfftpack file names.
DFFTPACKNMS = \
zfftb\
cfftb1\
zfftf\
cfftf1\
zffti\
cffti1\
dcosqb\
cosqb1\
dcosqf\
cosqf1\
dcosqi\
dcost\
dcosti\
ezfft1\
dzfftb\
dzfftf\
dzffti\
passb\
passb2\
passb3\
passb4\
passb5\
passf\
passf2\
passf3\
passf4\
passf5\
radb2\
radb3\
radb4\
radb5\
radbg\
radf2\
radf3\
radf4\
radf5\
radfg\
dfftb\
rfftb1\
dfftf\
rfftf1\
dffti\
rffti1\
dsinqb\
dsinqf\
dsinqi\
dsint\
sint1\
dsinti

# DFFTPACK source file names.
DFFTPACKSRCS = $(addsuffix .f,$(basename $(DFFTPACKNMS)))
# EZ_PARALLEL object file names.
DFFTPACKPTH = ./dfftpack/
DFFTPACKOBJS = $(addprefix $(DFFTPACKPTH), \
$(addsuffix .o,$(basename $(DFFTPACKNMS))))

# Compilation rules.
# Archive all ,.o files to library.
lib$(LIB).a: $(EZPAROBJS) $(DFFTPACKOBJS)
	ar crv lib$(LIB).a $(EZPAROBJS) $(DFFTPACKOBJS)

# Compile all DFFTPACK source code.
DFFTPACK/%.o: %.f
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)

# Compile all EZ_PARALLEL source code
%.o: %.f
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
	