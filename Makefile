EXEC   = darksage

OBJS   = 	./code/main.o \
			./code/core_read_parameter_file.o \
			./code/core_init.o \
			./code/core_io_tree.o \
			./code/core_cool_func.o \
			./code/core_build_model.o \
			./code/core_save.o \
			./code/core_mymalloc.o \
			./code/core_allvars.o \
			./code/model_infall.o \
			./code/model_cooling.o \
			./code/model_starformation_and_feedback.o \
			./code/model_disk_instability.o \
			./code/model_reincorporation.o \
			./code/model_mergers.o \
			./code/model_misc.o \

INCL   =	./code/core_allvars.h  \
			./code/core_proto.h  \
			./code/core_simulation.h  \
			./Makefile


#OPT += -DMPI
#OPT += -DMINIMIZE_IO

#SYSTYPE = "mac"
# SYSTYPE = "green"
# SYSTYPE = "gstar"

ifndef CC
CC       =   mpicc            # sets the C-compiler (default)

#ifdef USE-MPI
#OPT += -DMPI  #  This creates an MPI version that can be used to process files in parallel
#CC = mpicc  # sets the C-compiler
#else
#CC = cc  # sets the C-compiler
#endif

endif

OPTIMIZE =   -g -O2 -Wall -Wextra -Wshadow   # optimization and warning flags (default)

# GSL automatic detection
GSL_FOUND := $(shell gsl-config --version 2>/dev/null)
ifndef GSL_FOUND
$(warning GSL not found in path - please install GSL before installing SAGE (or, update the PATH environment variable such that "gsl-config" is found))
# if the automatic detection fails, set GSL_DIR appropriately
GSL_DIR := /opt/local
GSL_INCL := -I$(GSL_DIR)/include
GSL_LIBDIR := $(GSL_DIR)/lib
# since GSL is not in PATH, the runtime environment might not be setup correctly either
# therefore, adding the compiletime library path is even more important (the -Xlinker bit)
GSL_LIBS := -L$(GSL_LIBDIR) -lgsl -lgslcblas -Xlinker -rpath -Xlinker $(GSL_LIBDIR)
else
# GSL is probably configured correctly, pick up the locations automatically
GSL_INCL := $(shell gsl-config --cflags)
GSL_LIBDIR := $(shell gsl-config --prefix)/lib
GSL_LIBS   := $(shell gsl-config --libs) -Xlinker -rpath -Xlinker $(GSL_LIBDIR)
endif



#ifeq ($(SYSTYPE),"mac")
#CC       =  mpicc
#GSL_INCL = -I/opt/local/include
#GSL_LIBS = -L/opt/local/lib
#endif

#ifeq ($(SYSTYPE),"green")
#CC       = /usr/local/gnu/x86_64/openmpi-1.4/bin/mpicc
#OPTIMIZE = -O3 -Wall
#GSL_INCL = -I/usr/local/gnu/x86_64/gsl/include
#GSL_LIBS = -L/usr/local/gnu/x86_64/gsl/lib
#endif

#ifeq ($(SYSTYPE),"gstar")
#CC       = /usr/local/x86_64/gnu/openmpi-1.4.5/bin/mpicc
#OPTIMIZE = -O3 -Wall
#GSL_INCL = -I/usr/local/x86_64/gnu/gsl-1.9/include
#GSL_LIBS = -L/usr/local/x86_64/gnu/gsl-1.9/lib
#endif


LIBS   =   -g -lm  $(GSL_LIBS) -lgsl -lgslcblas

CFLAGS =   -g $(OPTIONS) $(OPT) $(OPTIMIZE) $(GSL_INCL)

default: all

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC) -g

$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS)

tidy:
	rm -f $(OBJS) ./$(EXEC)

all:  tidy $(EXEC) clean
