EXEC   = sage 

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
			./code/model_cooling_heating.o \
			./code/model_starformation_and_feedback.o \
			./code/model_disk_instability.o \
			./code/model_reincorporation.o \
			./code/model_mergers.o \
			./code/model_misc.o

INCL   =	./code/core_allvars.h  \
			./code/core_proto.h  \
			./code/core_simulation.h  \
			./Makefile


#OPT += -DNOUT=36	      # This sets the number of galaxy output times

OPT += -DMILLENNIUM         # Millennium simulation trees
# OPT += -DBOLSHOI            # Bolshoi simulation trees
# OPT += -DGIGGLEZ						# GiggleZ simulation trees

# OPT += -DMINIMIZE_IO        # tree files will be preloaded, galaxy data will be written in one go


#SYSTYPE = "mac"
#SYSTYPE = "green"
SYSTYPE = "gstar"
#SYSTYPE = "IMAC"
#SYSTYPE = "Mojtaba_Mac"

CC       =   mpicc            # sets the C-compiler (default)
OPTIMIZE =   -pg -O2 -Wall    # optimization and warning flags (default)
GSL_INCL =   -I $(GSL_DIR)/include
GSL_LIBS =   -L $(GSL_DIR)/lib

ifeq ($(SYSTYPE),"mac")
CC       =  mpicc
GSL_INCL = -I/opt/local/include
GSL_LIBS = -L/opt/local/lib
endif

ifeq ($(SYSTYPE),"green")
CC       = /usr/local/gnu/x86_64/openmpi-1.4/bin/mpicc
OPTIMIZE = -O3 -Wall
GSL_INCL = -I/usr/local/gnu/x86_64/gsl/include
GSL_LIBS = -L/usr/local/gnu/x86_64/gsl/lib
endif

ifeq ($(SYSTYPE),"gstar")
CC       = /usr/local/x86_64/gnu/openmpi-1.4.5/bin/mpicc
OPTIMIZE = -O3 -Wall
GSL_INCL = -I/usr/local/x86_64/gnu/gsl-1.9/include
GSL_LIBS = -L/usr/local/x86_64/gnu/gsl-1.9/lib
endif

# USE-MPI = yes  # set this if you want to run in parallel

ifdef USE-MPI
OPT += -DMPI  #  This creates an MPI version that can be used to process files in parallel
CC = mpicc  # sets the C-compiler
else
CC = cc  # sets the C-compiler
endif

ifeq ($(SYSTYPE),"IMAC")
OPTIMIZE = -g -O0 -Wall  # optimization and warning flags
GSL_INCL = -I/usr/local/include  # make sure your system knows where GSL_DIR is
GSL_LIBS = -L/usr/local/lib
endif

ifeq ($(SYSTYPE),"Mojtaba_Mac")
OPTIMIZE = -g -O0 -Wall  # optimization and warning flags
GSL_INCL = -I/Users/mojtabaraouf/Ureka/python/include
GSL_LIBS = -L/Users/mojtabaraouf/Ureka/python/lib -lgsl -lm
endif




LIBS   =   -pg -lm  $(GSL_LIBS) -lgsl -lgslcblas

CFLAGS =   -pg $(OPTIONS) $(OPT) $(OPTIMIZE) $(GSL_INCL)


default: all

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)

$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS)

tidy:
	rm -f $(OBJS) ./$(EXEC)

all:  tidy $(EXEC) clean
