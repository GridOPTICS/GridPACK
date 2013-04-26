TASCEL_HOME=/pic/projects/mca/openmpi/build-intel

CFLAGS=-I$(TASCEL_HOME)/include
CXXFLAGS=-I$(TASCEL_HOME)/include
LDFLAGS=-L$(TASCEL_HOME)/lib -ltascel
FFLAGS = -g -IPF_fltacc  -I$(TASCEL_HOME)/include
CDEFS = -DAdd_

MPICXX = mpicxx
#EXTRA_LIBS += $(MP_LIBS) $(COMM_LIBS) $(LUBLAS) -lm
EXTRA_LIBS += $(LUBLAS) -limf -lm
#LIBS += $(EXTRA_LIBS)


FC = mpif90
CC = mpicc

#Super_DIR = /pic/people/d3m998/SuperLU_4.3
Super_DIR = /pic/projects/mca/intel/SuperLU_3.0
HEADER = $(Super_DIR)/SRC

LUBLAS = -L/share/apps/intel/11.1/mkl/lib/em64t -lmkl_lapack95_lp64 -lmkl_blas95_lp64\
-lmkl_intel_lp64 -lmkl_sequential -lmkl_solver_lp64_sequential -lmkl_core -lmkl_cdft_core\
-liomp5 -lpthread 

SuperLU = $(Super_DIR)/libsuperlu_3.0.a



OBJS = main.o counters.o mpimod.o runca.o allo.o readconfigfile.o readsize.o alloc.o readmodel.o chkinput.o readinput.o saveini.o \
       savepf.o printcaheader.o createfilenames.o loadini.o loadpf.o elio.o stat.o makeYbus.o \
       ext2int.o bustypes.o makegbus.o flatstart.o makeSbus.o c_fortran_dgssv.o newtonpf.o pfsoln.o \
       int2ext.o printpf.o dropline.o printca.o allocca.o limit.o \
       reportviolation.o writecaheader.o updateYbus.o updatecnt.o chkYbus.o 
#loadieee.o loadpti.o 
#createoutnames.o bbCA.o bgCA.o ggCA.o m2.o m3.o m1.o\
	 

#OBJ_CNT = counters.o

ALL_OBJ = $(OBJS) $(SuperLU) $(SP_LIB)
#ALL_OBJ = $(OBJS) $(OBJ_CNT) $(SuperLU) $(SP_LIB)

all: allo runca


#c_fortran_zgssv.o : c_fortran_zgssv.c
#	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

#counters.o: counters.c
#	$(CC) $(CFLAGS) $(CDEFS) -I$(ARMCI_INCLUDE) -c $< 



allo:   allo.o 
	$(FC) -c -g allo.f

runca:  $(ALL_OBJ) 
#	$(FC) $(FFLAGS) -o  runca $(ALL_OBJ) -libverbs -L$(ARMCI_LIB) -larmci $(EXTRA_LIBS) 
#	$(FC) $(FFLAGS) -o  runca $(ALL_OBJ) -L$(ARMCI_LIB) -libverbs  -larmci $(EXTRA_LIBS) 
	$(MPICXX) $(CDEFS) -o runca $(ALL_OBJ) $(LDFLAGS) -lmpi_f90 -lmpi_f77 -lmpi $(EXTRA_LIBS) 
#	$(MPICXX) $(CDEFS) -o runca $(ALL_OBJ) $(LDFLAGS) -lmpi_f90 -lmpi_f77 -lmpi 

clean:
	rm -f *~
	rm -f *.o
	rm -f *.mod
	rm -f runca
	rm -f *.chk

.SUFFIXES: .c .cc .f

#.c.o: 
#	$(MPICC) $(CFLAGS) $(CDEFS) -c $<
.c.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

.cc.o:
	$(MPICXX) $(CXXFLAGS) $(CDEFS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<
