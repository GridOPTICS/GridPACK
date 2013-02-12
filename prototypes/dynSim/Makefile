TARGET = dynSim
OBJS = allocMainData.o dynSim.o readInput.o bcastGlobals.o admitBuild.o
MPICC = mpicc
MPICXX = mpicxx
CFLAGS = -03
CFLAGS = 

include $(PETSC_DIR)/conf/petscvariables
#include ${PETSC_DIR}/conf/rules

all: $(TARGET)

dynSim: $(OBJS)
	$(MPICXX) -o $(TARGET) $(OBJS) ${PETSC_MAT_LIB}

%.o: %.C
	$(MPICXX) ${PETSC_CC_INCLUDES} $(CFLAGS) -c -o $@ $<

.PHONY: clean

clean:
	rm -f *~ $(OBJS) $(TARGET)
