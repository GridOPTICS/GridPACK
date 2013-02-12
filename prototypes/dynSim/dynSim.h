#ifndef __DYNSIM_H
#define __DYNSIM_H

#include <cstring>

#include <petsc.h>

#include "admitBuild.h"
#include "vars.h"

class DynSim {

public:

  explicit DynSim(const int ome);
  ~DynSim();

  PetscErrorCode readInputSizes();
  PetscErrorCode readInputData();

  PetscErrorCode allocMainData();
  PetscErrorCode deallocMainData();

  PetscErrorCode buildAdmittanceMatrix();

private:

  const int me;

  SizeGlobals sizeGlobals;
  BusGlobals busGlobals;
  LineGlobals lineGlobals;
  GenGlobals genGlobals;
  SwGlobals swGlobals;
  ReducedYBusGlobals reducedYBusGlobals;
  YBusGlobals yBusGlobals;
  SimuTime simuTime;
  SimuResults simuResults;

  DynSim() : me(-1) { };
  DynSim(const DynSim &othis) : me(-1) { };

  PetscErrorCode broadcastGlobals();
};

inline DynSim::DynSim(const int ome)
  : me(ome)
{
} // DynSim

inline DynSim::~DynSim()
{
  deallocMainData();
} // ~DynSim

inline PetscErrorCode DynSim::buildAdmittanceMatrix()
{
  PetscErrorCode ierr;

  AdmitBuild ab(me, sizeGlobals.nbus, sizeGlobals.nbrch, busGlobals.ns, busGlobals.bus_i,
		busGlobals.bus_gs, busGlobals.bus_bs, lineGlobals.line_from,
		lineGlobals.line_to, lineGlobals.line_r, lineGlobals.line_x,
		lineGlobals.line_charge, lineGlobals.line_ratio, lineGlobals.line_shift);

  ierr = ab.buildMatrix(yBusGlobals.ybus);
  CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
} // buildAdmittanceMatrix

#endif // __DYNSIM_H
