#include <petsc.h>

#include "dynSim.h"
#include "vars.h"

static PetscErrorCode broadcastVector(const Vec &v);
static PetscErrorCode broadcastVector(const int *v, const int sz);

PetscErrorCode DynSim::broadcastGlobals()
{
  broadcastVector(&busGlobals.ns, 1);
  broadcastVector(busGlobals.bus_i, sizeGlobals.nbus);
  broadcastVector(busGlobals.bus_v);
  broadcastVector(busGlobals.bus_a);
  broadcastVector(busGlobals.bus_pg);
  broadcastVector(busGlobals.bus_qg);
  broadcastVector(busGlobals.bus_pl);
  broadcastVector(busGlobals.bus_ql);
  broadcastVector(busGlobals.bus_gs);
  broadcastVector(busGlobals.bus_bs);

  broadcastVector(lineGlobals.line_from, sizeGlobals.nbrch);
  broadcastVector(lineGlobals.line_to, sizeGlobals.nbrch);
  broadcastVector(lineGlobals.line_r);
  broadcastVector(lineGlobals.line_x);
  broadcastVector(lineGlobals.line_charge);
  broadcastVector(lineGlobals.line_ratio);
  broadcastVector(lineGlobals.line_shift);

  broadcastVector(genGlobals.gen_i, sizeGlobals.ngen);
  broadcastVector(genGlobals.gen_bus, sizeGlobals.ngen);
  broadcastVector(genGlobals.gen_mva);
  broadcastVector(genGlobals.gen_x);
  broadcastVector(genGlobals.gen_r);
  broadcastVector(genGlobals.gen_dsr);
  broadcastVector(genGlobals.gen_dtr);
  broadcastVector(genGlobals.gen_dstr);
  broadcastVector(genGlobals.gen_dtc);
  broadcastVector(genGlobals.gen_dstc);
  broadcastVector(genGlobals.gen_qsr);
  broadcastVector(genGlobals.gen_qtr);
  broadcastVector(genGlobals.gen_qstr);
  broadcastVector(genGlobals.gen_qtc);
  broadcastVector(genGlobals.gen_qstc);
  broadcastVector(genGlobals.gen_h);
  broadcastVector(genGlobals.gen_d0);
  broadcastVector(genGlobals.gen_d1);

  broadcastVector(swGlobals.sw1);
  broadcastVector(swGlobals.sw2);
  broadcastVector(swGlobals.sw3);
  broadcastVector(swGlobals.sw4);
  broadcastVector(swGlobals.sw5);
  broadcastVector(swGlobals.sw_type);
  broadcastVector(swGlobals.sw7);

  PetscFunctionReturn(0);
} // broadcastGlobals

PetscErrorCode broadcastVector(const Vec &v)
{
  PetscErrorCode ierr;
  PetscScalar *dat;
  PetscInt sz;

  ierr = VecGetSize(v, &sz);
  CHKERRQ(ierr);
  ierr = VecGetArray(v, &dat);
  CHKERRQ(ierr);

  MPI_Bcast(dat, sz, MPIU_SCALAR, 0, PETSC_COMM_WORLD);

  ierr = VecRestoreArray(v, &dat);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // broadcastVector

PetscErrorCode broadcastVector(const int *v, const int sz)
{
  MPI_Bcast(const_cast<int *>(v), sz, MPI_INT, 0, PETSC_COMM_WORLD);

  PetscFunctionReturn(0);
} // broadcastVector
