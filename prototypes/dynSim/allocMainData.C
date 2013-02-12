#include <cstring>

#include <petsc.h>

#include "dynSim.h"
#include "vars.h"

PetscErrorCode DynSim::allocMainData()
{
  PetscErrorCode ierr;

  const int nbus = sizeGlobals.nbus;

  busGlobals.bus_i = new int[nbus];
  // Replicated vectors
  ierr = VecCreateSeq(PETSC_COMM_SELF, nbus, &busGlobals.bus_v);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nbus, &busGlobals.bus_a);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nbus, &busGlobals.bus_pg);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nbus, &busGlobals.bus_qg);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nbus, &busGlobals.bus_pl);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nbus, &busGlobals.bus_ql);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nbus, &busGlobals.bus_gs);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nbus, &busGlobals.bus_bs);
  CHKERRQ(ierr);
  busGlobals.bus_type = new int[nbus];

  const int nbrch = sizeGlobals.nbrch;

  lineGlobals.line_from = new int[nbrch];
  lineGlobals.line_to = new int[nbrch];
  ierr = VecCreateSeq(PETSC_COMM_SELF, nbrch, &lineGlobals.line_r);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nbrch, &lineGlobals.line_x);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nbrch, &lineGlobals.line_charge);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nbrch, &lineGlobals.line_ratio);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nbrch, &lineGlobals.line_shift);
  CHKERRQ(ierr);

  const int ngen = sizeGlobals.ngen;

  genGlobals.gen_i = new int[ngen];
  genGlobals.gen_bus = new int[ngen];
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_mva);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_x);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_r);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_dsr);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_dtr);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_dstr);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_dtc);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_dstc);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_qsr);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_qtr);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_qstr);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_qtc);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_qstc);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_h);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_d0);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, ngen, &genGlobals.gen_d1);
  CHKERRQ(ierr);

  const int nswtch = sizeGlobals.nswtch;

  ierr = VecCreateSeq(PETSC_COMM_SELF, nswtch, &swGlobals.sw1);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nswtch, &swGlobals.sw2);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nswtch, &swGlobals.sw3);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nswtch, &swGlobals.sw4);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nswtch, &swGlobals.sw5);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nswtch, &swGlobals.sw_type);
  CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, nswtch, &swGlobals.sw7);
  CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD, &reducedYBusGlobals.prefy11);
  CHKERRQ(ierr);
  ierr = MatSetType(reducedYBusGlobals.prefy11, MATMPIDENSE);
  CHKERRQ(ierr);
  ierr = MatSetSizes(reducedYBusGlobals.prefy11, PETSC_DECIDE, PETSC_DECIDE, ngen, ngen);
  CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD, &reducedYBusGlobals.fy11);
  CHKERRQ(ierr);
  ierr = MatSetType(reducedYBusGlobals.fy11, MATMPIDENSE);
  CHKERRQ(ierr);
  ierr = MatSetSizes(reducedYBusGlobals.fy11, PETSC_DECIDE, PETSC_DECIDE, ngen, ngen);
  CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD, &reducedYBusGlobals.posfy11);
  CHKERRQ(ierr);
  ierr = MatSetType(reducedYBusGlobals.posfy11, MATMPIDENSE);
  CHKERRQ(ierr);
  ierr = MatSetSizes(reducedYBusGlobals.posfy11, PETSC_DECIDE, PETSC_DECIDE, ngen, ngen);
  CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD, &yBusGlobals.ybus);
  CHKERRQ(ierr);
  ierr = MatSetType(yBusGlobals.ybus, MATMPIAIJ);
  CHKERRQ(ierr);
  ierr = MatSetSizes(yBusGlobals.ybus, PETSC_DECIDE, PETSC_DECIDE, nbus + 2, nbus + 2);
  CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD, &yBusGlobals.ybusSave);
  CHKERRQ(ierr);
  ierr = MatSetType(yBusGlobals.ybusSave, MATMPIAIJ);
  CHKERRQ(ierr);
  ierr = MatSetSizes(yBusGlobals.ybusSave, PETSC_DECIDE, PETSC_DECIDE, nbus + 2, nbus + 2);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // allocMainData

PetscErrorCode DynSim::deallocMainData()
{
  PetscErrorCode ierr;

  delete [] busGlobals.bus_i;
  ierr = VecDestroy(&busGlobals.bus_v);
  CHKERRQ(ierr);
  ierr = VecDestroy(&busGlobals.bus_a);
  CHKERRQ(ierr);
  ierr = VecDestroy(&busGlobals.bus_pg);
  CHKERRQ(ierr);
  ierr = VecDestroy(&busGlobals.bus_qg);
  CHKERRQ(ierr);
  ierr = VecDestroy(&busGlobals.bus_pl);
  CHKERRQ(ierr);
  ierr = VecDestroy(&busGlobals.bus_ql);
  CHKERRQ(ierr);
  ierr = VecDestroy(&busGlobals.bus_gs);
  CHKERRQ(ierr);
  ierr = VecDestroy(&busGlobals.bus_bs);
  CHKERRQ(ierr);
  delete [] busGlobals.bus_type;

  delete [] lineGlobals.line_from;
  delete [] lineGlobals.line_to;
  ierr = VecDestroy(&lineGlobals.line_r);
  CHKERRQ(ierr);
  ierr = VecDestroy(&lineGlobals.line_x);
  CHKERRQ(ierr);
  ierr = VecDestroy(&lineGlobals.line_charge);
  CHKERRQ(ierr);
  ierr = VecDestroy(&lineGlobals.line_ratio);
  CHKERRQ(ierr);
  ierr = VecDestroy(&lineGlobals.line_shift);
  CHKERRQ(ierr);

  delete [] genGlobals.gen_i;
  delete [] genGlobals.gen_bus;
  ierr = VecDestroy(&genGlobals.gen_mva);
  CHKERRQ(ierr);
  ierr = VecDestroy(&genGlobals.gen_x);
  CHKERRQ(ierr);
  ierr = VecDestroy(&genGlobals.gen_r);
  CHKERRQ(ierr);
  ierr = VecDestroy(&genGlobals.gen_dsr);
  CHKERRQ(ierr);
  ierr = VecDestroy(&genGlobals.gen_dtr);
  CHKERRQ(ierr);
  ierr = VecDestroy(&genGlobals.gen_dstr);
  CHKERRQ(ierr);
  ierr = VecDestroy(&genGlobals.gen_dtc);
  CHKERRQ(ierr);
  ierr = VecDestroy(&genGlobals.gen_dstc);
  CHKERRQ(ierr);
  ierr = VecDestroy(&genGlobals.gen_qsr);
  CHKERRQ(ierr);
  ierr = VecDestroy(&genGlobals.gen_qtr);
  CHKERRQ(ierr);
  ierr = VecDestroy(&genGlobals.gen_qstr);
  CHKERRQ(ierr);
  ierr = VecDestroy(&genGlobals.gen_qtc);
  CHKERRQ(ierr);
  ierr = VecDestroy(&genGlobals.gen_qstc);
  CHKERRQ(ierr);
  ierr = VecDestroy(&genGlobals.gen_h);
  CHKERRQ(ierr);
  ierr = VecDestroy(&genGlobals.gen_d0);
  CHKERRQ(ierr);
  ierr = VecDestroy(&genGlobals.gen_d1);
  CHKERRQ(ierr);

  ierr = VecDestroy(&swGlobals.sw1);
  CHKERRQ(ierr);
  ierr = VecDestroy(&swGlobals.sw2);
  CHKERRQ(ierr);
  ierr = VecDestroy(&swGlobals.sw3);
  CHKERRQ(ierr);
  ierr = VecDestroy(&swGlobals.sw4);
  CHKERRQ(ierr);
  ierr = VecDestroy(&swGlobals.sw5);
  CHKERRQ(ierr);
  ierr = VecDestroy(&swGlobals.sw_type);
  CHKERRQ(ierr);
  ierr = VecDestroy(&swGlobals.sw7);
  CHKERRQ(ierr);

  ierr = MatDestroy(&reducedYBusGlobals.prefy11);
  CHKERRQ(ierr);
  ierr = MatDestroy(&reducedYBusGlobals.fy11);
  CHKERRQ(ierr);
  ierr = MatDestroy(&reducedYBusGlobals.posfy11);
  CHKERRQ(ierr);

  ierr = MatDestroy(&yBusGlobals.ybus);
  CHKERRQ(ierr);
  ierr = MatDestroy(&yBusGlobals.ybusSave);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // deallocMainData
