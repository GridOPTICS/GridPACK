#ifndef __VARS_H
#define __VARS_H

#include <cmath>
#include <complex>

const int swingBus = 1;
const int generatorBus = 2;
const int loadBus = 3;

const double eps = 2.22e-10;
const std::complex<double> jay(0.0, 1.0);
const double basmva = 100.0;
const double sysFreq = 60.0;
const double basrad = 2.0 * M_PI * sysFreq;

// Struct declarations for DynSim members

struct SizeGlobals {
  int nbus, nbrch, nSW, nPV, nPQ, nswtch, ngen;

  SizeGlobals() { memset(this, 0, sizeof(*this)); };
};

struct BusGlobals {
  int ns, *bus_i;
  Vec bus_v, bus_a, bus_pg, bus_qg, bus_pl, bus_ql, bus_gs, bus_bs;
  int *bus_type;

  BusGlobals() { memset(this, 0, sizeof(*this)); };
};

struct LineGlobals {
  int *line_from, *line_to;
  Vec line_r, line_x, line_charge, line_ratio, line_shift;

  LineGlobals() { memset(this, 0, sizeof(*this)); };
};

struct GenGlobals {
  int *gen_i, *gen_bus;
  Vec gen_mva, gen_x, gen_r, gen_dsr, gen_dtr, gen_dstr, gen_dtc, gen_dstc;
  Vec gen_qsr, gen_qtr, gen_qstr, gen_qtc, gen_qstc, gen_h, gen_d0, gen_d1;

  GenGlobals() { memset(this, 0, sizeof(*this)); };
};

struct SwGlobals {
  Vec sw1, sw2, sw3, sw4, sw5, sw_type, sw7;

  SwGlobals() { memset(this, 0, sizeof(*this)); };
};

struct ReducedYBusGlobals {
  Mat prefy11, fy11, posfy11; // PETSc dense matrices

  ReducedYBusGlobals() { memset(this, 0, sizeof(*this)); };
};

struct YBusGlobals {
  Mat ybus, ybusSave; // PETSc sparse matrices

  YBusGlobals() { memset(this, 0, sizeof(*this)); };
};

struct SimuTime {
  Vec simuTime;

  SimuTime() { memset(this, 0, sizeof(*this)); };
};

struct SimuResults {
  Mat pelect, pmech, psi_re, psi_im, mac_ang, mac_spd; // PETSc dense
  Mat dmac_ang, dmac_spd, edprime, eqprime; // PETSc dense

  Vec eterm, qelect, curdg, curqg, curd, curq;
  Mat cur_re, cur_im; // PETSc dense
  Vec ed, eq, vex, theta;
  Vec bus_volt;

  PetscScalar t_isimuloop, t_macem;

  SimuResults() { memset(this, 0, sizeof(*this)); };    
};

#endif // __VARS_H
