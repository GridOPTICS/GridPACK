#include <iostream>
#include <fstream>
#include <string>

#include <petsc.h>

#include "dynSim.h"
#include "vars.h"

static const int skipSize = 96;
static const double sectionEndMarker = 9999999.0;

static void readAndCountOneDouble(std::ifstream &ifs, int &cnt);

static PetscErrorCode readBusVars(std::ifstream &ifs, BusGlobals &busGlobals, const int i);
static PetscErrorCode assembleBusVecs(BusGlobals &busGlobals);
static PetscErrorCode readLineVars(std::ifstream &ifs, LineGlobals &lineGlobals, const int i);
static PetscErrorCode assembleLineVecs(LineGlobals &lineGlobals);
static PetscErrorCode readGenVars(std::ifstream &ifs, GenGlobals &genGlobals, const int i);
static PetscErrorCode assembleGenVecs(GenGlobals &genGlobals);
static PetscErrorCode readSwVars(std::ifstream &ifs, SwGlobals &swGlobals, const int i);
static PetscErrorCode assembleSwVecs(SwGlobals &swGlobals);

PetscErrorCode DynSim::readInputSizes()
{
  PetscErrorCode ierr;
  PetscBool flg;
  char file[PETSC_MAX_PATH_LEN];

  ierr = PetscOptionsGetString(PETSC_NULL, "-i", file, PETSC_MAX_PATH_LEN, &flg);
  CHKERRQ(ierr);

  if (!flg)
    SETERRQ(PETSC_COMM_WORLD, 1, "Must indicate input file with the -i option");

  if (me == 0) { // Master process read only
    std::ifstream ifs(file);

    if (!ifs.good())
      SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Could not open input file: %s", file);

    sizeGlobals.nbus = 0;
    sizeGlobals.nbrch = 0;
    sizeGlobals.nSW = 0;
    sizeGlobals.nPV = 0;
    sizeGlobals.nPQ = 0;
    sizeGlobals.nswtch = 0;
    sizeGlobals.ngen = 0;

    bool flagf = false;
    char rbuf[skipSize];
    double rl1, rl2;

    while (!flagf) {
      int endSect, type;

      ifs >> rl1;

      endSect = static_cast<int>(rl1);
      
      if (endSect == sectionEndMarker)
	break;

      ifs.read(rbuf, skipSize);
      ifs >> rl2;

      type = static_cast<int>(rl2);

      if (endSect != sectionEndMarker) {
	sizeGlobals.nbus++;

	switch (type) {
	case swingBus:
	  sizeGlobals.nSW++;
	  break;
	case generatorBus:
	  sizeGlobals.nPV++;
	  break;
	case loadBus:
	  sizeGlobals.nPQ++;
	  break;
	}
      }
      else
	flagf = true;
    }

    readAndCountOneDouble(ifs, sizeGlobals.nbrch);
    readAndCountOneDouble(ifs, sizeGlobals.ngen);
    readAndCountOneDouble(ifs, sizeGlobals.nswtch);

    std::cout << "Number of buses: " << sizeGlobals.nbus << std::endl;
    std::cout << "Number of branches: " << sizeGlobals.nbrch << std::endl;
    std::cout << "Number of swing buses: " << sizeGlobals.nSW << std::endl;
    std::cout << "Number of PQ buses: " << sizeGlobals.nPQ << std::endl;
    std::cout << "Number of PV buses: " << sizeGlobals.nPV << std::endl;
    std::cout << "Number of generators: " << sizeGlobals.ngen << std::endl;
    std::cout << "Number of switches: " << sizeGlobals.nswtch << std::endl;
    std::cout << std::flush;
  }

  // Broadcast input sizes to all processes
  MPI_Bcast(&sizeGlobals, sizeof(sizeGlobals), MPI_BYTE, 0, PETSC_COMM_WORLD);

  PetscFunctionReturn(0);
} // readInputSizes

PetscErrorCode DynSim::readInputData()
{
  PetscErrorCode ierr;
  PetscBool flg;
  char file[PETSC_MAX_PATH_LEN];

  ierr = PetscOptionsGetString(PETSC_NULL, "-i", file, PETSC_MAX_PATH_LEN, &flg);
  CHKERRQ(ierr);

  if (!flg)
    SETERRQ(PETSC_COMM_WORLD, 1, "Must indicate input file with the -i option");

  if (me == 0) { // Master process read only
    std::ifstream ifs(file);

    if (!ifs.good())
      SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Could not open input file: %s", file);

    std::string str;

    int cnt = sizeGlobals.nbus;
    int i = 0;

    while (cnt > 0) {
      readBusVars(ifs, busGlobals, i);
      i++;
      cnt--;
    }
    assembleBusVecs(busGlobals);

    getline(ifs, str); // Skip previous line before end-of-data marker
    getline(ifs, str); // Skip line with end-of-data marker

    cnt = sizeGlobals.nbrch;
    i = 0;

    while (cnt > 0) {
      readLineVars(ifs, lineGlobals, i);
      i++;
      cnt--;
    }
    assembleLineVecs(lineGlobals);

    getline(ifs, str); // Skip previous line before end-of-data marker
    getline(ifs, str); // Skip line with end-of-data marker

    cnt = sizeGlobals.ngen;
    i = 0;

    while (cnt > 0) {
      readGenVars(ifs, genGlobals, i);
      i++;
      cnt--;
    }
    assembleGenVecs(genGlobals);

    getline(ifs, str); // Skip previous line before end-of-data marker
    getline(ifs, str); // Skip line with end-of-data marker

    cnt = sizeGlobals.nswtch;
    i = 0;

    while (cnt > 0) {
      readSwVars(ifs, swGlobals, i);
      i++;
      cnt--;
    }

    assembleSwVecs(swGlobals);

    // All data has been read
  }

  // Broadcast all read data from process 0 to all other processes
  broadcastGlobals();

  PetscFunctionReturn(0);
} // readInputData

void readAndCountOneDouble(std::ifstream &ifs, int &cnt)
{
  bool flagf = false;

  while (!flagf) {
    double rl1;
    std::string str;

    ifs >> rl1;
    getline(ifs, str);

    int endSect = static_cast<int>(rl1);

    if (!ifs.good())
      break;

    if (endSect != sectionEndMarker)
      cnt++;
    else
      flagf = true;
  }
} // readAndCountOneDouble

PetscErrorCode readBusVars(std::ifstream &ifs, BusGlobals &busGlobals, const int i)
{
  const int numBusVars = 8;

  double dbi;
  double dr[numBusVars];
  PetscScalar sc[numBusVars];
  PetscErrorCode ierr;

  ifs >> dbi;

  busGlobals.bus_i[i] = static_cast<int>(dbi) - 1; // Modify data to be zero-based

  for (int j = 0; j < numBusVars; j++) {
    ifs >> dr[j];
    sc[j] = dr[j];
  }

  ierr = VecSetValues(busGlobals.bus_v, 1, &i, &sc[0], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(busGlobals.bus_a, 1, &i, &sc[1], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(busGlobals.bus_pg, 1, &i, &sc[2], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(busGlobals.bus_qg, 1, &i, &sc[3], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(busGlobals.bus_pl, 1, &i, &sc[4], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(busGlobals.bus_ql, 1, &i, &sc[5], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(busGlobals.bus_gs, 1, &i, &sc[6], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(busGlobals.bus_bs, 1, &i, &sc[7], INSERT_VALUES);
  CHKERRQ(ierr);

  double dbt;

  ifs >> dbt;

  busGlobals.bus_type[i] = static_cast<int>(dbt);

  if (sc[6] != 0.0 || sc[7] != 0.0)
    busGlobals.ns++;

  PetscFunctionReturn(0);
} // readBusVars

PetscErrorCode assembleBusVecs(BusGlobals &busGlobals)
{
  PetscErrorCode ierr;

  ierr = VecAssemblyBegin(busGlobals.bus_v);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(busGlobals.bus_a);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(busGlobals.bus_pg);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(busGlobals.bus_qg);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(busGlobals.bus_pl);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(busGlobals.bus_ql);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(busGlobals.bus_gs);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(busGlobals.bus_bs);
  CHKERRQ(ierr);

  ierr = VecAssemblyEnd(busGlobals.bus_v);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(busGlobals.bus_a);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(busGlobals.bus_pg);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(busGlobals.bus_qg);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(busGlobals.bus_pl);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(busGlobals.bus_ql);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(busGlobals.bus_gs);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(busGlobals.bus_bs);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // assembleBusVecs

PetscErrorCode readLineVars(std::ifstream &ifs, LineGlobals &lineGlobals, const int i)
{
  const int numLineVars = 5;

  double dr[numLineVars];
  PetscScalar sc[numLineVars];
  PetscErrorCode ierr;

  double dlf, dlt;

  ifs >> dlf;
  ifs >> dlt;

  lineGlobals.line_from[i] = static_cast<int>(dlf) - 1; // Modify data to be zero-based
  lineGlobals.line_to[i] = static_cast<int>(dlt) - 1; // Modify data to be zero-based

  for (int j = 0; j < numLineVars; j++) {
    ifs >> dr[j];
    sc[j] = dr[j];
  }

  ierr = VecSetValues(lineGlobals.line_r, 1, &i, &sc[0], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(lineGlobals.line_x, 1, &i, &sc[1], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(lineGlobals.line_charge, 1, &i, &sc[2], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(lineGlobals.line_ratio, 1, &i, &sc[3], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(lineGlobals.line_shift, 1, &i, &sc[4], INSERT_VALUES);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // readLineVars

PetscErrorCode assembleLineVecs(LineGlobals &lineGlobals)
{
  PetscErrorCode ierr;

  ierr = VecAssemblyBegin(lineGlobals.line_r);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(lineGlobals.line_x);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(lineGlobals.line_charge);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(lineGlobals.line_ratio);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(lineGlobals.line_shift);
  CHKERRQ(ierr);

  ierr = VecAssemblyEnd(lineGlobals.line_r);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(lineGlobals.line_x);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(lineGlobals.line_charge);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(lineGlobals.line_ratio);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(lineGlobals.line_shift);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // assembleLineVecs

PetscErrorCode readGenVars(std::ifstream &ifs, GenGlobals &genGlobals, const int i)
{
  const int numGenVars = 17; // One variable in the data file is not used

  double dr[numGenVars];
  PetscScalar sc[numGenVars];
  PetscErrorCode ierr;

  double dgi, dgb;

  ifs >> dgi;
  ifs >> dgb;

  genGlobals.gen_i[i] = static_cast<int>(dgi) - 1;  // Modify data to be zero-based
  genGlobals.gen_bus[i] = static_cast<int>(dgb) - 1; // Modify data to be zero-based

  for (int j = 0; j < numGenVars; j++) {
    ifs >> dr[j];
    sc[j] = dr[j];
  }

  ierr = VecSetValues(genGlobals.gen_mva, 1, &i, &sc[0], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(genGlobals.gen_x, 1, &i, &sc[1], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(genGlobals.gen_r, 1, &i, &sc[2], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(genGlobals.gen_dsr, 1, &i, &sc[3], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(genGlobals.gen_dtr, 1, &i, &sc[4], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(genGlobals.gen_dstr, 1, &i, &sc[5], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(genGlobals.gen_dtc, 1, &i, &sc[6], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(genGlobals.gen_dstc, 1, &i, &sc[7], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(genGlobals.gen_qsr, 1, &i, &sc[8], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(genGlobals.gen_qtr, 1, &i, &sc[9], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(genGlobals.gen_qstr, 1, &i, &sc[10], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(genGlobals.gen_qtc, 1, &i, &sc[11], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(genGlobals.gen_qstc, 1, &i, &sc[12], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(genGlobals.gen_h, 1, &i, &sc[13], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(genGlobals.gen_d0, 1, &i, &sc[14], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(genGlobals.gen_d1, 1, &i, &sc[15], INSERT_VALUES);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // readGenVars

PetscErrorCode assembleGenVecs(GenGlobals &genGlobals)
{
  PetscErrorCode ierr;

  ierr = VecAssemblyBegin(genGlobals.gen_mva);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(genGlobals.gen_x);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(genGlobals.gen_r);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(genGlobals.gen_dsr);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(genGlobals.gen_dtr);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(genGlobals.gen_dstr);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(genGlobals.gen_dtc);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(genGlobals.gen_dstc);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(genGlobals.gen_qsr);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(genGlobals.gen_qtr);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(genGlobals.gen_qstr);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(genGlobals.gen_qtc);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(genGlobals.gen_qstc);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(genGlobals.gen_h);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(genGlobals.gen_d0);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(genGlobals.gen_d1);
  CHKERRQ(ierr);

  ierr = VecAssemblyEnd(genGlobals.gen_mva);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(genGlobals.gen_x);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(genGlobals.gen_r);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(genGlobals.gen_dsr);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(genGlobals.gen_dtr);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(genGlobals.gen_dstr);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(genGlobals.gen_dtc);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(genGlobals.gen_dstc);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(genGlobals.gen_qsr);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(genGlobals.gen_qtr);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(genGlobals.gen_qstr);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(genGlobals.gen_qtc);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(genGlobals.gen_qstc);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(genGlobals.gen_h);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(genGlobals.gen_d0);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(genGlobals.gen_d1);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // assembleGenVecs

PetscErrorCode readSwVars(std::ifstream &ifs, SwGlobals &swGlobals, const int i)
{
  const int numSwVars = 7;

  double dr[numSwVars];
  PetscScalar sc[numSwVars];
  PetscErrorCode ierr;

  for (int j = 0; j < numSwVars; j++) {
    ifs >> dr[j];
    sc[j] = dr[j];
  }

  ierr = VecSetValues(swGlobals.sw1, 1, &i, &sc[0], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(swGlobals.sw2, 1, &i, &sc[1], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(swGlobals.sw3, 1, &i, &sc[2], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(swGlobals.sw4, 1, &i, &sc[3], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(swGlobals.sw5, 1, &i, &sc[4], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(swGlobals.sw_type, 1, &i, &sc[5], INSERT_VALUES);
  CHKERRQ(ierr);
  ierr = VecSetValues(swGlobals.sw7, 1, &i, &sc[6], INSERT_VALUES);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // readSwVars

PetscErrorCode assembleSwVecs(SwGlobals &swGlobals)
{
  PetscErrorCode ierr;

  ierr = VecAssemblyBegin(swGlobals.sw1);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(swGlobals.sw2);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(swGlobals.sw3);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(swGlobals.sw4);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(swGlobals.sw5);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(swGlobals.sw_type);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin(swGlobals.sw7);
  CHKERRQ(ierr);

  ierr = VecAssemblyEnd(swGlobals.sw1);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(swGlobals.sw2);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(swGlobals.sw3);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(swGlobals.sw4);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(swGlobals.sw5);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(swGlobals.sw_type);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd(swGlobals.sw7);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // assembleSwVecs
