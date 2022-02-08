/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   equivalent_matrix.hpp
 * @author Bruce Palmer
 * @Last modified:   March 24, 2022
 * 
 * @brief
 * This routine is designed to construct a reduced matrix from the values of a
 * collection of solution vectors. Each vector is added to the matrix one column
 * (or row) at a time. The entire matrix can then be accessed from any process.
 * 
 * 
 */

#include "gridpack/applications/modules/dynamic_simulation_full_y/equivalent_matrix.hpp"
#include "gridpack/mapper/full_map.hpp"

namespace gridpack {
namespace dynamic_simulation {
/**
 * Basic constructor
 * @param[in] config pointer to configuration object
 * @param[in] network pointer to dynamic simulation network. Network should
 *            already be distributed and initialized.
 * @param[in] flag: 1 for generating negative sequence equivalent matrix
              flag: 2 for generating zero sequence equivalent matrix
 */
EquivalentMatrix::EquivalentMatrix(gridpack::utility::Configuration *config,
    boost::shared_ptr<DSFullNetwork> &network, int flag) : p_initialized(false)
{
  gridpack::utility::StringUtils util;
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.Dynamic_simulation.Equivalent_matrix");
  // Read in buses from input deck
  std::string eq_buses;
  cursor->get("equivalentBusSet", &eq_buses);
  util.trim(eq_buses);
  if (eq_buses.length() > 0) {
    std::vector<int> buses;
    gridpack::utility::StringUtils util;
    std::vector<std::string> str_buses = util.blankTokenizer(eq_buses);
    int i;
    for (i=0; i<str_buses.size(); i++) {
      buses.push_back(atoi(str_buses[i].c_str()));
    }

    p_dim = buses.size();
    if (p_dim <= 0) {
      std::cout<<"Illegal number of buses used to create EquivalentMatrix: "
        <<p_dim<<std::endl;
    }
    for (i=0; i<p_dim; i++) {
      p_map.insert(std::pair<int,int>(buses[i],i));
    }
    p_network = network;

    // Create global array representing equivalent matrix
    int dims[2] = {p_dim,p_dim};
    p_ga = NGA_Create_handle();
    NGA_Set_pgroup(p_ga, network->communicator().getGroup());
    NGA_Set_data(p_ga,2,dims,C_DCPL);
    NGA_Allocate(p_ga);
    NGA_Zero(p_ga);
    p_initialized = true;

    // Add values to global array. Create Y-matrix
    boost::shared_ptr<DSFullFactory> factory;
    factory.reset(new DSFullFactory(network));
    if (flag==1){
      factory->setMode(Ybus_neg);
    } else if (flag==2){
      factory->setMode(Ybus_zero);
    } else {
      factory->setMode(YBUS);
    }
	
    boost::shared_ptr < gridpack::mapper::FullMatrixMap<DSFullNetwork> > ybusMap;
    ybusMap.reset(new gridpack::mapper::FullMatrixMap<DSFullNetwork> (p_network));
    boost::shared_ptr<gridpack::math::Matrix> Ybus;
    Ybus = ybusMap->mapToMatrix();

    // Initialize first right hand side vector and create solution vector
    factory->setMode(EQBUS);
    p_ivecMap.reset(new gridpack::mapper::BusVectorMap<DSFullNetwork>(p_network));
    boost::shared_ptr<gridpack::math::Vector> iVec;
    factory->setEquivalentBus(buses[0]);
    iVec = p_ivecMap->mapToVector();
    boost::shared_ptr<gridpack::math::Vector> voltage;
    voltage.reset(iVec->clone());

    // Create linear solver and solve for first injection voltages
    boost::shared_ptr<gridpack::math::LinearSolver> solver;
    solver.reset(new gridpack::math::LinearSolver (*Ybus));
    solver->configure(cursor);
    voltage->zero();
    solver->solve(*iVec,*voltage);
    addColumn(voltage.get(), buses[0]);

    // Evaluate remaining injection vectors
    for (i=1; i<buses.size(); i++) {
      factory->setEquivalentBus(buses[i]);
      iVec = p_ivecMap->mapToVector();
      voltage->zero();
      solver->solve(*iVec,*voltage);
      addColumn(voltage.get(), buses[i]);
    }
  }
}

/**
 * Basic destructor
 */
EquivalentMatrix::~EquivalentMatrix()
{
  if (p_initialized) NGA_Destroy(p_ga);
  p_map.clear();
}

/**
 * Add column to matrix
 * @param vector[in] add distributed vector to matrix
 * @param bus_col[in] ID of bus corresponding to column location
 */
void EquivalentMatrix::addColumn(gridpack::math::Vector *vector, int bus_col)
{
  if (!p_initialized) return;
  p_ivecMap->mapToBus(*vector);
  int Jdx = p_map.find(bus_col)->second;
  std::map<int,int>::iterator it = p_map.begin();
  int lo, hi;
  int one = 1;
  vector->localIndexRange(lo, hi);
  while (it != p_map.end()) {
    int busID = it->first;
    int Idx = it->second;
    std::vector<int> buses = p_network->getLocalBusIndices(busID);
    // if bus has some elements, then check to see if bus is locally owned
    if (buses.size() > 0) {
      int i;
      for (i=0; i<buses.size(); i++) {
        if (p_network->getActiveBus(buses[i])) {
          ComplexType value;
          value = p_network->getBus(buses[i])->getInjectionVoltage();
          int tlo[2], thi[2];
          tlo[0] = Idx;
          thi[0] = Idx;
          tlo[1] = Jdx;
          thi[1] = Jdx;
          NGA_Put(p_ga,tlo,thi,&value,&one);
        }
      }
    }
    it++;
  }
  NGA_Sync();
}

/**
 * Add row to matrix
 * @param vector[in] add distributed vector to matrix
 * @param bus_row[in] ID of bus corresponding to row location
 */
void EquivalentMatrix::addRow(gridpack::math::Vector *vector, int bus_row)
{
  if (!p_initialized) return;
  p_ivecMap->mapToBus(*vector);
  int Idx = p_map.find(bus_row)->second;
  std::map<int,int>::iterator it = p_map.begin();
  int lo, hi;
  int one = 1;
  vector->localIndexRange(lo, hi);
  while (it != p_map.end()) {
    int busID = it->first;
    int Jdx = it->second;
    std::vector<int> buses = p_network->getLocalBusIndices(busID);
    // if bus has some elements, then check to see if bus is locally owned
    if (buses.size() > 0) {
      int i;
      for (i=0; i<buses.size(); i++) {
        if (p_network->getActiveBus(buses[i])) {
          ComplexType value;
          value = p_network->getBus(buses[i])->getInjectionVoltage();
          int tlo[2], thi[2];
          tlo[0] = Idx;
          thi[0] = Idx;
          tlo[1] = Jdx;
          thi[1] = Jdx;
          NGA_Put(p_ga,tlo,thi,&value,&one);
        }
      }
    }
    it++;
  }
  NGA_Sync();
}

/**
 * Get a complete copy of matrix and copy it into a vector
 * @param[out] matrix vector data structure containing full matrix
 * @param[out] n dimension of matrix
 */
void EquivalentMatrix::getMatrix(std::vector<ComplexType> &matrix, int &n)
{
  if (!p_initialized) return;
  int lo[2], hi[2], ld;
  int nelem = p_dim*p_dim;
  n = p_dim;
  ld = p_dim;
  lo[0] = 0;
  lo[1] = 0;
  hi[0] = p_dim-1;
  hi[1] = p_dim-1;
  matrix.clear();
  matrix.resize(nelem);
  NGA_Get(p_ga,lo,hi,&matrix[0],&ld);
}

}  // dynamic_simulation
}  // gridpack
