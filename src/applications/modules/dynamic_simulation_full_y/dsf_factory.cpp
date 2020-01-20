/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dsf_factory.cpp
 * @author Shuangshuang Jin 
 * @date   Feb 04, 2015
 * @last modified date   May 13, 2015
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "dsf_factory.hpp"

namespace gridpack {
namespace dynamic_simulation {

// Dynamic simulation factory class implementations

/**
 * Basic constructor
 * @param network: network associated with factory
 */
DSFullFactory::DSFullFactory(DSFullFactory::NetworkPtr network)
  : gridpack::factory::BaseFactory<DSFullNetwork>(network)
{
  p_network = network;

  int i;
  p_numBus = p_network->numBuses();
  p_numBranch = p_network->numBranches();

  p_buses = new gridpack::dynamic_simulation::DSFullBus*[p_numBus];
  for (i=0; i<p_numBus; i++) {
    p_buses[i] = dynamic_cast<DSFullBus*>(p_network->getBus(i).get());
  }
  p_branches = new gridpack::dynamic_simulation::DSFullBranch*[p_numBranch];
  for (i=0; i<p_numBranch; i++) {
    p_branches[i] = dynamic_cast<DSFullBranch*>(p_network->getBranch(i).get());
  }
#ifdef USE_FNCS
  p_hash = new gridpack::hash_distr::HashDistribution<DSFullNetwork,
               gridpack::ComplexType,gridpack::ComplexType>
               (network);
  int dsize = sizeof(gridpack::dynamic_simulation::DSFullBus::voltage_data);
  p_busIO = new gridpack::serial_io::SerialBusIO<DSFullNetwork>(dsize,
                network);
#endif
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation::DSFullFactory::~DSFullFactory()
{
  delete [] p_buses;
  delete [] p_branches;
#ifdef USE_FNCS
  delete p_busIO;
  delete p_hash;
#endif
}

/**
 * Create the admittance (Y-Bus) matrix
 */
void gridpack::dynamic_simulation::DSFullFactory::setYBus(void)
{
  int i;

  // Invoke setYBus method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->setYBus();
  }

  // Invoke setYBus method on all branch objects
  for (i=0; i<p_numBranch; i++) {
    p_branches[i]->setYBus();
  }
}

/**
 * Create the admittance (Y-Bus) matrix for Branch only
 */
void gridpack::dynamic_simulation::DSFullFactory::setYBranch(void)
{
  int i;

  // Invoke setYBus method on all branch objects
  for (i=0; i<p_numBranch; i++) {
    p_branches[i]->setYBus();
  }
}

/**
 * Get the updating factor for posfy11 stage ybus
 */
gridpack::ComplexType
gridpack::dynamic_simulation::DSFullFactory::setFactor(int sw2_2, int sw3_2)
{
  gridpack::ComplexType dummy(-999.0, -999.0);

  int i;

  // Invoke getPosfy11YbusUpdateFactor method on all branch objects
  for (i=0; i<p_numBranch; i++) {
    gridpack::ComplexType ret = p_branches[i]
      ->getPosfy11YbusUpdateFactor(sw2_2, sw3_2);
    if (ret != dummy) {
      return ret;
    }
  }
}

/**
 * Apply an event to all branches in the system
 * @param event a struct describing a fault
 */
void gridpack::dynamic_simulation::DSFullFactory::setEvent(const
    gridpack::dynamic_simulation::Event &event)
{
  int i;
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->clearEvent();
  }
  for (i=0; i<p_numBranch; i++) {
    p_branches[i]->setEvent(event);
  }
}

/**
 * Check network to see if there is a process with no generators
 * @return true if all processors have at least on generator
 */
bool gridpack::dynamic_simulation::DSFullFactory::checkGen(void)
{
  int i, count;
  count = 0;
  for (i=0; i<p_numBus; i++) {
    if (p_network->getActiveBus(i)) {
      count += p_buses[i]->getNumGen();
    }
  }
  //printf("p[%d] number of generators: %d\n",p_network->communicator().rank(),count);
  int iok = 0;
  if (count > 0) iok = 1;
  int ok;
  boost::mpi::all_reduce(p_network->communicator(),iok,ok,std::plus<int>());
  bool ret = false;
  if (ok == p_network->communicator().size()) ret = true;
  return ret;
}

/**
 * Initialize init vectors for integration
 * @param ts time step
 */
void gridpack::dynamic_simulation::DSFullFactory::initDSVect(double ts)
{
  int i;

  // Invoke initDSVect method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->initDSVect(ts);
  }
}

/**
 * Update vectors in each integration time step (Predictor)
 */
void gridpack::dynamic_simulation::DSFullFactory::predictor_currentInjection(bool flag)
{
  int i;

  // Invoke method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->predictor_currentInjection(flag);
  }
}

/**
 * Update vectors in each integration time step (Predictor)
 */
void gridpack::dynamic_simulation::DSFullFactory::predictor(double t_inc, bool flag)
{
  int i;

  // Invoke updateDSVect method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->predictor(t_inc,flag);
  }
}

/**
 * Update vectors in each integration time step (Corrector)
 */
void gridpack::dynamic_simulation::DSFullFactory::corrector_currentInjection(bool flag)
{
  int i;

  // Invoke method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->corrector_currentInjection(flag);
  }
}

/**
 * Update vectors in each integration time step (Corrector)
 */
void gridpack::dynamic_simulation::DSFullFactory::corrector(double t_inc, bool flag)
{
  int i;

  // Invoke updateDSVect method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->corrector(t_inc,flag);
  }
}

/**
 * Update dynamic load internal relays action
 */
void gridpack::dynamic_simulation::DSFullFactory::dynamicload_post_process(double t_inc, bool flag)
{
	int i;

    for (i=0; i<p_numBus; i++) {
      p_buses[i]->dynamicload_post_process(t_inc,flag);
    }
}

/**
 * Set volt from volt_full
 */
void gridpack::dynamic_simulation::DSFullFactory::setVolt(bool flag)
{
  int i;

  // Invoke updateDSVect method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->setVolt(flag);
  }
}

/**
 * update bus frequecy
 */
void gridpack::dynamic_simulation::DSFullFactory::updateBusFreq(double delta_t)
{
  int i;

  // Invoke updateDSVect method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->updateFreq(delta_t);
  }
}

void gridpack::dynamic_simulation::DSFullFactory::setWideAreaFreqforPSS(double freq)
{
  int i;

  // Invoke updateDSVect method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->setWideAreaFreqforPSS(freq);
  }
}

std::vector<double> gridpack::dynamic_simulation::DSFullFactory::grabWideAreaFreq( )  //renke hardcoded
{
  int i;

  // Invoke updateDSVect method on all bus objects
  double freq1, freq2, freqdiff;
  freq1 = 60.0;
  freq2 = 60.0;
  freqdiff = 0.0;
  int busnumtmp;
  
  // the following code grabs specific bus frequency
  /*
  for (i=0; i<p_numBus; i++) {
	  
    busnumtmp = p_buses[i]->getOriginalIndex();
	if (busnumtmp == 34){
		freq1 = p_buses[i]->getBusVolFrequency();
	}
	if (busnumtmp == 30){
		freq2 = p_buses[i]->getBusVolFrequency();
	}
	
  }
  */
  
  freqdiff = freq1 - freq2;
  
  //printf("-----------!!renke debug DSFullFactory::grabWideAreaFreq( ): bus 30: %12.6f,  bus 34: %12.6f,  diff 34-30: %12.6f \n", freq2, freq1, freqdiff);
  
  std::vector <double> vbusfreq;
  vbusfreq.push_back(freq2);
  vbusfreq.push_back(freq1);
  vbusfreq.push_back(freqdiff);
  
  //return freqdiff;
  return vbusfreq;
}

bool gridpack::dynamic_simulation::DSFullFactory::updateBusRelay(bool flag,double delta_t)
{
	int i;
	bool bflag;
	
	bflag = false;
	
	for (i=0; i<p_numBus; i++) {
		bflag = bflag || p_buses[i]->updateRelay(flag,delta_t);
	
        }
        return checkTrueSomewhere(bflag); 
	//return bflag;
	
}

bool gridpack::dynamic_simulation::DSFullFactory::updateBranchRelay(bool flag,double delta_t)
{
	int i;
	bool bflag;
	
	bflag = false;
	
	for (i=0; i<p_numBranch; i++) {
		bflag = bflag || p_branches[i]->updateRelay(flag,delta_t);
	}
	
	return bflag;
	
}

/**
 * update old bus voltage, renke add
 */
void gridpack::dynamic_simulation::DSFullFactory::updateoldbusvoltage()
{
	int i;
	
	for (i=0; i<p_numBus; i++) {
    p_buses[i]->updateoldbusvoltage();
	
    }
}

/**
 * print all bus voltage, renke add
 */
void gridpack::dynamic_simulation::DSFullFactory::printallbusvoltage()
{
	int i;
	
	for (i=0; i<p_numBus; i++) {
    p_buses[i]->printbusvoltage();
	
    }
}

/**
 * Add constant impedance load admittance to diagonal elements of
 * Y-matrix
 */
void gridpack::dynamic_simulation::DSFullFactory::addLoadAdmittance()
{
  int i;

  // Invoke addLoadAdmittance method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->addLoadAdmittance();
  }
}

bool gridpack::dynamic_simulation::DSFullFactory::securityCheck()
{
  int i;
  bool secure;

  double maxAngle = -999999.0;
  double minAngle = 999999.0;
  double angle;
  for (i = 0; i < p_numBus; i++) {
    int nGen = p_buses[i]->getNumGen();
    if (nGen > 0) {
      angle = p_buses[i]->getAngle();
      //printf("angle = %f\n", angle);
      if (angle > maxAngle) {
        maxAngle = angle;
        //printf("   max = %f\n", maxAngle);
      }
      if (angle < minAngle) {
        minAngle = angle;
        //printf("   min = %f\n", minAngle);
      }
      double pi = 4.0*atan(1.0);
      if ((maxAngle - minAngle) > 360 * pi / 180)
        secure = false;
      else
        secure = true; 
    }
  }
  //printf("max = %f, min = %f\n", maxAngle, minAngle);
  return secure;
}

/**
 * Scale generator real power. If zone less than 1 then scale all
 * generators in the area
 * @param scale factor to scale real power generation
 * @param area index of area for scaling generation
 * @param zone index of zone for scaling generation
 */
void gridpack::dynamic_simulation::DSFullFactory::scaleGeneratorRealPower(
    double scale, int area, int zone)
{
  int nbus = p_network->numBuses();
  int i, izone;
  for (i=0; i<nbus; i++) {
    gridpack::dynamic_simulation::DSFullBus *bus = p_network->getBus(i).get();
    boost::shared_ptr<gridpack::component::DataCollection> data
      = p_network->getBusData(i);
    if (zone > 0) {
      izone = bus->getZone();
    } else {
      izone = zone;
    }
    if (bus->getArea() == area && zone == izone) {
      std::vector<std::string> tags = bus->getGenerators();
      int j;
      for (j=0; j<tags.size(); j++) {
        bus->scaleGeneratorRealPower(tags[j], scale, data);
      }
    }
  }
}

/**
 * Scale load power. If zone less than 1 then scale all
 * loads in the area
 * @param scale factor to scale load real power
 * @param area index of area for scaling load
 * @param zone index of zone for scaling load
 * @return false if there is not enough capacity to change generation
 *         by requested amount
 */
void gridpack::dynamic_simulation::DSFullFactory::scaleLoadPower(
    double scale, int area, int zone)
{
  int nbus = p_network->numBuses();
  int i, izone;
  for (i=0; i<nbus; i++) {
    gridpack::dynamic_simulation::DSFullBus *bus = p_network->getBus(i).get();
    if (zone > 0) { 
      izone = bus->getZone();
    } else { 
      izone = zone;
    } 
    if (bus->getArea() == area && zone == izone) {
      std::vector<std::string> tags = bus->getLoads();
      int j;
      for (j=0; j<tags.size(); j++) {
        bus->scaleLoadPower(tags[j],scale);
      }
    }
  }
}

/**
 * Return the total real power load for all loads in the zone. If zone
 * less than 1, then return the total load for the area
 * @param area index of area
 * @param zone index of zone
 * @return total load
 */
double gridpack::dynamic_simulation::DSFullFactory::getTotalLoadRealPower(
    int area, int zone)
{
  double ret = 0.0;
  int nbus = p_network->numBuses();
  int i, j, izone;
  for (i=0; i<nbus; i++) {
    gridpack::dynamic_simulation::DSFullBus *bus = p_network->getBus(i).get();
    if (zone > 0) {
      izone = bus->getZone();
    } else {
      izone = zone;
    }
    if (bus->getArea() == area && zone == izone) {
      std::vector<std::string> tags;
      std::vector<double> pl;
      std::vector<double> ql;
      std::vector<int> status;
      bus->getLoadPower(tags,pl,ql,status);
      for (j=0; j<pl.size(); j++) {
        if (static_cast<bool>(status[j])) {
          ret += pl[j];
        }
      }
    }
  }
  p_network->communicator().sum(&ret,1);
  return ret;
}

/**
 * Return the current real power generation and the maximum and minimum total
 * power generation for all generators in the zone. If zone is less than 1
 * then return values for all generators in the area
 * @param area index of area
 * @param zone index of zone
 * @param total total real power generation
 * @param pmin minimum allowable real power generation
 * @param pmax maximum available real power generation
 */
void gridpack::dynamic_simulation::DSFullFactory::getGeneratorMargins(int area,
    int zone, double *total, double *pmin, double *pmax)
{
  *total = 0.0;
  *pmin = 0.0;
  *pmax = 0.0;
  int nbus = p_network->numBuses();
  int i, j, izone;
  for (i=0; i<nbus; i++) {
    gridpack::dynamic_simulation::DSFullBus *bus = p_network->getBus(i).get();
    if (zone > 0) {
      izone = bus->getZone();
    } else {
      izone = zone;
    }
    if (bus->getArea() == area && zone == izone) {
      std::vector<std::string> tags;
      std::vector<double> tcurrent;
      std::vector<double> tpmin;
      std::vector<double> tpmax;
      std::vector<int> status;
      bus->getGeneratorMargins(tags,tcurrent,tpmin,tpmax,status);
      for (j=0; j<tcurrent.size(); j++) {
        if (static_cast<bool>(status[j])) {
          *total += tcurrent[j];
          *pmin += tpmin[j];
          *pmax += tpmax[j];
        }
      }
    }
  }
}

/**
 * Reset real power of loads and generators to original values
 */
void gridpack::dynamic_simulation::DSFullFactory::resetPower()
{
  int nbus = p_network->numBuses();
  int i;
  for (i=0; i<nbus; i++) {
    boost::shared_ptr<gridpack::component::DataCollection> data
      = p_network->getBusData(i);
    p_network->getBus(i)->resetPower(data);
  }
}


/**
 * load parameters for the extended buses from composite load model
 */
void gridpack::dynamic_simulation::DSFullFactory::setExtendedCmplBusVoltage()
{
	int i;

  // Invoke updateDSVect method on all bus objects
  for (i=0; i<p_numBus; i++) {
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(p_buses[i])->setExtendedCmplBusVoltage(p_network->getBusData(i));
  }
	
}

/**
 * load parameters for the extended buses from composite load model
 */
void gridpack::dynamic_simulation::DSFullFactory::LoadExtendedCmplBus()
{
	int i;

  // Invoke updateDSVect method on all bus objects
  for (i=0; i<p_numBus; i++) {
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(p_buses[i])->LoadExtendedCmplBus(p_network->getBusData(i));
  }
	
}

#ifdef USE_FNCS
/**
 * Scatter load from FNCS framework to buses
 */
void gridpack::dynamic_simulation::DSFullFactory::scatterLoad()
{
  std::vector<int> keys;
  std::vector<gridpack::ComplexType> values;

  if (p_network->communicator().rank() == 0) {
    /**
     * Add FNCS code here that puts bus IDs in keys and corresponding LOAD_PL and
     * LOAD_QL as the real and imaginary parts of a single complex number in
     * values
     */
  }

  // distribute loads to all processors
  p_hash->distributeBusValues(keys, values);

  // assign loads to individual processors
  int nvals = keys.size();
  int i;
  for (i=0; i<nvals; i++) {
    p_buses[keys[i]]->setLoad(real(values[i]), imag(values[i]));
  }
}

/**
 * Gather voltage magnitude and phase angle to root node and transmit to
 * FNCS framework
 */
void gridpack::dynamic_simulation::DSFullFactory::gatherVoltage()
{
  // Get data from network buses
  std::vector<gridpack::dynamic_simulation::DSFullBus::voltage_data> data;
  p_busIO->gatherData(data);
  
  // Loop over all data items
  int i;
  if (p_network->communicator().rank() == 0) {
    int nvals = data.size();
    for (i=0; i<nvals; i++) {
      int id = data[i].busID;
      gridpack::ComplexType voltage = data[i].voltage;
      /**
       * Add FNCS code here that maps complex voltage to bus id
       */
    }
  }
}

#endif

/**
 * Set parameters for real time path rating diagnostics
 * @param src_area generation area
 * @param src_zone generation zone
 * @param load_area load area
 * @param load_zone load zone
 * @param gen_scale scale factor for generation
 * @param load_scale scale factor for loads
 */
void gridpack::dynamic_simulation::DSFullFactory::setRTPRParams(
    int src_area, int src_zone, int load_area,
    int load_zone, double gen_scale, double load_scale)
{
  int nbus = p_network->numBuses();
  int i, j, izone;
  for (i=0; i<nbus; i++) {
    gridpack::dynamic_simulation::DSFullBus *bus = p_network->getBus(i).get();
    int tarea = bus->getArea();
    int tzone = bus->getZone();
    printf("RTPR Parameters bus: %d Area: %d Zone: %d\n",
        bus->getOriginalIndex(),tarea,tzone);
    if (src_zone > 0) {
      izone = tzone;
    } else {
      izone = src_zone;
    }
    bus->setScale(1.0);
    if (tarea == src_area && src_zone == izone) {
      bus->setSource(true);
      bus->setScale(gen_scale);
    } else {
      bus->setSource(false);
    }
    if (load_zone > 0) {
      izone = tzone;
    } else {
      izone = load_zone;
    }
    if (tarea == load_area && load_zone == izone) {
    printf("RTPR Parameters set sink bus: %d Area: %d Zone: %d\n",
        bus->getOriginalIndex(),tarea,tzone);
      bus->setSink(true);
      bus->setScale(load_scale);
    } else {
      bus->setSink(false);
    }
  }
}

} // namespace dynamic_simulation
} // namespace gridpack
