/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   rtpr_driver.cpp
 * @author Bruce Palmer
 * @date   2019-10-09 13:12:46 d3g293
 *
 * @brief Driver for real-time path rating calculation that make use of the
 *        powerflow module to implement individual power flow simulations for
 *        each contingency. The different contingencies are distributed across
 *        separate communicators using the task manager.
 *
 *
 */
// -------------------------------------------------------------

#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"
#include "rtpr_driver.hpp"

#define RTPR_DEBUG
#define USE_STATBLOCK

//#define USE_SUCCESS
// Sets up multiple communicators so that individual contingency calculations
// can be run concurrently

/**
 * Basic constructor
 */
gridpack::rtpr::RTPRDriver::RTPRDriver(void)
{
}

/**
 * Basic destructor
 */
gridpack::rtpr::RTPRDriver::~RTPRDriver(void)
{
}

/**
 * Get list of contingencies from external file
 * @param cursor pointer to contingencies in input deck
 * @return vector of contingencies
 */
std::vector<gridpack::powerflow::Contingency>
  gridpack::rtpr::RTPRDriver::getContingencies(
      gridpack::utility::Configuration::ChildCursors &contingencies)
{
  // The contingencies ChildCursors argument is a vector of configuration
  // pointers. Each element in the vector is pointing at a seperate Contingency
  // block within the Contingencies block in the input file.
  std::vector<gridpack::powerflow::Contingency> ret;
  int size = contingencies.size();
  int i, idx;
  // Create string utilities object to help parse file
  gridpack::utility::StringUtils utils;
  // Loop over all child cursors
  printf("p[%d] Reading contingencies size: %d\n",p_world.rank(),size);
  for (idx = 0; idx < size; idx++) {
    std::string ca_type;
    contingencies[idx]->get("contingencyType",&ca_type);
    // Contingency name is used to direct output to different files for each
    // contingency
    std::string ca_name;
    contingencies[idx]->get("contingencyName",&ca_name);
    if (ca_type == "Line") {
      printf("p[%d] Found line contingency\n",p_world.rank());
      std::string buses;
      contingencies[idx]->get("contingencyLineBuses",&buses);
      std::string names;
      contingencies[idx]->get("contingencyLineNames",&names);
      // Tokenize bus string to get a list of individual buses
      std::vector<std::string> string_vec = utils.blankTokenizer(buses);
      // Convert buses from character strings to ints
      std::vector<int> bus_ids;
      for (i=0; i<string_vec.size(); i++) {
        bus_ids.push_back(atoi(string_vec[i].c_str()));
      }
      string_vec.clear();
      // Tokenize names string to get a list of individual line tags
      string_vec = utils.blankTokenizer(names);
      std::vector<std::string> line_names;
      // clean up line tags so that they are exactly two characters
      for (i=0; i<string_vec.size(); i++) {
        line_names.push_back(utils.clean2Char(string_vec[i]));
      }
      // Check to make sure we found everything
      if (bus_ids.size() == 2*line_names.size()) {
        // Add contingency parameters to contingency struct
        gridpack::powerflow::Contingency contingency;
        contingency.p_name = ca_name;
        contingency.p_type = Branch;
        int i;
        for (i = 0; i < line_names.size(); i++) {
          contingency.p_from.push_back(bus_ids[2*i]);
          contingency.p_to.push_back(bus_ids[2*i+1]);
          contingency.p_ckt.push_back(line_names[i]);
          contingency.p_saveLineStatus.push_back(true);
        }
        // Add branch contingency to contingency list
        ret.push_back(contingency);
      }
    } else if (ca_type == "Generator") {
      printf("p[%d] Found generator contingency\n",p_world.rank());
      std::string buses;
      contingencies[idx]->get("contingencyBuses",&buses);
      std::string gens;
      contingencies[idx]->get("contingencyGenerators",&gens);
      // Tokenize bus string to get a list of individual buses
      std::vector<std::string> string_vec = utils.blankTokenizer(buses);
      std::vector<int> bus_ids;
      // Convert buses from character strings to ints
      for (i=0; i<string_vec.size(); i++) {
        bus_ids.push_back(atoi(string_vec[i].c_str()));
      }
      string_vec.clear();
      // Tokenize gens string to get a list of individual generator tags
      string_vec = utils.blankTokenizer(gens);
      std::vector<std::string> gen_ids;
      // clean up generator tags so that they are exactly two characters
      for (i=0; i<string_vec.size(); i++) {
        gen_ids.push_back(utils.clean2Char(string_vec[i]));
      }
      // Check to make sure we found everything
      if (bus_ids.size() == gen_ids.size()) {
        gridpack::powerflow::Contingency contingency;
        contingency.p_name = ca_name;
        contingency.p_type = Generator;
        int i;
        for (i = 0; i < bus_ids.size(); i++) {
          contingency.p_busid.push_back(bus_ids[i]);
          contingency.p_genid.push_back(gen_ids[i]);
          contingency.p_saveGenStatus.push_back(true);
        }
        // Add generator contingency to contingency list
        ret.push_back(contingency);
      }
    }
  }
  return ret;
}

/**
 * Dummy struct to get around some problems with vectors of characters
 */
struct char2 {
  char str[3];
};

/**
 * Create a list of all N-1 generator contingencies for a given area and
 * zone. If zone less than 1, then just use the entire area.
 * @param network power grid network on which contingencies are defined 
 * @param area index of area that will generate contingencies
 * @param zone index of zone that will generate contingencies
 * @return vector of contingencies
 */
std::vector<gridpack::powerflow::Contingency> 
  gridpack::rtpr::RTPRDriver::createGeneratorContingencies(
  boost::shared_ptr<gridpack::powerflow::PFNetwork> network, int area,
  int zone)
{
  std::vector<gridpack::powerflow::Contingency> ret;
  gridpack::utility::StringUtils util;
  int nbus = network->numBuses();
  int i,j,iarea, izone;
  std::vector<int> bus_ids;
  std::vector<char2> tags;
  char2 buf;
  int nproc = network->communicator().size();
  int me = network->communicator().rank();
  // Loop over all buses
  for (i=0; i<nbus; i++) {
    if (network->getActiveBus(i)) {
      // Get data collection object
      gridpack::component::DataCollection *data
        = network->getBusData(i).get();
      int ngen=0; 
      data->getValue(BUS_AREA, &iarea);
      if (zone > 0) {
        data->getValue(BUS_ZONE, &izone);
      } else {
        izone = zone;
      }
      if (iarea == area && izone == zone) {
        if (data->getValue(GENERATOR_NUMBER, &ngen)) {
          for (j=0; j<ngen; j++) {
            // Check generator status
            int status = 0;
            data->getValue(GENERATOR_STAT, &status, j);
            if (status != 0) {
              // Generator is active in base case. Add a contingency
              bus_ids.push_back(network->getOriginalBusIndex(i));
              std::string tag, clean_tag;
              data->getValue(GENERATOR_ID,&tag,j);
              clean_tag = util.clean2Char(tag);
              strncpy(buf.str,clean_tag.c_str(),2);
              buf.str[2] = '\0';
              tags.push_back(buf);
            }
          }
        }
      }
    }
  }
  // Have a list of local faults. Find out total number of faults.
  int nflt = bus_ids.size();
  std::vector<int> sizes(nproc);
  for (i=0; i<nproc; i++) sizes[i] = 0;
  sizes[me] = nflt;
  network->communicator().sum(&sizes[0], nproc);
  // create global list containing buses and generator tags
  nflt = 0;
  for (i=0; i<nproc; i++) nflt += sizes[i];
  if (nflt == 0) {
    return ret;
  }
  // create a list of indices for local data values
  int offset = 0;
  for (i=0; i<me; i++) offset += sizes[i];
  std::vector<int> idx;
  for (i=0; i < bus_ids.size(); i++) idx.push_back(offset+i);
  // upload all values to a global vector
  gridpack::parallel::GlobalVector<int> busIDs(network->communicator());
  gridpack::parallel::GlobalVector<char2> genIDs(network->communicator());
  busIDs.addElements(idx,bus_ids);
  genIDs.addElements(idx,tags);
  busIDs.upload();
  genIDs.upload();
  busIDs.getAllData(bus_ids);
  genIDs.getAllData(tags);
  // create a list of contingencies
  char sbuf[128];
  for (i=0; i<nflt; i++) {
    gridpack::powerflow::Contingency fault;
    fault.p_type = Generator;
    fault.p_busid.push_back(bus_ids[i]);
    std::string tag = tags[i].str;
    fault.p_genid.push_back(util.clean2Char(tag));
    fault.p_saveGenStatus.push_back(true);
    sprintf(sbuf,"GenCTG%d",i+1);
    fault.p_name = sbuf;
    ret.push_back(fault);
  }
  return ret;
}

/**
 * Create a list of all N-1 branch contingencies for a given area and
 * zone. If zone less than 1, then just use the entire area.
 * @param network power grid network on which contingencies are defined 
 * @param area index of area that will generate contingencies
 * @param zone index of zone that will generate contingencies
 * @return vector of contingencies
 */
std::vector<gridpack::powerflow::Contingency>
  gridpack::rtpr::RTPRDriver::createBranchContingencies(
  boost::shared_ptr<gridpack::powerflow::PFNetwork> network, int area,
  int zone)
{
  std::vector<gridpack::powerflow::Contingency> ret;
  gridpack::utility::StringUtils util;
  int nbranch = network->numBranches();
  int i,j,idx1,idx2,i1,i2,area1,area2,zone1,zone2;
  std::vector<int> bus_from;
  std::vector<int> bus_to;
  std::vector<char2> tags;
  char2 buf;
  // Loop over all branches
  for (i=0; i<nbranch; i++) {
    if (network->getActiveBranch(i)) {
      // Get data collection object
      gridpack::component::DataCollection *data
        = network->getBranchData(i).get();
      int nline=0; 
      network->getBranchEndpoints(i,&i1,&i2);
      network->getBusData(i1)->getValue(BUS_AREA,&area1);
      network->getBusData(i2)->getValue(BUS_AREA,&area2);
      if (zone > 0) {
        network->getBusData(i1)->getValue(BUS_ZONE,&zone1);
        network->getBusData(i2)->getValue(BUS_ZONE,&zone2);
      } else {
        zone1 = zone;
        zone2 = zone;
      }
      if (area == area1 && area == area2 && zone == zone1 && zone == zone2) {
        if (data->getValue(BRANCH_NUM_ELEMENTS, &nline)) {
          for (j=0; j<nline; j++) {
            // Check generator status
            int status = 0;
            data->getValue(BRANCH_STATUS, &status, j);
            if (status != 0) {
              // Generator is active in base case. Add a contingency
              network->getOriginalBranchEndpoints(i,&idx1,&idx2);
              bus_from.push_back(idx1);
              bus_to.push_back(idx2);
              std::string tag, clean_tag;
              data->getValue(BRANCH_CKT,&tag,j);
              clean_tag = util.clean2Char(tag);
              strncpy(buf.str,clean_tag.c_str(),2);
              buf.str[2] = '\0';
              tags.push_back(buf);
            }
          }
        }
      }
    }
  }
  // Have a list of local faults. Find out total number of faults.
  int nflt = bus_from.size();
  int nproc = network->communicator().size();
  int me = network->communicator().rank();
  std::vector<int> sizes(nproc);
  for (i=0; i<nproc; i++) sizes[i] = 0;
  sizes[me] = nflt;
  network->communicator().sum(&sizes[0], nproc);
  // create global list containing buses and generator tags
  nflt = 0;
  for (i=0; i<nproc; i++) nflt += sizes[i];
  if (nflt == 0) {
    return ret;
  }
  // create a list of indices for local data values
  int offset = 0;
  for (i=0; i<me; i++) offset += sizes[i];
  std::vector<int> idx;
  for (i=0; i < bus_from.size(); i++) idx.push_back(offset+i);
  // upload all values to a global vector
  gridpack::parallel::GlobalVector<int> fromIDs(network->communicator());
  gridpack::parallel::GlobalVector<int> toIDs(network->communicator());
  gridpack::parallel::GlobalVector<char2> branchIDs(network->communicator());
  fromIDs.addElements(idx,bus_from);
  toIDs.addElements(idx,bus_to);
  branchIDs.addElements(idx,tags);
  fromIDs.upload();
  toIDs.upload();
  branchIDs.upload();
  fromIDs.getAllData(bus_from);
  toIDs.getAllData(bus_to);
  branchIDs.getAllData(tags);
  // create a list of contingencies
  char sbuf[128];
  for (i=0; i<nflt; i++) {
    gridpack::powerflow::Contingency fault;
    fault.p_type = Branch;
    fault.p_from.push_back(bus_from[i]);
    fault.p_to.push_back(bus_to[i]);
    std::string tag = tags[i].str;
    fault.p_ckt.push_back(util.clean2Char(tag));
    fault.p_saveLineStatus.push_back(true);
    sprintf(sbuf,"LineCTG%d",i+1);
    fault.p_name = sbuf;
    ret.push_back(fault);
  }
  return ret;
}

/**
 * Get list of tie lines
 * @param cursor pointer to tie lines in input deck
 * @return vector of contingencies
 */
std::vector<gridpack::rtpr::TieLine>
  gridpack::rtpr::RTPRDriver::getTieLines(
      gridpack::utility::Configuration::ChildCursors &tielines)
{
  // The ChildCursors argument is a vector of configuration
  // pointers. Each element in the vector is pointing at a seperate tieLine
  // block within the TieLines block in the input file.
  std::vector<gridpack::rtpr::TieLine> ret;
  int size = tielines.size();
  int i, idx;
  // Create string utilities object to help parse file
  gridpack::utility::StringUtils utils;
  // Loop over all child cursors
  for (idx = 0; idx < size; idx++) {
    std::string branchIDs;
    tielines[idx]->get("Branch",&branchIDs);
    std::vector<std::string> from_to = utils.blankTokenizer(branchIDs);
    gridpack::rtpr::TieLine tieline;
    tieline.from = atoi(from_to[0].c_str());
    tieline.to = atoi(from_to[1].c_str());
    std::string tag;
    tielines[idx]->get("Tag",&tag);
    std::string clean_tag = utils.clean2Char(tag);
    strncpy(tieline.tag,clean_tag.c_str(),2);
    tieline.tag[2] = '\0';
    ret.push_back(tieline);
  }
  return ret;
}

/**
 * Automatically create list of tie lines between area1,zone1 and
 * area2,zone2
 * @param area1 index of source area
 * @param zone1 index of source zone
 * @param area2 index of destination area
 * @param zone2 index of destination zone
 * @param tielines list of tie lines between area1,zone1 and area2,zone2
 */
std::vector<gridpack::rtpr::TieLine> 
  gridpack::rtpr::RTPRDriver::getTieLines(int area1, int zone1,
    int area2, int zone2)
{
  std::vector<gridpack::rtpr::TieLine> ret;
  // Loop through all branches and find the ones that connect area1,zone1 with
  // area2,zone2
  int nbranch = p_pf_network->numBranches();
  int i, j, idx1, idx2;
  int iarea1, iarea2, izone1, izone2, tzone1, tzone2;
  bool found;
  for (i=0; i<nbranch; i++) {
    if (p_pf_network->getActiveBranch(i)) {
      gridpack::powerflow::PFBranch *branch = p_pf_network->getBranch(i).get();
      p_pf_network->getBranchEndpoints(i,&idx1,&idx2);
      gridpack::powerflow::PFBus *bus1 = p_pf_network->getBus(idx1).get();
      gridpack::powerflow::PFBus *bus2 = p_pf_network->getBus(idx2).get();
      iarea1 = bus1->getArea();
      iarea2 = bus2->getArea();
      izone1 = bus1->getZone();
      izone2 = bus2->getZone();
      // Check bus1 corresponds to area1,zone1, bus 2 corresponds to area2/zone2
      found = false;
      if (zone1 > 0) {
        tzone1 = izone1;
      } else {
        tzone1 = zone1;
      }
      if (zone2 > 0) {
        tzone2 = izone2;
      } else {
        tzone2 = zone2;
      }
      if (iarea1 == area1 && tzone1 == zone1 &&
          iarea2 == area2 && tzone2 == zone2) {
        found = true;
        int from, to;
        p_pf_network->getOriginalBranchEndpoints(i,&from,&to);
        std::vector<std::string> tags = branch->getLineIDs();
        for (j=0; j<tags.size(); j++) {
          gridpack::rtpr::TieLine tieline;
          tieline.from = from;
          tieline.to = to;
          strncpy(tieline.tag,tags[j].c_str(),2);
          tieline.tag[2] = '\0';
          ret.push_back(tieline);
        }
      }
      // Check bus2 corresponds to area1,zone1, bus 1 corresponds to area2/zone2
      if (!found) {
        if (zone1 > 0) {
          tzone2 = izone1;
        } else {
          tzone2 = zone1;
        }
        if (zone2 > 0) {
          tzone1 = izone2;
        } else {
          tzone1 = zone2;
        }
        if (iarea2 == area1 && tzone1 == zone1 &&
            iarea1 == area2 && tzone2 == zone2) {
          int from, to;
          p_pf_network->getOriginalBranchEndpoints(i,&from,&to);
          std::vector<std::string> tags = branch->getLineIDs();
          for (j=0; j<tags.size(); j++) {
            gridpack::rtpr::TieLine tieline;
            tieline.from = from;
            tieline.to = to;
            strncpy(tieline.tag,tags[j].c_str(),2);
            tieline.tag[2] = '\0';
            ret.push_back(tieline);
          }
        }
      }
    }
  }
  // Have a list of local tie lines. Find out total number of tie lines.
  int ntie = ret.size();
  int nproc = p_pf_network->communicator().size();
  int me = p_pf_network->communicator().rank();
  std::vector<int> sizes(nproc);
  for (i=0; i<nproc; i++) sizes[i] = 0;
  sizes[me] = ntie;
  p_pf_network->communicator().sum(&sizes[0], nproc);
  // create global list of tie lines
  ntie = 0;
  for (i=0; i<nproc; i++) ntie += sizes[i];
  if (ntie == 0) {
    return ret;
  }
  // create a list of indices for local data values
  int offset = 0;
  for (i=0; i<me; i++) offset += sizes[i];
  std::vector<int> idx;
  for (i=0; i < ret.size(); i++) idx.push_back(offset+i);
  // upload all values to a global vector
  gridpack::parallel::GlobalVector<gridpack::rtpr::TieLine>
    tieLines(p_pf_network->communicator());
  tieLines.addElements(idx,ret);
  tieLines.upload();
  tieLines.getAllData(ret);
  return ret;
}

/**
 * Create a list of generators for a given area and zone. These
 * generators will be monitored to see if they exceed operating
 * specifications.
 * @param network power grid used for DS simulation
 * @param area index of area that will generate contingencies
 * @param zone index of zone that will generate contingencies
 * @param buses list of bus IDs that contain generators
 * @param tags list of generator IDs
 */
void  gridpack::rtpr::RTPRDriver::findWatchedGenerators(
  boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork> network, int area,
  int zone, std::vector<int> &buses, std::vector<std::string> &tags)
{
  gridpack::utility::StringUtils util;
  int nbus = network->numBuses();
  int i,j,iarea, izone;
  std::vector<int> bus_ids;
  std::vector<char2> ctags;
  char2 buf;
  int nproc = network->communicator().size();
  int me = network->communicator().rank();
  // Loop over all buses
  for (i=0; i<nbus; i++) {
    if (network->getActiveBus(i)) {
      // Get data collection object
      gridpack::component::DataCollection *data
        = network->getBusData(i).get();
      int ngen=0; 
      data->getValue(BUS_AREA, &iarea);
      if (zone > 0) {
        data->getValue(BUS_ZONE, &izone);
      } else {
        izone = zone;
      }
      if (iarea == area && izone == zone) {
        if (data->getValue(GENERATOR_NUMBER, &ngen)) {
          for (j=0; j<ngen; j++) {
            // Check generator status
            int status = 0;
            data->getValue(GENERATOR_STAT, &status, j);
            if (status != 0) {
              // Generator is active in base case. Add a contingency
              bus_ids.push_back(network->getOriginalBusIndex(i));
              std::string tag, clean_tag;
              data->getValue(GENERATOR_ID,&tag,j);
              clean_tag = util.clean2Char(tag);
              strncpy(buf.str,clean_tag.c_str(),2);
              buf.str[2] = '\0';
              ctags.push_back(buf);
            }
          }
        }
      }
    }
  }

  // Have a list of local faults. Find out total number of faults.
  int ngen = bus_ids.size();
  std::vector<int> sizes(nproc);
  for (i=0; i<nproc; i++) sizes[i] = 0;
  sizes[me] = ngen;
  network->communicator().sum(&sizes[0], nproc);
  ngen = 0;
  for (i=0; i<nproc; i++) ngen += sizes[i];
  if (ngen == 0) {
    return;
  }
  // create a list of indices for local data values
  int offset = 0;
  for (i=0; i<me; i++) offset += sizes[i];
  std::vector<int> idx;
  for (i=0; i < bus_ids.size(); i++) idx.push_back(offset+i);
  // upload all values to a global vector
  gridpack::parallel::GlobalVector<int> busIDs(network->communicator());
  gridpack::parallel::GlobalVector<char2> genIDs(network->communicator());
  busIDs.addElements(idx,bus_ids);
  genIDs.addElements(idx,ctags);
  busIDs.upload();
  genIDs.upload();
  busIDs.getAllData(bus_ids);
  genIDs.getAllData(ctags);

  ngen = bus_ids.size();
  buses.clear();
  tags.clear();
  for (i=0; i<ngen; i++) {
    buses.push_back(bus_ids[i]);
    tags.push_back(ctags[i].str);
  }
}

/**
 * Scale loads by rating parameter and adjust generation to match the change
 * in load
 * @param scale scale factor on loads
 * @param flag signal system that should be scaled
 *             0: powerflow
 *             1: dynamic simulation
 * @return true if generator capacity is sufficent to match change in load,
 *         false otherwise
 */
bool gridpack::rtpr::RTPRDriver::adjustRating(double rating, int flag)
{
  bool ret = true;
  int idx;
  double ltotal = p_pf_app.getTotalLoadRealPower(p_dstArea,p_dstZone);
  double gtotal, pmin, pmax;
  p_pf_app.getGeneratorMargins(p_srcArea,p_srcZone,&gtotal,&pmin,&pmax);
  double extra, g_rating;
  if (rating > 1.0) {
    extra = (rating-1.0)*ltotal;
    if (extra <= pmax-gtotal) {
      // Sufficient capacity exists to meet adjustment
      if (pmax > gtotal) {
        g_rating = extra/(pmax-gtotal);
      } else {
        g_rating = 0.0;
      }
    } else {
      extra = pmax-gtotal;
      rating = (ltotal+extra)/ltotal;
      if (pmax > gtotal) {
        g_rating = extra/(pmax-gtotal);
      } else {
        g_rating = 0.0;
      }
      ret = false;
    }
  } else {
    extra = (1.0-rating)*ltotal;
    if (extra <= gtotal-pmin) {
      // Sufficient capacity exists to meet adjustment
      if (gtotal > pmin) {
        g_rating = -extra/(gtotal-pmin);
      } else {
        g_rating = 0.0;
      }
    } else {
      extra = gtotal-pmin;
      rating = (ltotal-extra)/ltotal;
      if (gtotal > pmin) {
        g_rating = -extra/(gtotal-pmin);
      } else {
        g_rating = 0.0;
      }
      ret = false;
    }
  }
  p_pf_app.scaleLoadPower(rating,p_dstArea,p_dstZone);
  p_pf_app.scaleGeneratorRealPower(g_rating,p_srcArea,p_srcZone);
#if 0
  char file[128];
  if (flag == 0) {
    sprintf(file,"pf_diagnostic_%f.dat",rating);
  } else {
    sprintf(file,"ds_diagnostic_%f.dat",rating);
  }
  p_pf_app.writeRTPRDiagnostics(p_srcArea,p_srcZone,p_dstArea,p_dstZone,
      g_rating,rating,file);
#endif
  return ret;
}


/**
 * Execute application. argc and argv are standard runtime parameters
 */
void gridpack::rtpr::RTPRDriver::execute(int argc, char** argv)
{
  int i, j;
  // Get timer instance for timing entire calculation
  gridpack::utility::CoarseTimer *timer =
    gridpack::utility::CoarseTimer::instance();
  int t_total = timer->createCategory("Total Application");
  timer->start(t_total);

  // Read configuration file (user specified, otherwise assume that it is
  // call input.xml)
  gridpack::utility::Configuration *config
    = gridpack::utility::Configuration::configuration();
  if (argc >= 2 && argv[1] != NULL) {
    char inputfile[256];
    sprintf(inputfile,"%s",argv[1]);
    config->open(inputfile,p_world);
  } else {
    config->open("input.xml",p_world);
  }

  // Get size of group (communicator) that individual contingency calculations
  // will run on and create a task communicator. Each process is part of only
  // one task communicator, even though the world communicator is broken up into
  // many task communicators
  gridpack::utility::Configuration::CursorPtr cursor;
  cursor = config->getCursor("Configuration.RealTimePathRating");
  int grp_size;
  // Check to find out if files should be printed for individual power flow
  // calculations
  std::string tmp_bool;
  gridpack::utility::StringUtils util;
  if (!cursor->get("printCalcFiles",&tmp_bool)) {
    p_print_calcs = true;
  } else {
    util.toLower(tmp_bool);
    if (tmp_bool == "false") {
      p_print_calcs = false;
    } else {
      p_print_calcs = true;
    }
  }
  if (!cursor->get("groupSize",&grp_size)) {
    grp_size = 1;
  }
  bool foundArea = true;
  bool found;
  found = cursor->get("sourceArea", &p_srcArea);
  if (!found) {
    p_srcArea = 0;
  }
  foundArea = found && foundArea;
  found = cursor->get("destinationArea", &p_dstArea);
  if (!found) {
    p_dstArea = 0;
  }
  foundArea = found && foundArea;
  if (!cursor->get("sourceZone",&p_srcZone)) {
    p_srcZone = 0;
  }
  if (!cursor->get("destinationZone",&p_dstZone)) {
    p_dstZone = 0;
  }
  bool calcGenCntngcy;
  if (!cursor->get("calculateGeneratorContingencies",&tmp_bool)) {
    calcGenCntngcy = false;
  } else {
    util.toLower(tmp_bool);
    if (tmp_bool == "false") {
      calcGenCntngcy = false;
    } else {
      calcGenCntngcy = true;
    }
  }
  bool calcLineCntngcy;
  if (!cursor->get("calculateLineContingencies",&tmp_bool)) {
    calcLineCntngcy = false;
  } else {
    util.toLower(tmp_bool);
    if (tmp_bool == "false") {
      calcLineCntngcy = false;
    } else {
      calcLineCntngcy = true;
    }
  }
  bool calcTieLines;
  if (!cursor->get("calculateTieLines",&tmp_bool)) {
    calcTieLines = false;
  } else {
    util.toLower(tmp_bool);
    if (tmp_bool == "false") {
      calcTieLines = false;
    } else {
      calcTieLines = true;
    }
  }
  if (!cursor->get("minVoltage",&p_Vmin)) {
    p_Vmin = 0.9;
  }
  if (!cursor->get("maxVoltage",&p_Vmax)) {
    p_Vmax = 1.1;
  }
  // Check for Q limit violations
  if (!cursor->get("checkQLimit",&p_check_Qlim)) {
    p_check_Qlim = false;
  }

  // Monitor generators for frequency violations
  p_monitorGenerators = cursor->get("monitorGenerators",true);
  p_maximumFrequency = cursor->get("frequencyMaximum",61.8);

  // Use branch rating B parameter
  p_useRateB = cursor->get("useBranchRatingB",false);
  if (p_useRateB && p_world.rank() == 0) {
    printf("Using Branch Rating B parameter for checking line overloads\n");
  }

  // TODO: Set these values from input deck
  double start;
  if (!cursor->get("contingencyDSStart",&start)) {
    start = 1.0;
  }
  double end;
  if (!cursor->get("contingencyDSEnd",&end)) {
    end = 1.03;
  }
  double tstep;
  if (!cursor->get("contingencyDSTimeStep",&tstep)) {
    tstep = 0.005;
  }
  if (p_world.rank() == 0) {
    printf("Start time for dynamic simulation contingencies: %f\n",start);
    printf("End time for dynamic simulation contingencies:   %f\n",end);
    printf("Time step for dynamic simulation contingencies:  %f\n",tstep);
  }

  // Create task communicators from world communicator
  p_task_comm = p_world.divide(grp_size);
  // Create powerflow applications on each task communicator
  p_pf_network.reset(new gridpack::powerflow::PFNetwork(p_task_comm));
  // Read in the network from an external file and partition it over the
  // processors in the task communicator. This will read in power flow
  // parameters from the Powerflow block in the input
  p_pf_app.readNetwork(p_pf_network,config);
  // Finish initializing the network
  p_pf_app.initialize();

  if (!calcGenCntngcy && !calcLineCntngcy) {
    // Read in contingency file name
    std::string contingencyfile;
    cursor = config->getCursor(
        "Configuration.RealTimePathRating");
    if (!cursor->get("contingencyList",&contingencyfile)) {
      contingencyfile = "contingencies.xml";
    }
    if (p_world.rank() == 0) {
      printf("Contingency File: %s\n",contingencyfile.c_str());
    }
    // Open contingency file
    bool ok = config->open(contingencyfile,p_world);

    // Get a list of contingencies. Set cursor so that it points to the
    // Contingencies block in the contingency file
    cursor = config->getCursor(
        "ContingencyList.RealTimePathRating.Contingencies");
    if (!cursor) {
      cursor = config->getCursor(
          "ContingencyList.Contingency_analysis.Contingencies");
    }
    gridpack::utility::Configuration::ChildCursors contingencies;
    if (cursor) cursor->children(contingencies);
    p_events = getContingencies(contingencies);
  } else {
    std::vector<gridpack::powerflow::Contingency> genContingencies;
    if (calcGenCntngcy) {
      if (p_world.rank()==0) {
        printf("Calculating generator contingencies automatically"
            " for area %d",p_srcArea);
        if (p_srcZone > 0) {
          printf(" and zone %d\n",p_srcZone);
        } else {
          printf("\n");
        }
      }
      genContingencies = 
        createGeneratorContingencies(p_pf_network, p_srcArea, p_srcZone);
    }
    std::vector<gridpack::powerflow::Contingency> branchContingencies;
    if (calcLineCntngcy) {
      if (p_world.rank()==0) {
        printf("Calculating line contingencies automatically"
            " for area %d",p_srcArea);
        if (p_srcZone > 0) {
          printf(" and zone %d\n",p_srcZone);
        } else {
          printf("\n");
        }
      }
      branchContingencies =
        createBranchContingencies(p_pf_network, p_srcArea, p_srcZone);
    }
    p_events = genContingencies;
    int nsize = branchContingencies.size();
    int ic;
    for (ic = 0; ic<nsize; ic++) {
      p_events.push_back(branchContingencies[ic]);
    }
  }
  // copy power flow continencies to dynamic simulation events
  p_eventsDS.clear();
  for (i=0; i<p_events.size(); i++) {
    gridpack::dynamic_simulation::Event event;
    if (p_events[i].p_type == gridpack::powerflow::Generator) {
      event.from_idx = p_events[i].p_busid[0];
      event.to_idx = p_events[i].p_busid[0];
      event.bus_idx = p_events[i].p_busid[0];
      strncpy(event.tag,p_events[i].p_genid[0].c_str(),2);
      event.tag[2] = '\0';
      event.isGenerator = true;
      event.isLine = false;
    } else {
      event.from_idx = p_events[i].p_from[0];
      event.to_idx = p_events[i].p_to[0];
      strncpy(event.tag,p_events[i].p_ckt[0].c_str(),2);
      event.tag[2] = '\0';
      event.isGenerator = false;
      event.isLine = true;
    }
    event.start = start;
    event.end = end;
    event.step = tstep;
    p_eventsDS.push_back(event);
  }
  // Check for tie lines. Set cursor so that it points to the
  // TieLines block in the input file
  std::vector<gridpack::rtpr::TieLine>  ties;
  if (!calcTieLines) {
    cursor = config->getCursor(
        "Configuration.RealTimePathRating.tieLines");
    gridpack::utility::Configuration::ChildCursors
      tielines;
    if (cursor) cursor->children(tielines);
    ties = getTieLines(tielines);
    // no tie lines found in file but areas (at least) were specified
    // so try calculating tie lines
    if (ties.size() == 0 && foundArea) {
      ties = getTieLines(p_srcArea,p_srcZone,p_dstArea,p_dstZone);
    }
  } else {
    ties = getTieLines(p_srcArea,p_srcZone,p_dstArea,p_dstZone);
  }
  int p_numTies = ties.size();
  p_from_bus.resize(p_numTies);
  p_to_bus.resize(p_numTies);
  p_tags.resize(p_numTies);
  for (i=0; i<p_numTies; i++) {
    p_from_bus[i] = ties[i].from;
    p_to_bus[i] = ties[i].to;
    p_tags[i] = ties[i].tag;
  }
  // Print out list of tie lines
  if (p_world.rank() == 0) {
    printf("List of Tie Lines\n");
    for (i=0; i<p_numTies; i++) {
      printf("From bus: %8d To bus: %8d Line ID: %s\n",p_from_bus[i],
          p_to_bus[i],p_tags[i].c_str());
    }
  }
  // Print out a list of contingencies from process 0
  // (the list is replicated on all processors)
  if (p_world.rank() == 0) {
    int idx;
    for (idx = 0; idx < p_events.size(); idx++) {
      printf("Name: %s\n",p_events[idx].p_name.c_str());
      if (p_events[idx].p_type == Branch) {
        int nlines = p_events[idx].p_from.size();
        int j;
        for (j=0; j<nlines; j++) {
          printf(" Line: (from) %d (to) %d (line) \'%s\'\n",
              p_events[idx].p_from[j],p_events[idx].p_to[j],
              p_events[idx].p_ckt[j].c_str());
        }
      } else if (p_events[idx].p_type == Generator) {
        int nbus = p_events[idx].p_busid.size();
        int j;
        for (j=0; j<nbus; j++) {
          printf(" Generator: (bus) %d (generator ID) \'%s\'\n",
              p_events[idx].p_busid[j],p_events[idx].p_genid[j].c_str());
        }
      }
    }
    if (p_events.size() == 0) {
      printf("***** No contingencies found *****\n");
    }
  }

  p_rating = 1.0;
  bool checkTie = runContingencies();
  if (checkTie) {
    // Tie lines are secure for all contingencies. Increase loads and generation
    while (checkTie) {
      p_rating += 0.05;
      if (!adjustRating(p_rating,0)) {
        if (p_world.rank() == 0) {
          printf("Rating capacity exceeded: %f\n",p_rating);
        }
        p_rating -= 0.05;
        p_pf_app.resetPower();
        break;
      }
      checkTie = runContingencies();
      if (!checkTie) p_rating -= 0.05;
      p_pf_app.resetPower();
    }
    // Refine estimate of rating
    checkTie = true;
    while (checkTie) {
      p_rating += 0.01;
      if (!adjustRating(p_rating,0)) {
        p_rating -= 0.01;
        if (p_world.rank() == 0) {
          printf("Real power generation for power flow"
              " is capacity-limited for Rating: %f\n",
              p_rating);
        }
        p_pf_app.resetPower();
        break;
      }
      checkTie = runContingencies();
      if (!checkTie) p_rating -= 0.01;
      p_pf_app.resetPower();
    }
  } else {
    // Tie lines are insecure for some contingencies. Decrease loads and generation
    while (!checkTie && p_rating >= 0.0) {
      p_rating -= 0.05;
      if (!adjustRating(p_rating,0)) {
        if (p_world.rank() == 0) {
          printf("Rating capacity exceeded: %f\n",p_rating);
        }
        p_rating += 0.05;
        p_pf_app.resetPower();
        break;
      }
      checkTie = runContingencies();
      if (checkTie) p_rating += 0.05;
      p_pf_app.resetPower();
    }
    // Refine estimate of rating
    checkTie = false;
    while (!checkTie && p_rating >= 0.0) {
      p_rating -= 0.01;
      if (!adjustRating(p_rating,0)) {
        p_rating += 0.01;
        if (p_world.rank() == 0) {
          printf("Real power generation for power flow"
              " is capacity-limited for Rating: %f\n",
              p_rating);
        }
        p_pf_app.resetPower();
        break;
      }
      checkTie = runContingencies();
      p_pf_app.resetPower();
    }
  }
  if (p_world.rank() == 0) {
    printf("Final Power Flow Rating: %f\n",p_rating);
  }

  if (p_monitorGenerators) {
    // Power flow estimate of path rating is complete. Now check to see if the
    // rating is stable using dynamic simulation. Start by cloning the dynamic
    // simulation network from the power flow network
    p_ds_network.reset(new gridpack::dynamic_simulation::DSFullNetwork(p_task_comm));
    p_pf_network->clone<gridpack::dynamic_simulation::DSFullBus,
      gridpack::dynamic_simulation::DSFullBranch>(p_ds_network);
    // initialize dynamic simulation
    p_ds_app.setNetwork(p_ds_network, config);
    p_ds_app.readGenerators();
    p_ds_app.initialize();
    p_ds_app.setFrequencyMonitoring(p_monitorGenerators, p_maximumFrequency);

    // find the generators that will be monitored
    findWatchedGenerators(p_ds_network,p_srcArea,p_srcZone,
        p_watch_busIDs,p_watch_genIDs);
    if (p_world.rank() == 0) {
      printf("Generators being monitored in dynamic simulations\n");
      for (i=0; i<p_watch_busIDs.size(); i++) {
        printf("    Host bus: %8d   Generator ID: %s\n",
            p_watch_busIDs[i],p_watch_genIDs[i].c_str());
      }
    }

//    p_ds_app.scaleLoadPower(p_rating,p_dstArea,p_dstZone);
//    p_ds_app.scaleGeneratorRealPower(p_rating,p_srcArea,p_srcZone);
    p_pf_app.resetPower();
    adjustRating(p_rating,1);
    checkTie = runDSContingencies();
    p_pf_app.resetPower();
    p_ds_app.resetPower();

    // Current rating is an upper bound. Only check lower values.
    while (!checkTie && p_rating >= 0.0) {
      p_rating -= 0.05;
      if (!adjustRating(p_rating,1)) {
        p_rating += 0.05;
        p_pf_app.resetPower();
        p_ds_app.resetPower();
        break;
      }
      checkTie = runDSContingencies();
      p_pf_app.resetPower();
      p_ds_app.resetPower();
    }
    // Refine estimate of rating
    checkTie = false;
    while (!checkTie && p_rating >= 0.0) {
      p_rating -= 0.01;
      if (!adjustRating(p_rating,1)) {
        p_rating += 0.01;
        if (p_world.rank() == 0) {
          printf("Real power generation for dynamic simulation"
              " is capacity-limited for Rating: %f\n",
              p_rating);
        }
        p_pf_app.resetPower();
        p_ds_app.resetPower();
        break;
      }
      checkTie = runDSContingencies();
      p_pf_app.resetPower();
      p_ds_app.resetPower();
    }

    if (p_world.rank() == 0) {
      printf("Final Dynamic Simulation Rating: %f\n",p_rating);
    }
  }

  timer->stop(t_total);
  // If all processors executed at least one task, then print out timing
  // statistics (this printout does not work if some processors do not define
  // all timing variables)
  if (p_events.size()*grp_size >= p_world.size()) {
    timer->dump();
  }
}

/**
 * Run complete set of contingencies
 * @return true if no tie line violations found
 */
bool gridpack::rtpr::RTPRDriver::runContingencies()
{
  bool ret = true;
  bool chkLines;
  bool chkSolve;
  // Keep track of failed calculations
  std::vector<int> contingency_idx;
  std::vector<bool> violations;
  std::vector<bool> contingency_success;
  // Keep track of which calculations completed successfully
  gridpack::parallel::GlobalVector<bool> ca_success(p_world);
  // Keep track of violation status for completed calculations
  // 1: no violations
  // 2: voltage violation
  // 3: line overload violation
  // 4: both voltage violation and line overload violation
  std::vector<int> contingency_violation;
  gridpack::parallel::GlobalVector<int> ca_violation(p_world);

  //  Set minimum and maximum voltage limits on all buses
  p_pf_app.setVoltageLimits(p_Vmin, p_Vmax);
  // Solve the base power flow calculation. This calculation is replicated on
  // all task communicators
  char sbuf[128];
  sprintf(sbuf,"base_%f.out",fabs(p_rating)); // absolute function to get rid of
                                              // possible minus sign for zero
  if (p_print_calcs) p_pf_app.open(sbuf);
  sprintf(sbuf,"\nRunning base case on %d processes\n",p_task_comm.size());
  if (p_print_calcs) p_pf_app.writeHeader(sbuf);
  chkSolve = p_pf_app.solve();
  // Check for Qlimit violations
  if (p_check_Qlim && !p_pf_app.checkQlimViolations()) {
    chkSolve = p_pf_app.solve();
  }
  if (!chkSolve) printf("Failed solution on base case\n");
  // Some buses may violate the voltage limits in the base problem. Flag these
  // buses to ignore voltage violations on them.
  p_pf_app.ignoreVoltageViolations();
#ifdef RTPR_DEBUG
  std::vector<std::string> violationDesc = p_pf_app.getContingencyFailures();
#endif

  p_pf_app.useRateB(p_useRateB);
  // Check for tie line violations
  chkLines = p_pf_app.checkLineOverloadViolations(p_from_bus, p_to_bus,
      p_tags, violations);
  ret = ret && chkLines;
  if (!chkLines) printf("Line overload on base case\n");
  // Write out voltages and currents
  if (p_print_calcs) p_pf_app.write();
  if (p_print_calcs) p_pf_app.close();

  // Set up task manager on the world communicator. The number of tasks is
  // equal to the number of contingencies
  gridpack::parallel::TaskManager taskmgr(p_world);
  int ntasks = p_events.size();
  if (ntasks == 0) {
    return chkSolve;
  }
  taskmgr.set(ntasks);
#ifdef USE_STATBLOCK
  gridpack::utility::StringUtils util;
  std::vector<std::string> v_vals = p_pf_app.writeBranchString("flow_str");
  int nsize = v_vals.size();
  std::vector<int> ids;
  std::vector<std::string> tags;
  std::vector<int> mask;
  std::vector<int> id1;
  std::vector<int> id2;
  std::vector<double> pmin, pmax;
  std::vector<double> pflow;
  int i,j;
  for (i=0; i<nsize; i++) {
    std::vector<std::string> tokens = util.blankTokenizer(v_vals[i]);
    if (tokens.size()%8 != 0) {
      printf("Incorrect branch power flow listing\n");
      continue;
    }
    int nline = tokens.size()/8;
    for (j=0; j<nline; j++) {
      id1.push_back(atoi(tokens[j*8].c_str()));
      id2.push_back(atoi(tokens[j*8+1].c_str()));
      tags.push_back(tokens[j*8+2]);
      pflow.push_back(atof(tokens[j*8+3].c_str()));
      pmin.push_back(-atof(tokens[j*8+6].c_str()));
      pmax.push_back(atof(tokens[j*8+6].c_str()));
      if (atoi(tokens[j*8+7].c_str()) == 0) {
        mask.push_back(1);
      } else {
        mask.push_back(2);
      }
    }
  }
  gridpack::analysis::StatBlock pflow_stats(p_world,nsize,ntasks+1);
  if (p_world.rank() == 0) {
    pflow_stats.addRowLabels(id1, id2, tags);
    pflow_stats.addColumnValues(0,pflow,mask);
    pflow_stats.addRowMinValue(pmin);
    pflow_stats.addRowMaxValue(pmax);
  }
#endif

  // Get bus voltage information for base case
  if (p_check_Qlim) p_pf_app.clearQlimViolations();

  // Evaluate contingencies using the task manager
  int task_id;
  // nextTask returns the same task_id on all processors in task_comm. When the
  // calculation runs out of task, nextTask will return false.
  while (taskmgr.nextTask(p_task_comm, &task_id)) {
#ifdef RTPR_DEBUG
    int nsize = violationDesc.size();
    if (task_id == 0 && nsize>0) {
      sprintf(sbuf,"baseCTG_%f.desc",p_rating);
      if (p_task_comm.rank() == 0) {
        int i;
        std::ofstream fout;
        fout.open(sbuf);
        for (i=0; i<nsize; i++) {
          fout << violationDesc[i] << std::endl;
        }
        fout.close();
      }
    }
#endif
    printf("Executing task %d on process %d\n",task_id,p_world.rank());
    sprintf(sbuf,"%s_%f.out",p_events[task_id].p_name.c_str(),fabs(p_rating));
    // Open a new file, based on the contingency name, to store results from
    // this particular contingency calculation
    if (p_print_calcs) p_pf_app.open(sbuf);
    // Write out information to the top of the output file providing some
    // information on the contingency
    sprintf(sbuf,"\nRunning task on %d processes\n",p_task_comm.size());
    if (p_print_calcs) p_pf_app.writeHeader(sbuf);
    if (p_events[task_id].p_type == Branch) {
      int nlines = p_events[task_id].p_from.size();
      int j;
      for (j=0; j<nlines; j++) {
        sprintf(sbuf," Line: (from) %d (to) %d (line) \'%s\'\n",
            p_events[task_id].p_from[j],p_events[task_id].p_to[j],
            p_events[task_id].p_ckt[j].c_str());
        printf("p[%d] Line: (from) %d (to) %d (line) \'%s\'\n",
            p_pf_network->communicator().rank(),
            p_events[task_id].p_from[j],p_events[task_id].p_to[j],
            p_events[task_id].p_ckt[j].c_str());
      }
    } else if (p_events[task_id].p_type == Generator) {
      int nbus = p_events[task_id].p_busid.size();
      int j;
      for (j=0; j<nbus; j++) {
        sprintf(sbuf," Generator: (bus) %d (generator ID) \'%s\'\n",
            p_events[task_id].p_busid[j],p_events[task_id].p_genid[j].c_str());
        printf("p[%d] Generator: (bus) %d (generator ID) \'%s\'\n",
            p_pf_network->communicator().rank(),
            p_events[task_id].p_busid[j],p_events[task_id].p_genid[j].c_str());
      }
    }
    if (p_print_calcs) p_pf_app.writeHeader(sbuf);
    // Reset all voltages back to their original values
    p_pf_app.resetVoltages();
    // Set contingency
    p_pf_app.setContingency(p_events[task_id]);
    // Solve power flow equations for this system
#ifdef USE_SUCCESS
    contingency_idx.push_back(task_id);
#endif
    if (p_pf_app.solve()) {
      chkSolve = true;
#ifdef USE_SUCCESS
      contingency_success.push_back(true);
#endif
      if (p_check_Qlim && !p_pf_app.checkQlimViolations()) {
        chkSolve = p_pf_app.solve();
      }
      if (!chkSolve) printf("Failed solution on continency %d\n",task_id+1);
      // If power flow solution is successful, write out voltages and currents
      if (p_print_calcs) p_pf_app.write();
      chkLines = p_pf_app.checkLineOverloadViolations(p_from_bus, p_to_bus, p_tags,
          violations);
      ret = ret && chkLines;
      // Check for violations
      bool ok1 = p_pf_app.checkVoltageViolations();
      bool ok2 = p_pf_app.checkLineOverloadViolations();
      bool ok = ok1 && ok2;
      // Include results of violation checks in output
      if (ok) {
        sprintf(sbuf,"\nNo violation for contingency %s\n",
            p_events[task_id].p_name.c_str());
#ifdef USE_SUCCESS
        contingency_violation.push_back(1);
#endif
      } 
      if (!ok1) {
        sprintf(sbuf,"\nBus Violation for contingency %s\n",
            p_events[task_id].p_name.c_str());
      }
      if (p_print_calcs) p_pf_app.print(sbuf);
      if (p_print_calcs) p_pf_app.writeCABus();
      if (!ok2) {
        sprintf(sbuf,"\nBranch Violation for contingency %s\n",
            p_events[task_id].p_name.c_str());
      }

#ifdef USE_SUCCESS
      if (!ok1 && !ok2) {
        contingency_violation.push_back(4);
      } else if (!ok1) {
        contingency_violation.push_back(2);
      } else if (!ok2) {
        contingency_violation.push_back(3);
      }
#endif
#ifdef USE_STATBLOCK
      pflow.clear();
      mask.clear();
      v_vals.clear();
      v_vals = p_pf_app.writeBranchString("flow_str");
      nsize = v_vals.size();
      for (i=0; i<nsize; i++) {
        std::vector<std::string> tokens = util.blankTokenizer(v_vals[i]);
        if (tokens.size()%8 != 0) {
          printf("Incorrect branch power flow listing\n");
          continue;
        }
        int nline = tokens.size()/8;
        for (j=0; j<nline; j++) {
          pflow.push_back(atof(tokens[j*8+3].c_str()));
          if (atoi(tokens[j*8+7].c_str()) == 0) {
            mask.push_back(1);
          } else {
            mask.push_back(2);
          }
        }
      }
      if (p_task_comm.rank() == 0) {
        pflow_stats.addColumnValues(task_id+1,pflow,mask);
      }
#endif

#ifdef RTPR_DEBUG
      violationDesc = p_pf_app.getContingencyFailures();
      nsize = violationDesc.size();
      if (nsize>0) {
        sprintf(sbuf,"%s_%f.desc",p_events[task_id].p_name.c_str(),p_rating);
        if (p_task_comm.rank() == 0) {
          int i;
          std::ofstream fout;
          fout.open(sbuf);
          for (i=0; i<nsize; i++) {
            fout << violationDesc[i] << std::endl;
          }
          fout.close();
        }
      }
#endif
      if (p_print_calcs) p_pf_app.print(sbuf);
      if (p_print_calcs) p_pf_app.writeCABranch();
      if (p_check_Qlim) p_pf_app.clearQlimViolations();
    } else {
      if (!chkSolve) printf("Failed solution on continency %d\n",task_id+1);
#ifdef USE_SUCCESS
      contingency_success.push_back(false);
      contingency_violation.push_back(0);
#endif
      sprintf(sbuf,"\nDivergent for contingency %s\n",
          p_events[task_id].p_name.c_str());
      if (p_print_calcs) p_pf_app.print(sbuf);
#ifdef USE_STATBLOCK
      pflow.clear();
      mask.clear();
      v_vals.clear();
      v_vals = p_pf_app.writeBranchString("fail_str");
      nsize = v_vals.size();
      for (i=0; i<nsize; i++) {
        std::vector<std::string> tokens = util.blankTokenizer(v_vals[i]);
        if (tokens.size()%8 != 0) {
          printf("Incorrect branch power flow listing\n");
          continue;
        }
        int nline = tokens.size()/8;
        for (j=0; j<nline; j++) {
          pflow.push_back(0.0);
          mask.push_back(0);
        }
      }
      if (p_task_comm.rank() == 0) {
        pflow_stats.addColumnValues(task_id+1,pflow,mask);
      }
#endif
    } 
    // Return network to its original base case state
    p_pf_app.unSetContingency(p_events[task_id]);
    // Close output file for this contingency
    if (p_print_calcs) p_pf_app.close();
  }
  // Print statistics from task manager describing the number of tasks performed
  // per processor
  taskmgr.printStats();

  // Gather stats on successful contingency calculations
#ifdef USE_SUCCESS
  if (p_task_comm.rank() == 0) {
    ca_success.addElements(contingency_idx, contingency_success);
    ca_violation.addElements(contingency_idx, contingency_violation);
  }
  ca_success.upload();
  ca_violation.upload();
  // Write out stats on successful calculations
  if (p_world.rank() == 0) {
    char sbuf[128];
    contingency_idx.clear();
    contingency_success.clear();
    contingency_violation.clear();
    int i;
    for (i=0; i<ntasks; i++) contingency_idx.push_back(i);
    ca_success.getData(contingency_idx, contingency_success);
    contingency_success.clear();
    ca_violation.getData(contingency_idx, contingency_violation);
    std::ofstream fout;
    sprintf(sbuf,"success_%f.dat",p_rating);
    fout.open(sbuf);
    for (i=0; i<ntasks; i++) {
      if (contingency_success[i]) {
        fout << "contingency: " << p_events[i].p_name << " success: true";
        if (contingency_violation[i] == 1) {
          fout << " violation: none" << std::endl;
        } else if (contingency_violation[i] == 2) {
          fout << " violation: bus" << std::endl;
        } else if (contingency_violation[i] == 3) {
          fout << " violation: branch" << std::endl;
        } else if (contingency_violation[i] == 4) {
          fout << " violation: bus and branch" << std::endl;
        }
      } else {
        fout << "contingency: " << p_events[i].p_name
          << " success: false" << std::endl;
      }
    }
    fout.close();
  }
#endif
#ifdef USE_STATBLOCK
  sprintf(sbuf,"line_flt_cnt_%f.txt",p_rating);
  pflow_stats.writeMaskValueCount(sbuf,2);
#endif
  // Check to see if ret is true on all processors
  int ok;
  if (ret) {
    ok = 1;
  } else {
    ok = 0;
  }
  p_world.sum(&ok,1);
  if (ok == p_world.size()) {
    ret = true;
  } else {
    ret = false;
  }
  return ret;
}

/**
 * Transfer data from power flow to dynamic simulation
 * @param pf_network power flow network
 * @param ds_network dynamic simulation network
 */
void gridpack::rtpr::RTPRDriver::transferPFtoDS(
    boost::shared_ptr<gridpack::powerflow::PFNetwork> pf_network,
    boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork>
    ds_network)
{
  int numBus = pf_network->numBuses();
  int i;
  gridpack::component::DataCollection *pfData;
  gridpack::component::DataCollection *dsData;
  double rval;
  for (i=0; i<numBus; i++) {
    pfData = pf_network->getBusData(i).get();
    dsData = ds_network->getBusData(i).get();
    pfData->getValue("BUS_PF_VMAG",&rval);
    dsData->setValue(BUS_VOLTAGE_MAG,rval);
    //printf("Step0 bus%d mag = %f\n", i+1,rval);
    pfData->getValue("BUS_PF_VANG",&rval);
    dsData->setValue(BUS_VOLTAGE_ANG,rval);
    int ngen = 0;
    if (pfData->getValue(GENERATOR_NUMBER, &ngen)) {
      int j;
      for (j=0; j<ngen; j++) {
        pfData->getValue("GENERATOR_PF_PGEN",&rval,j);
        dsData->setValue(GENERATOR_PG,rval,j);
        //printf("p[%d] save PGEN: %f\n", rval,p_world.rank());
        pfData->getValue("GENERATOR_PF_QGEN",&rval,j);
        dsData->setValue(GENERATOR_QG,rval,j);
        //printf("p[%d] save QGEN: %f\n", rval,p_world.rank());
      }
    }
    gridpack::powerflow::PFBus *bus = pf_network->getBus(i).get();
    int izone;
    if (p_dstZone > 0) {
      izone = bus->getZone();
    } else {
      izone = p_dstZone;
    }
    if (bus->getArea() == p_dstArea && p_dstZone == izone) {
      int nld = 0;
      if (pfData->getValue(LOAD_NUMBER, &ngen)) {
        std::vector<double> pl, ql;
        std::vector<int> status;
        std::vector<std::string> tags;
        bus->getLoadPower(tags,pl,ql,status);
        if (tags.size() != ngen) {
          printf("Mismatch in loads on bus %d\n",
              bus->getOriginalIndex());
          continue;
        }
        int j;
        for (j=0; j<ngen; j++) {
          dsData->setValue(LOAD_PL,p_rating*pl[j],j);
          dsData->setValue(LOAD_QL,p_rating*ql[j],j);
        }
      }
    }
  }
}


/**
 * Run dynamic simulations over full set of contingencies
 * @return true if no violations found on complete set of contingencies
 */
bool gridpack::rtpr::RTPRDriver::runDSContingencies()
{
  bool ret = true;
  // Set up task manager on the world communicator. The number of tasks is
  // equal to the number of contingencies
  gridpack::parallel::TaskManager taskmgr(p_world);
  int ntasks = p_eventsDS.size();
  taskmgr.set(ntasks);

  // Evaluate contingencies using the task manager
  int task_id;
  // nextTask returns the same task_id on all processors in task_comm. When the
  // calculation runs out of task, nextTask will return false.
#if 0
  char file[128];
  sprintf(file,"pf2_diagnostic_%f_%d.dat",p_rating,p_world.rank());
  p_pf_app.writeRTPRDiagnostics(p_srcArea,p_srcZone,p_dstArea,p_dstZone,
      p_rating,p_rating,file);
#endif
  bool chkSolve = p_pf_app.solve();
  p_pf_app.useRateB(p_useRateB);
  // Check for Qlimit violations
  if (p_check_Qlim && !p_pf_app.checkQlimViolations()) {
    chkSolve = p_pf_app.solve();
  }
  p_pf_app.saveData();
#if 0
  sprintf(file,"pf3_diagnostic_%f_%d.dat",p_rating,p_world.rank());
  p_pf_app.writeRTPRDiagnostics(p_srcArea,p_srcZone,p_dstArea,p_dstZone,
      p_rating,p_rating,file);
#endif
  while (taskmgr.nextTask(p_task_comm, &task_id)) {
#ifdef RTPR_DEBUG
    if (task_id == 0) {
      char sbuf[128];
      sprintf(sbuf,"baseDS_%f.out",fabs(p_rating)); // absolute function to get rid of
      // possible minus sign for zero
      p_pf_app.open(sbuf);
      sprintf(sbuf,"\nRunning initialization for Rating: %f\n",fabs(p_rating));
      p_pf_app.writeHeader(sbuf);
      p_pf_app.write();
      p_pf_app.close();
    }
#endif
    // Print out results of power flow calculation
    printf("Executing dynamic simulation task %d on process %d\n",
        task_id,p_world.rank());
    // reinitialize dynamic simulation from powerflow calculation
    transferPFtoDS(p_pf_network,p_ds_network);
    p_ds_app.reload();
    p_ds_app.setGeneratorWatch(p_watch_busIDs,p_watch_genIDs,false);
    try {
      p_ds_app.solve(p_eventsDS[task_id]);
    } catch (const gridpack::Exception e) {
      printf("Failed to execute DS task %d on process %d\n",
          task_id,p_world.rank());
    }
#ifdef RTPR_DEBUG
    if (!p_ds_app.frequencyOK()) {
      ret = false;
      std::vector<int> violations = p_ds_app.getFrequencyFailures();
      int nsize = violations.size();
      char sbuf[128];
      if (nsize>0) {
        sprintf(sbuf,"DS_%s_%f.desc",p_events[task_id].p_name.c_str(),p_rating);
        if (p_task_comm.rank() == 0) {
          int i;
          std::ofstream fout;
          fout.open(sbuf);
          for (i=0; i<nsize; i++) {
            fout << "Frequency violation on bus "<<violations[i] << std::endl;
          }
          fout.close();
        }
      }
    }
#else
    ret = ret && p_ds_app.frequencyOK();
#endif
  
  }
  int iret = static_cast<int>(ret);
  p_world.sum(&iret,1);
  if (iret == p_world.size()) {
    ret = true;
  } else {
    ret = false;
  }
  return ret;
}
