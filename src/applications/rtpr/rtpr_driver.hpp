/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   rtpr_driver.hpp
 * @author Bruce Palmer
 * @date   October 9, 2019
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------

#ifndef _rtpr_driver_h_
#define _rtpr_driver_h_

#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/modules/powerflow/pf_app_module.hpp"

namespace gridpack {
namespace rtpr {

enum ContingencyType{Generator, Branch};

/* Defininition of contingency data structure (from powerflow module)
struct Contingency
{
  int p_type;
  std::string p_name;
  // Line contingencies
  std::vector<int> p_from;
  std::vector<int> p_to;
  std::vector<std::string> p_ckt;
  // Status of line before contingency
  std::vector<bool> p_saveLineStatus;
  // Generator contingencies
  std::vector<int> p_busid;
  std::vector<std::string> p_genid;
  // Status of generator before contingency
  std::vector<bool> p_saveGenStatus;
};
*/

struct TieLine {
  int from;
  int to;
  char tag[3];
};

// Calling program for real-time path rating application
class RTPRDriver
{
  public:
    /**
     * Basic constructor
     */
    RTPRDriver(void);

    /**
     * Basic destructor
     */
    ~RTPRDriver(void);

    /**
     * Get list of contingencies from external file
     * @param cursor pointer to contingencies in input deck
     * @return vector of contingencies
     */
    std::vector<gridpack::powerflow::Contingency> getContingencies(
        gridpack::utility::Configuration::ChildCursors &contingencies);

    /**
     * Create a list of all N-1 generator contingencies for a given area and
     * zone
     * @param network power grid network on which contingencies are defined
     * @param area index of area that will generate contingencies
     * @param zone index of zone that will generate contingencies
     * @return vector of contingencies
     */
    std::vector<gridpack::powerflow::Contingency> createGeneratorContingencies(
        boost::shared_ptr<gridpack::powerflow::PFNetwork> network, int area,
        int zone);

    /**
     * Create a list of all N-1 branch contingencies for a given area and
     * zone
     * @param network power grid network on which contingencies are defined
     * @param area index of area that will generate contingencies
     * @param zone index of zone that will generate contingencies
     * @return vector of contingencies
     */
    std::vector<gridpack::powerflow::Contingency> createBranchContingencies(
        boost::shared_ptr<gridpack::powerflow::PFNetwork> network, int area,
        int zone);

    /**
     * Get list of tie lines
     * @param cursor pointer to tie lines in input deck
     * @return vector of contingencies
     */
    std::vector<gridpack::rtpr::TieLine> getTieLines(
          gridpack::utility::Configuration::ChildCursors &tielines);

    /**
     * Automatically create list of tie lines between area1,zone1 and
     * area2,zone2
     * @param area1 index of source area
     * @param zone1 index of source zone
     * @param area2 index of destination area
     * @param zone2 index of destination zone
     * @param tielines list of tie lines between area1,zone1 and area2,zone2
     */
    std::vector<gridpack::rtpr::TieLine> getTieLines(int area1, int zone1,
        int area2, int zone2);

    /**
     * Create a list of generator for a given area and zone. These
     * generators will be monitored to see if they exceed operating
     * specifications.
     * @param network power grid used for DS simulation
     * @param area index of area that will generate contingencies
     * @param zone index of zone that will generate contingencies
     * @param buses list of bus IDs that contain generators
     * @param tags list of generator IDs
     */
    void  findWatchedGenerators(
        boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork> network,
        int area, int zone, std::vector<int> &buses,
        std::vector<std::string> &tags);

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
    bool adjustRating(double rating, int flag);

    /**
     * Execute application
     * @param argc number of arguments
     * @param argv list of character strings
     */
    void execute(int argc, char** argv);

    /**
     * Run complete set of contingencies
     * @return true if no tie line violations found
     */
    bool runContingencies();

    /**
     * Run dynamic simulations over full set of contingencies
     * @return true if no violations found on complete set of contingencies
     */
    bool runDSContingencies();

    /**
     * Transfer data from power flow to dynamic simulation
     * @param pf_network power flow network
     * @param ds_network dynamic simulation network
     */
    void transferPFtoDS(
        boost::shared_ptr<gridpack::powerflow::PFNetwork> pf_network,
        boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork>
        ds_network);

    private:

    boost::shared_ptr<gridpack::powerflow::PFNetwork> p_pf_network;

    boost::shared_ptr<gridpack::dynamic_simulation::DSFullNetwork>
      p_ds_network;

    gridpack::powerflow::PFAppModule p_pf_app;
    gridpack::dynamic_simulation::DSFullApp p_ds_app;

    std::vector<gridpack::powerflow::Contingency> p_events;
    std::vector<gridpack::dynamic_simulation::Event> p_eventsDS;

    int p_srcArea, p_dstArea, p_srcZone, p_dstZone;

    double p_rating;

    double p_Vmin, p_Vmax;
    
    bool p_monitorGenerators;

    double p_maximumFrequency;

    bool p_check_Qlim, p_print_calcs;

    std::vector<int> p_from_bus, p_to_bus;

    std::vector<std::string> p_tags;

    std::vector<bool> p_violations;

    gridpack::parallel::Communicator p_world;
    gridpack::parallel::Communicator p_task_comm;

    int p_numTies;

    bool p_useRateB;

    std::vector<int> p_watch_busIDs;
    std::vector<std::string> p_watch_genIDs;

};

} // rtpr
} // gridpack
#endif
