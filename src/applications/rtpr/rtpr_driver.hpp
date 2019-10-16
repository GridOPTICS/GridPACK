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
     * Create a list of all N-1 generator contingencies for a given area
     * @param network power grid network on which contingencies are defined
     * @param area index of area that will generate contingencies
     * @return vector of contingencies
     */
    std::vector<gridpack::powerflow::Contingency> createGeneratorContingencies(
        boost::shared_ptr<gridpack::powerflow::PFNetwork> network, int area);

    /**
     * Create a list of all N-1 branch contingencies for a given area
     * @param network power grid network on which contingencies are defined
     * @param area index of area that will generate contingencies
     * @return vector of contingencies
     */
    std::vector<gridpack::powerflow::Contingency> createBranchContingencies(
        boost::shared_ptr<gridpack::powerflow::PFNetwork> network, int area);

    /**
     * Get list of tie lines
     * @param cursor pointer to tie lines in input deck
     * @return vector of contingencies
     */
    std::vector<gridpack::rtpr::TieLine> getTieLines(
          gridpack::utility::Configuration::ChildCursors &tielines);

    /**
     * Scale generation in a specified area
     * @param scale value to scale real power generation
     * @param area index of area
     */
    void scaleAreaGeneration(double scale, int area);

    /**
     * Scale loads in a specified area
     * @param scale value to scale real power load
     * @param area index of area
     */
    void scaleAreaLoads(double scale, int area);

    /**
     * Execute application
     * @param argc number of arguments
     * @param argv list of character strings
     */
    void execute(int argc, char** argv);

    private:

    boost::shared_ptr<gridpack::powerflow::PFNetwork> p_pf_network;
};

} // rtpr
} // gridpack
#endif
