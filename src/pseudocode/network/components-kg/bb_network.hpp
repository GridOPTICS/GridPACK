// -------------------------------------------------------------
/**
 * @file   network.hpp
 * @author Kevin A. Glass
 * @date   Fri Apr  19 13:36:28 2013
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  19, 2013 by Kevin A. Glass
// Last Change:
// -------------------------------------------------------------

#ifndef _pf_network_hpp_
#define _pf_network_hpp_

namespace gridpack {
namespace network {

class MatrixInterface;
class MapConstructor;
class MapSize;
class MapData;
#include <iostream>
#include <vector>

// -------------------------------------------------------------
//  class BBNetwork
// -------------------------------------------------------------
/**
 * 
 */

template <typename BUS, typename BRANCH, typename MEASUREMENT>
class BBNetwork : public Network
{
public:
  BBNetwork(){};
  virtual ~BBNetwork(void){};

  // The factory creates a branch object and passes it to the BBNetwork
  virtual void addBranch(BRANCH  * branch)
  {
    branches_.push_back(new Component<BRANCH>(branch));
  }

  virtual void addBuses(BUS  * bus) {
    buses_.push_back(new Component<BUS>(bus));
  }

  virtual void addMeasurement(MEASUREMENT * measurement) {
    measurements_.push_back(new Component<MEASUREMENT>(measurement));
  }

  // this constructs the objects required by the analysis
  void networkMap(MapConstructor * map) {
    networkMap_(map);
  }
protected:
  virtual void networkMap_(MapConstructor * size){};
  virtual void networkMap_(MapSize * size)
  {
    for () {
      buses_[i]->setAnalysisData(size)
    }
    for () {
      branches_[i]->setAnalysisData(size)
    }
  }

  virtual void networkMap_(MapData * data)
  {
    for () {
      buses_[i]->setAnalysisData(data)
    }
    for () {
      branches_[i]->setAnalysisData(data)
    }
  }

private:
  std::vector<Component<BUS> * >        buses_;
  std::vector<Component<BRANCH> * >     branches_;
  std::vector<Component<MEASUREMENT> *> measurements_;

  int connection[];
};

} // namespace math
} // namespace gridpack

#endif
