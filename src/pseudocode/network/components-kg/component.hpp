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

#ifndef _component_hpp_
#define _component_hpp_

namespace gridpack {
namespace network {

#include <iostream>
#include <vector>

// -------------------------------------------------------------
//  class Network
// -------------------------------------------------------------
/**
 * 
 * To be safe, these should be used simultaneously on all processes.  
 */
class MapConstructor;
class MapSize;
class MapData;

template <typename COMPONENT>
class Component
{
public:
  Component<COMPONENT>(COMPONENT * component) : component_(component){};
  virtual ~Component<COMPONENT>(void){};

  void setAnalysisData (MapConstructor * map){setAnalysisData_(map);};
protected:
private:
  void setAnalysisData_(MapConstructor * map){/* throw an exception */};
  void setAnalysisData_(MapSize * map){component_->increment(map);};
  void setAnalysisData_(MapData * map){component_->mapData(map);};

  COMPONENT             * component_;
};

} // namespace math
} // namespace gridpack

#endif
