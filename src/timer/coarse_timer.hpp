/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#ifndef _coarse_timer_h
#define _coarse_timer_h

#include <map>
#include <string>
#include <vector>

#include <boost/serialization/export.hpp>

// Simple outline of data collection object

namespace gridpack{
namespace utility{

class CoarseTimer {
public:
  
  /**
   * Retrieve instance of the CoarseTimer object
   */
  static CoarseTimer *instance();

  /**
   * Create a new timer category and return a handle to the category. It is up
   * to the application to keep track of this handle.
   * @param title the title is the name that will be used to label the timing
   *        statistics in the output
   * @return an integer handle that can be used to refer to this category
   */
  int createCategory(const std::string title);

  /**
   * Start timing the category
   * @param idx category handle
   */
  void start(const int idx);

  /**
   * Stop timing the category
   * @param idx category handle
   */
  void stop(const int idx);

  /**
   * Write all timing statistics to standard out
   */
  void dump(void) const;

  /**
   * Write out profile for all processors for requested handle
   * @param idx category handle
   */
  void dumpProfile(const int idx) const;

  /**
   * Write out profile for all processors for requested title
   * @param idx category title
   */
  void dumpProfile(const std::string title) const;

  /**
   * Return current time. Can be used to solve timing problems that can't be
   * handled using the regular timing capabilities
   * @return current time in seconds according to internal clock
   */
  double currentTime();

  /**
   * Turn timing on and off. If timing is off, no data is collected.
   * @param flag turn timer on (true) or off (false)
   */
  void configTimer(bool flag);

protected:
  /**
   * Constructor
   */
  CoarseTimer();

  /**
   * Destructor
   */
  ~CoarseTimer();

private:

  std::map<std::string, int> p_title_map; 
  std::vector<std::string> p_title;
  std::vector<double> p_start;
  std::vector<double> p_time;
  std::vector<int>    p_istart;
  std::vector<int>    p_istop;

  static CoarseTimer *p_instance;

  bool                p_profile;
};


}    // utility
}    // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::utility::CoarseTimer)

#endif // _coarse_timer_h
