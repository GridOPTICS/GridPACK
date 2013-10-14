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
  int createCategory(string title);

  /**
   * Start timing the category
   * @param idx category handle
   */
  void start(int idx);

  /**
   * Stop timing the category
   * @param idx category handle
   */
  void stop(int idx);

  /**
   * Write all timing statistics to standard out
   */
  void dump(void);

protected
  /**
   * Constructor
   */
  CoarseTimer();

  /**
   * Destructor
   */
  CoarseTimer();

private:

  std::map<std::string, int> p_title_map; 
  std::vector<std::string> p_title;
  std::vector<double> p_start;
  std::vector<double> p_time;
  std::vector<int>    p_istart;
  std::vector<int>    p_istop;

  static CoarseTimer *p_instance;
};


}    // utility
}    // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::utility::CoarseTimer);

#endif // _coarse_timer_h
