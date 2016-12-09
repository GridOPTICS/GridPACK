// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   line_wrapping_output_filter.hpp
 * @author William A. Perkins
 * @date   2016-12-07 13:32:29 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December  7, 2016 by William A. Perkins
// Last Change: 2013-05-03 12:23:12 d3g096
// -------------------------------------------------------------


#ifndef _line_wrapping_output_filter_hpp_
#define _line_wrapping_output_filter_hpp_

#include <boost/iostreams/filtering_stream.hpp>

namespace io = boost::iostreams;

#include <boost/iostreams/filtering_stream.hpp>

// -------------------------------------------------------------
// line_wrapping_output_filter
//
// slightly modified from the Boost iostreams example
// -------------------------------------------------------------
class line_wrapping_output_filter : public io::output_filter {
public:
  explicit line_wrapping_output_filter(int line_length = 80, int margin = 8)
    : do_new_line_(false), line_length_(line_length), end_margin_(margin), col_no_(0) 
  { }

  template<typename Sink>
  bool put(Sink& dest, int c)
  {
    // if the column is past the goal end of line, find the next white space, eat it up, then
    if (do_new_line_) {
      if (c == '\n') {
        do_new_line_ = false;
      } else if (!std::isspace(c)) {
        if (!put_char(dest, '\n')) return false;
        do_new_line_ = false;
      } else {
        // eating white space
        return true;
      }
    } else if (col_no_ >= (line_length_ - end_margin_)) {
      // wait until we see some white space before doing a new line
      if (std::isspace(c)) {
        do_new_line_ = true;
      }
    }
    return put_char(dest, c);
  }

  template<typename Sink>
  void close(Sink&)
  { col_no_ = 0; }

private:
  template<typename Sink>
  bool put_char(Sink& dest, int c)
  {
    if (!io::put(dest, c))
      return false;
    if (c != '\n')
      ++col_no_;
    else
      col_no_ = 0;
    return true;
  }
  int  do_new_line_;
  int  line_length_;
  int  end_margin_;
  int  col_no_;
};

#endif
