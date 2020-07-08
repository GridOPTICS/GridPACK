// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dae_event.hpp
 * @author Perkins
 * @date   2020-07-08 12:45:52 d3g096
 * 
 * @brief  Encapsulate an "Event" that would affect a time integration problem
 * 
 * This is based primarily on PETSc's TSEvent model.  It may not be
 * general enough.  
 * 
 */
// -------------------------------------------------------------
// Created November 18, 2019 by Perkins
// Last Change: 2018-07-24 09:27:04 d3g096
// -------------------------------------------------------------


#ifndef _dae_event_hpp_
#define _dae_event_hpp_

#include <vector>
#include <list>
#include <boost/scoped_array.hpp>
#include <boost/foreach.hpp>
#include <gridpack/math/vector.hpp>
#include <gridpack/utilities/uncopyable.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
// enum DAEEventDirection
// -------------------------------------------------------------
enum DAEEventDirection
  {
   CrossZeroNegative = -1,
   CrossZeroPositive = 1,
   CrossZeroEither = 0
  };



// -------------------------------------------------------------
// abstract class DAEEventT
// -------------------------------------------------------------
/// Description of zero events and how to handle them
/**
 * This type knows how to compute one or more event values and adjust
 * the integration problem when any of those values are triggered
 * 
 */

template <typename T, typename I = int>
class DAEEventT
  : private utility::Uncopyable
{
public:
  /// Default constructor.
  /** 
   * A single instance can handle multiple "events". 
   * 
   * @param n number of events handled by this instance
   * @param dirs the directions of each event
   * 
   * @return new instance
   */
  DAEEventT(const int& nevent)
    : utility::Uncopyable(),
      p_stateIndex(0),
      p_size(nevent), 
      p_dir(p_size, CrossZeroEither), // filled by children
      p_term(p_size, false),          // filled by children, but all false initially
      p_current(p_size, -9999.0)      // filled by children, not zero initially
  {
  }

  /// Destructor
  virtual ~DAEEventT(void)
  {}

  /// get the size of this event
  int size(void) const
  {
    return p_size;
  }

  /// get the required index into the @em local state vector
  int stateIndex(void) const
  {
    return p_stateIndex;
  }

  /// get the directions of this event
  const std::vector<DAEEventDirection>& dirs(void) const
  {
    return p_dir;
  }

  /// get the termination flags for this instance
  const std::vector<bool>& termFlags(void) const
  {
    return p_term;
  }

  /// get the current values
  const std::vector<T>& values(void) const
  {
    return p_current;
  }

  /// update and return event values, given the state
  const std::vector<T>& values(const double& t, T *state)
  {
    this->p_update(t, state);
    return this->values();
  }

  /// handle any triggered events
  void handle(const bool *triggered, const double& t, T *state)
  {
    this->p_handle(triggered, t, state);
  }

protected:

  /// The number of event values handled by this instance
  const int p_size;

  /// The base index into the @em local solution index
  int p_stateIndex;

  /// The direction which event value passes zero
  std::vector<DAEEventDirection> p_dir;

  /// A flag indicating which event values terminate the solver
  std::vector<bool> p_term;

  /// Current event values
  std::vector<T> p_current;
  
  /// update and return event values, given the state (specialized)
  /** 
   * Fill #p_current with values base on the current #state
   * 
   * @param t 
   * @param T 
   * @param state 
   */
  virtual void p_update(const double& t, T *state) = 0;

  /// handle any triggered events  (specialized)
  /** 
   * This does whatever is necessary to adjust the internals of this
   * instance based on the events that were triggered.
   * 
   * @param triggered array of length #p_size that indicate which events occurred
   * @param t solver time
   * @param state current solver state
   */
  virtual void p_handle(const bool *triggered, const double& t, T *state) = 0;
};

// -------------------------------------------------------------
// class DAEEventManagerT
// -------------------------------------------------------------
template <typename T, typename I = int>
class DAEEventManagerT
  : private utility::Uncopyable
{
public:

  typedef DAEEventT< T, I > Event; 
  typedef boost::shared_ptr< Event > EventPtr;

  /// Default constructor.
  DAEEventManagerT(void)
    : utility::Uncopyable(),
      p_events(), p_size(), p_dirs(), p_terms(), p_trigger(), p_terminated()
  {}

  /// Destructor
  virtual ~DAEEventManagerT(void)
  {}

  /// Add an event
  void add(EventPtr e)
  {
    int nextidx;
    if (!p_events.empty()) {
      nextidx = p_events.back().evidx;
      nextidx += p_events.back().nval;
    } else {
      nextidx = 0;
    }
    
    p_EventInfo rec;
    rec.event = e;
    rec.nval = e->size();
    rec.evidx = nextidx;
    rec.solidx = e->stateIndex();
    p_events.push_back(rec);
    p_size = -1;
    p_dirs.reset();
    p_terms.reset();
    p_values.reset();
    p_trigger.reset();
  }

  /// Get the number of event values handled
  int size(void) const
  {
    int n(0);
    if (p_size >= 0) {
      n = p_size;
    } else {
      BOOST_FOREACH(const p_EventInfo rec, p_events) {
        n += rec.nval;
      }
      p_size = n;
    }
    return n;
  }

  /// Get all the event directions
  DAEEventDirection const *
  directions(void) const
  {
    int n(this->size());
    if (!p_dirs) {
      p_dirs.reset(new DAEEventDirection[n]);
    }

    n = 0;
    BOOST_FOREACH(const p_EventInfo rec, p_events) {
      const std::vector<DAEEventDirection> edirs((rec.event)->dirs());
      BOOST_FOREACH(const DAEEventDirection d, edirs) {
        p_dirs[n++] = d;
      }
    }
    return p_dirs.get();
  }

  bool const *
  terminateFlags(void) const
  {
    int n;
    if (!p_terms) {
      n = this->size();
      p_terms.reset(new bool[n]);
    }

    n = 0;
    BOOST_FOREACH(const p_EventInfo rec, p_events) {
      const std::vector<bool>& eterm((rec.event)->termFlags());
      BOOST_FOREACH(const bool t, eterm) {
        p_terms[n++] = t;
      }
    }
    return p_terms.get();
  }


  /// Update and get the current event values
  T const *
  values(const double t, VectorT<T, I>& state)
  {
    int n;
    if (!p_values) {
      n = this->size();
      p_values.reset(new T[n]);
    }

    // Get the local part of the state vector into a plain arrary
    typename VectorT<T, I>::IdxType lo, hi;
    state.localIndexRange(lo, hi);
    boost::scoped_array<T> lstate(new T[hi-lo]);
    state.getElementRange(lo, hi, &lstate[0]);

    n = 0;
    BOOST_FOREACH(const p_EventInfo rec, p_events) {
      int lidx(rec.solidx);
      const std::vector<T>& evals((rec.event)->values(t, &lstate[lidx]));
      BOOST_FOREACH(const T v, evals) {
        p_values[n++] = v;
      }
    }
    return p_values.get();
  }

  /// Handle triggered events
  void
  handle(const int nevent, int *eventidx, const double t, VectorT<T, I>& state)
  {
    if (!p_trigger) {
      p_trigger.reset(new bool[this->size()]);
    }
    for (int i = 0; i < this->size(); ++i) {
      p_trigger[i] = false;
    }
    for (int e = 0; e < nevent; ++e) {
      p_trigger[eventidx[e]] = true;
    }

    const bool *term = this->terminateFlags();

    T *lstate = state.getLocalElements();

    BOOST_FOREACH(p_EventInfo rec, p_events) {
      for (int i = 0; i < rec.nval; ++i) {
        int idx(rec.evidx + i);
        if (p_trigger[idx]) {
          rec.event->handle(&p_trigger[rec.evidx], t, &lstate[rec.evidx]);
          if (term[idx]) { p_terminated = true; }
          break;
        }
      }
    }
    state.releaseLocalElements(lstate);
    return;
  }

  /// Has a termination event occurred
  bool terminated(void) const
  {
    return p_terminated;
  }

  /// (Re)Set the termination event flag
  void terminated(const bool& flag)
  {
    p_terminated = flag;
  }

  
  
protected:

  /// A thing to hold event information
  struct p_EventInfo {
    EventPtr event;
    int evidx;                  /**< Index of first event value in this record */
    int nval;                   /**< The number of values generated by event */
    I solidx;                   /**< The base index into the (local) solution vector */
  };

  /// All the DAEEventT instances managed by this instance
  std::list<p_EventInfo> p_events;

  /// Cache the overall number of event values
  mutable int p_size;

  /// An array of all event directions
  mutable boost::scoped_array<DAEEventDirection> p_dirs;

  /// An array of all termination flags
  mutable boost::scoped_array<bool> p_terms;

  /// An array of all event values
  mutable boost::scoped_array<T> p_values;

  /// An array of trigger flags
  mutable boost::scoped_array<bool> p_trigger;

  /// A flag to indicate the solver has been terminated due to an event
  bool p_terminated;
};



} // namespace math
} // namespace gridpack


#endif
