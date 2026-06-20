/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_components.cpp
 * @author Bruce Palmer
 * @date   2023-11-29 14:57:06 d3g096
 * 
 * @updated Shri Abhyankar
 * Conversion of constant current, constant admittance values from raw file
 * to constant power
 * @date  2022-12-23
 *
 * @updated Yousu Chen
 * - Improved Q-limit handling with iterative PV-PQ conversion
 * - Added island detection function
 * - Automatic slack bus transfer for contingency analysis
 * - RMPCT-based reactive power distribution for multi-generator buses
 * @date  2026-01-31
 *
 * @updated Yousu Chen
 * - Added initStart option for power flow initialization (warm/flat start)
 * - Treated generators with Qmax == Qmin (zero Q capability) as PQ buses
 * - Q distribution uses RMPCT when available, otherwise uses relative reactive
 *   capability (Qmax) to share reactive power among generators
 * @date  2026-02-02
 *
 * @updated Yousu Chen
 * - Implemented voltage-dependent ZIP load model
 * @date  2026-02-19
 *
 * @updated Yousu Chen
 * - Added switched shunt control (MODSW=1 discrete, MODSW=2 continuous)
 * @date  2026-02-24
 *
 * @updated Yousu Chen
 * - Added LTC (load tap changer) control on PFBranch
 * @date  2026-03-28
 *
 * @updated Yousu Chen
 * - IREG PV bus swap for remote voltage regulation
 * - Star bus filtering for 3-winding transformer output
 * - LTC tap direction fix for to-bus controlled transformers
 * @date  2026-04-04
 *
 * @brief Methods used in power flow application
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include <iostream>
#include <cstring>
#include <stdio.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "pf_components.hpp"
#include "gridpack/parser/dictionary.hpp"

//#define LARGE_MATRIX

// Static member initialization
gridpack::powerflow::InitStartMode gridpack::powerflow::PFBus::p_initStartMode = INIT_START_WARM;
bool gridpack::powerflow::PFBus::p_qlim = true;
double gridpack::powerflow::PFBus::p_qlim_deadband = 0.1;
std::vector<std::string> gridpack::powerflow::PFBus::p_qlimWarnings;

/**
 * Set the initial start mode for power flow solver
 */
void gridpack::powerflow::PFBus::setInitStartMode(InitStartMode mode)
{
  p_initStartMode = mode;
}

/**
 * Set the qlim flag
 */
void gridpack::powerflow::PFBus::setQlim(bool qlim)
{
  p_qlim = qlim;
}

void gridpack::powerflow::PFBus::setQlimDeadband(double db)
{
  p_qlim_deadband = db;
}

/**
 * Clear accumulated Q limit warning messages
 */
void gridpack::powerflow::PFBus::clearQlimWarnings()
{
  p_qlimWarnings.clear();
}

/**
 * Get accumulated Q limit warning messages
 */
std::vector<std::string>& gridpack::powerflow::PFBus::getQlimWarnings()
{
  return p_qlimWarnings;
}

/**
 *  Simple constructor
 */
gridpack::powerflow::PFBus::PFBus(void)
{
  p_shunt_gs = 0.0;
  p_shunt_bs = 0.0;
  p_v = 0.0;
  p_a = 0.0;
  p_theta = 0.0;
  p_angle = 0.0;
  p_voltage = 0.0;
  p_vmin = 0.0;
  p_vmax = 0.0;
  /*p_pl = 0.0;
  p_ql = 0.0;
  p_ip = 0.0;
  p_iq = 0.0;
  p_yp = 0.0;
  p_yq = 0.0;*/
  p_sbase = 0.0;
  p_mode = YBus;
  setReferenceBus(false);
  p_ngen = 0;
  p_data = NULL;
  p_ignore = false;
  p_isStarBus = false;
  p_isIREG_PV = false;
  p_ireg_vs = 0.0;
  p_ireg_remote_bus = 0;
  p_ireg_remote_v_ptr = NULL;
  p_vMag_ptr = NULL;
  p_vAng_ptr = NULL;
  p_PV_ptr = NULL;
  // Switched shunt defaults
  p_hasSwitchedShunt = false;
  p_swshunt_modsw = 0;
  p_swshunt_vswhi = 1.0;
  p_swshunt_vswlo = 1.0;
  p_swshunt_swrem = 0;
  p_swshunt_binit = 0.0;
  p_swshunt_adjm = 0;
  p_swshunt_nblocks = 0;
  p_swshunt_bcurrent = 0.0;
  p_swshunt_bmin = 0.0;
  p_swshunt_bmax = 0.0;
  p_swshunt_locked = false;
  p_swshunt_adj_count = 0;
  p_swshunt_bprev = 0.0;
}

/**
 *  Simple destructor
 */
gridpack::powerflow::PFBus::~PFBus(void)
{
}

/**
 *  Return size of matrix block contributed by the component
 *  @param isize, jsize: number of rows and columns of matrix block
 *  @return: false if network component does not contribute matrix element
 */
bool gridpack::powerflow::PFBus::matrixDiagSize(int *isize, int *jsize) const
{
  if (p_mode == Jacobian) {
    if (!isIsolated()) {
#ifdef LARGE_MATRIX
      *isize = 2;
      *jsize = 2;
      return true;
#else
      if (getReferenceBus()) {
        return false;
      } else if (p_isPV && !p_isIREG_PV) {
        *isize = 1;
        *jsize = 1;
        return true;
      } else {
        *isize = 2;
        *jsize = 2;
        return true;
      }
#endif
    } else {
      return false;
    }
  } else if (p_mode == YBus) {
    return YMBus::matrixDiagSize(isize,jsize);
  }
  return true;
}

/**
 * Return the values of the matrix block. The values are
 * returned in row-major order.
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::powerflow::PFBus::matrixDiagValues(ComplexType *values)
{
  if (p_mode == YBus) {
    return YMBus::matrixDiagValues(values);
  } else if (p_mode == Jacobian) {
    double rvals[4];
    int nvals = diagonalJacobianValues(rvals);
    for (int i=0; i<nvals; i++) values[i] = rvals[i];
    if (nvals == 0) {
      return false;
    } else  {
      return true;
    }
  }
  return false;
}

bool gridpack::powerflow::PFBus::matrixDiagValues(RealType *values)
{
  if (p_mode == Jacobian) {
    int nvals = diagonalJacobianValues(values);
    if (nvals == 0) {
      return false;
    } else  {
      return true;
    }
  }
  return false;
}

/**
 * Return the size of the block that this component contributes to the
 * vector
 * @param size: size of vector block
 * @return: false if component does not contribute to vector
 */
bool gridpack::powerflow::PFBus::vectorSize(int *size) const
{
  if (p_mode == RHS || p_mode == State) {
    if (!isIsolated()) {
#ifdef LARGE_MATRIX
      *size = 2;
      return true;
#else
      if (getReferenceBus()) {
        return false;
      } else if (p_isPV && !p_isIREG_PV) {
        *size = 1;
      } else {
        *size = 2;
      }
      return true;
#endif
    } else {
      return false;
    }
  } else if (p_mode == S_Cal){
    *size = 1;
  } else {
    *size = 2;
  }
  return true;
}

/**
 * Return the values of the vector block
 * @param values: pointer to vector values
 * @return: false if network component does not contribute
 *        vector element
 */
bool gridpack::powerflow::PFBus::vectorValues(ComplexType *values)
{
  if (p_mode == S_Cal)  {
    double retr = p_v * cos(p_a);
    double reti = p_v * sin(p_a);
    gridpack::ComplexType ret(retr, reti);
    values[0] = ret;
    return true;
  } else if (p_mode == State) {
    values[0] = p_v;
    values[1] = p_a;
    return true;
  } else if (p_mode == RHS) {
    double rvals[2];
    int nvals = rhsValues(rvals);
    for (int i=0; i<nvals; i++) values[i] = rvals[i];
    if (nvals == 0) {
      return false;
    } else {
      return true;
    }
  }
  return false;
}

bool gridpack::powerflow::PFBus::vectorValues(RealType *values)
{
  if (p_mode == State) {
    values[0] = p_v;
    values[1] = p_a;
    return true;
  }
  if (p_mode == RHS) {
    int nvals = rhsValues(values);
    if (nvals == 0) {
      return false;
    } else {
      return true;
    }
  }
  return false;
}

/**
 * Check QLIM with RMPCT-based iterative Q distribution 
 *
 * Algorithm:
 * 1. Calculate required Q for the bus
 * 2. Distribute Q among generators using RMPCT ratios
 * 3. Check each generator against its limits
 * 4. If any generator exceeds limits, clamp it and redistribute to others
 * 5. Only convert PV->PQ if all generators are at limits and Q_req still not met
 *
 * @return false: violations exist (Q limits hit, bus converted to PQ)
 * @return true:  no violations (all generators within limits)
 */
bool gridpack::powerflow::PFBus::chkQlim(double q_deadband)
{
  if (!p_isPV) {
    if (p_PV_ptr) *p_PV_ptr = p_isPV;
    return false;
  }

  int ngen = p_gstatus.size();
  if (ngen == 0) {
    if (p_PV_ptr) *p_PV_ptr = p_isPV;
    return false;
  }

  // Calculate total load Q (constant-power part)
  double ql = 0.0;
  for (int i = 0; i < p_lstatus.size(); i++) {
    if (p_lstatus[i] == 1) {
      ql += p_ql[i];
    }
  }
  // Add voltage-dependent Q load (IQ*V - YQ*V^2)
  double pzip, qzip;
  getZIPLoadPower(p_v, pzip, qzip);
  ql += qzip;

  // Required Q from generators = Q_injection + Q_load
  double Q_required = p_Qinj * p_sbase + ql;

  // Track which generators are still active (not at limits)
  std::vector<bool> at_limit(ngen, false);
  std::vector<double> q_assigned(ngen, 0.0);

  // Save original state for restoration
  p_save2isPV = p_isPV;

  // Iterative distribution loop
  const int MAX_ITER = 20;  // Prevent infinite loops
  bool converged = false;

  for (int iter = 0; iter < MAX_ITER && !converged; iter++) {
    // Calculate total RMPCT and Qmax for active (non-limited) generators
    double total_rmpct = 0.0;
    double total_qmax = 0.0;
    int active_count = 0;

    for (int i = 0; i < ngen; i++) {
      if (p_gstatus[i] == 1 && !at_limit[i]) {
        total_rmpct += p_rmpct[i];
        total_qmax += p_qmax[i];
        active_count++;
      }
    }

    // If no active generators left, all are at limits
    if (active_count == 0) {
      break;
    }

    // Calculate remaining Q to distribute (subtract already-assigned Q from limited gens)
    double Q_remaining = Q_required;
    for (int i = 0; i < ngen; i++) {
      if (p_gstatus[i] == 1 && at_limit[i]) {
        Q_remaining -= q_assigned[i];
      }
    }

    // Choose distribution basis:
    // 1. RMPCT if available (total_rmpct > 0)
    // 2. Otherwise, use relative reactive capability (Qmax)
    bool use_rmpct = (total_rmpct > 0.0);
    double total_basis;
    if (use_rmpct) {
      total_basis = total_rmpct;
    } else {
      total_basis = total_qmax;
    }

    if (total_basis <= 0.0) {
      // No basis for distribution - all generators have zero capability
      break;
    }

    // Distribute Q_remaining among active generators
    converged = true;
    for (int i = 0; i < ngen; i++) {
      if (p_gstatus[i] != 1 || at_limit[i]) continue;

      double basis;
      if (use_rmpct) {
        basis = p_rmpct[i];
      } else {
        basis = p_qmax[i];
      }
      double share = (basis / total_basis) * Q_remaining;

      // Check against limits
      if (share > p_qmax[i]) {
        q_assigned[i] = p_qmax[i];
        at_limit[i] = true;
        converged = false;  // Need another iteration
      } else if (share < p_qmin[i]) {
        q_assigned[i] = p_qmin[i];
        at_limit[i] = true;
        converged = false;  // Need another iteration
      } else {
        q_assigned[i] = share;
      }
    }
  }

  // Calculate total Q that can be supplied
  double Q_supplied = 0.0;
  double Q_max_total = 0.0;
  double Q_min_total = 0.0;
  for (int i = 0; i < ngen; i++) {
    if (p_gstatus[i] == 1) {
      Q_supplied += q_assigned[i];
      Q_max_total += p_qmax[i];
      Q_min_total += p_qmin[i];
    }
  }

  // Check if Q requirement can be met.
  // q_deadband avoids switching buses that are only marginally over their Q
  // limit due to floating-point differences. Configurable via XML qlimDeadband
  // (default 0.1 Mvar)
  bool need_pv_to_pq = false;
  char warnBuf[256];
  if (Q_required > Q_max_total + q_deadband) {
    // Exceeds total Qmax - need to convert to PQ
    snprintf(warnBuf, sizeof(warnBuf),
             "\nWarning: Bus %d Q requirement (%8.3f) exceeds total QMAX (%8.3f), converting to PQ\n",
             getOriginalIndex(), Q_required, Q_max_total);
    printf("%s", warnBuf);  // Keep screen output
    p_qlimWarnings.push_back(std::string(warnBuf));  // Store for file output
    need_pv_to_pq = true;
    // Clamp all generators to Qmax
    for (int i = 0; i < ngen; i++) {
      if (p_gstatus[i] == 1) {
        q_assigned[i] = p_qmax[i];
      }
    }
  } else if (Q_required < Q_min_total - q_deadband) {
    // Below total Qmin - need to convert to PQ
    snprintf(warnBuf, sizeof(warnBuf),
             "\nWarning: Bus %d Q requirement (%8.3f) below total QMIN (%8.3f), converting to PQ\n",
             getOriginalIndex(), Q_required, Q_min_total);
    printf("%s", warnBuf);  // Keep screen output
    p_qlimWarnings.push_back(std::string(warnBuf));  // Store for file output
    need_pv_to_pq = true;
    // Clamp all generators to Qmin
    for (int i = 0; i < ngen; i++) {
      if (p_gstatus[i] == 1) {
        q_assigned[i] = p_qmin[i];
      }
    }
  }

  // Apply the Q assignments to generators
  for (int i = 0; i < ngen; i++) {
    p_gstatus_save.push_back(p_gstatus[i]);
    if (p_gstatus[i] == 1) {
      p_qg[i] = q_assigned[i];
    }
  }

  // Convert PV to PQ if needed
  if (need_pv_to_pq) {
    p_isPV = false;
    p_type = 1;  // Change bus type from PV(2) to PQ(1)
    if (p_PV_ptr) *p_PV_ptr = false;
    return true;  // Violation found
  }

  // Bus stays PV - Q was successfully distributed within limits
  if (p_PV_ptr) *p_PV_ptr = p_isPV;
  return false;  // No violation
}

/**
 * Clear changes that were made for Q limit violations and reset
 * bus to its original state
 *
 * Note: chkQlim() does NOT modify p_gstatus - generators stay online when Q limits
 * are hit. It only converts PV buses to PQ and clamps Q to limits. Therefore,
 * clearQlim() should NOT restore p_gstatus. Generator status restoration is handled
 * by unSetContingency() which uses setGenStatus().
 *
 * The p_gstatus_save mechanism was designed for a different Q limit implementation
 * that turned off generators. In the current implementation, we only need to:
 * 1. Clear p_gstatus_save to reset the tracking state
 * 2. Restore p_isPV if there are online generators
 */
void gridpack::powerflow::PFBus::clearQlim()
{
  // Clear p_gstatus_save to reset Q limit violation tracking for next iteration.
  // Do NOT restore p_gstatus from p_gstatus_save - that would undo the generator
  // restoration done by unSetContingency().
  p_gstatus_save.clear();

  // Only restore p_isPV to true (PV bus) if there's at least one online generator.
  // This prevents an inconsistent state when a generator contingency tripped a generator
  // and caused a Q limit violation due to qmax=0. Without this check, the bus would be
  // restored to PV status but with no online generator, causing subsequent power flows
  // to immediately hit Q limit violations and potentially diverge.
  bool hasOnlineGen = false;
  int ngen = p_gstatus.size();
  for (int i = 0; i < ngen; i++) {
    if (p_gstatus[i] == 1) {
      hasOnlineGen = true;
      break;
    }
  }

  if (hasOnlineGen) {
    p_isPV = p_save2isPV;
    p_type = p_save_type;  // Restore original bus type (e.g., from PQ back to PV)
  }
  // If no online generators, keep p_isPV and p_type unchanged

  // Always restore generator Q to original values.
  // During Q limit handling, p_qg was clamped to limits. This caused setSBus()
  // to compute a different scheduled Q (p_Q0) than the original base case.
  // Restoring p_qg ensures contingencies start with the same scheduled power.
  for (int i = 0; i < ngen; i++) {
    p_qg[i] = p_saveQg[i];
  }

  if (p_PV_ptr) *p_PV_ptr = p_isPV;
}

/**
 * Adjust switched shunt B to bring controlled voltage within deadband.
 * @param v_controlled voltage at controlled bus (local or remote)
 * @return true if an adjustment was made
 */
bool gridpack::powerflow::PFBus::adjustSwitchedShunt(double v_controlled)
{
  if (!p_hasSwitchedShunt || p_swshunt_locked) return false;

  // Skip adjustment for generator buses
  // Shunt control only becomes active after PV->PQ conversion
  if (p_isPV || getReferenceBus()) return false;

  // Check if voltage is within deadband.
  // When VSWHI == VSWLO (zero-width deadband), add a minimum tolerance to
  // prevent cycling due to floating-point representation.
  double vswlo = p_swshunt_vswlo;
  double vswhi = p_swshunt_vswhi;
  if (vswhi - vswlo < 0.002) {
    double vmid = (vswhi + vswlo) / 2.0;
    vswlo = vmid - 0.001;
    vswhi = vmid + 0.001;
  }
  if (v_controlled >= vswlo && v_controlled <= vswhi) {
    return false;  // No action needed
  }

  double delta_b = 0.0;

  if (p_swshunt_modsw == 1) {
    // MODSW=1: Discrete adjustment - one step per controller iteration
    if (v_controlled < p_swshunt_vswlo) {
      // Voltage too low - increase B (add capacitor step or remove reactor step)
      bool adjusted = false;
      // Try adding a capacitor step (positive B blocks, left to right)
      for (int ib = 0; ib < p_swshunt_nblocks; ib++) {
        if (p_swshunt_b[ib] > 0.0 && p_swshunt_steps_on[ib] < p_swshunt_n[ib]) {
          p_swshunt_steps_on[ib]++;
          delta_b = p_swshunt_b[ib];
          adjusted = true;
          break;
        }
      }
      // If no capacitor step available, try removing a reactor step (negative B blocks)
      if (!adjusted) {
        for (int ib = p_swshunt_nblocks - 1; ib >= 0; ib--) {
          if (p_swshunt_b[ib] < 0.0 && p_swshunt_steps_on[ib] > 0) {
            p_swshunt_steps_on[ib]--;
            delta_b = -p_swshunt_b[ib];  // Removing reactor increases B
            adjusted = true;
            break;
          }
        }
      }
      if (!adjusted) return false;
    } else {
      // Voltage too high - decrease B (remove capacitor step or add reactor step)
      bool adjusted = false;
      // Try removing a capacitor step (positive B blocks, right to left)
      for (int ib = p_swshunt_nblocks - 1; ib >= 0; ib--) {
        if (p_swshunt_b[ib] > 0.0 && p_swshunt_steps_on[ib] > 0) {
          p_swshunt_steps_on[ib]--;
          delta_b = -p_swshunt_b[ib];
          adjusted = true;
          break;
        }
      }
      // If no capacitor step to remove, try adding a reactor step (negative B blocks)
      if (!adjusted) {
        for (int ib = 0; ib < p_swshunt_nblocks; ib++) {
          if (p_swshunt_b[ib] < 0.0 && p_swshunt_steps_on[ib] < p_swshunt_n[ib]) {
            p_swshunt_steps_on[ib]++;
            delta_b = p_swshunt_b[ib];  // Adding reactor decreases B
            adjusted = true;
            break;
          }
        }
      }
      if (!adjusted) return false;
    }
  } else if (p_swshunt_modsw == 2) {
    // MODSW=2: Continuous adjustment toward deadband midpoint
    double v_target = (p_swshunt_vswhi + p_swshunt_vswlo) / 2.0;
    double k = p_swshunt_bmax - p_swshunt_bmin;
    if (k <= 0.0) k = 1.0;  // Fallback
    double b_target = p_swshunt_bcurrent + k * (v_target - v_controlled);
    // Clamp to [Bmin, Bmax]
    if (b_target < p_swshunt_bmin) b_target = p_swshunt_bmin;
    if (b_target > p_swshunt_bmax) b_target = p_swshunt_bmax;
    delta_b = b_target - p_swshunt_bcurrent;
    if (fabs(delta_b) < 1.0e-6) return false;  // No meaningful change
  } else {
    return false;
  }

  // Cycle detection: if this adjustment would return to previous B, we are
  // oscillating between two discrete steps that straddle the deadband.
  // Choose the better of the two states (less deviated from the deadband)
  // before locking.  If the current voltage is noticeably outside the deadband,
  // apply the reverting step so we lock at the less-deviated setting.
  double b_new = p_swshunt_bcurrent + delta_b;
  if (p_swshunt_adj_count > 0 && fabs(b_new - p_swshunt_bprev) < 1.0e-6) {
    double outside = (v_controlled > vswhi) ? (v_controlled - vswhi)
                                             : (vswlo - v_controlled);
    if (outside > 0.01) {
      // Current state is far outside deadband; accept the reverting step
      // (which moves toward the deadband) and then lock.
      p_swshunt_bprev = p_swshunt_bcurrent;
      p_swshunt_bcurrent = b_new;
      addShuntValues(0.0, delta_b / p_sbase);
      p_swshunt_locked = true;
      return true;
    }
    p_swshunt_locked = true;
    return false;
  }

  // Apply the change to YMBus::p_shunt_bs via addShuntValues()
  p_swshunt_bprev = p_swshunt_bcurrent;
  p_swshunt_bcurrent = b_new;
  addShuntValues(0.0, delta_b / p_sbase);

  // Increment adjustment count and lock out after 10 adjustments
  p_swshunt_adj_count++;
  if (p_swshunt_adj_count >= 10) {
    p_swshunt_locked = true;
  }

  return true;
}

/**
 * Reset switched shunt to BINIT state
 */
void gridpack::powerflow::PFBus::resetSwitchedShunt()
{
  if (!p_hasSwitchedShunt) return;

  // Remove current B from shunt, restore to BINIT
  double delta = p_swshunt_binit - p_swshunt_bcurrent;
  addShuntValues(0.0, delta / p_sbase);
  p_swshunt_bcurrent = p_swshunt_binit;

  // Re-decompose BINIT into step counts
  double b_remaining = p_swshunt_binit;
  for (int ib = 0; ib < p_swshunt_nblocks; ib++) {
    if (p_swshunt_b[ib] == 0.0) {
      p_swshunt_steps_on[ib] = 0;
      continue;
    }
    int steps = static_cast<int>(b_remaining / p_swshunt_b[ib]);
    if (steps < 0) steps = 0;
    if (steps > p_swshunt_n[ib]) steps = p_swshunt_n[ib];
    p_swshunt_steps_on[ib] = steps;
    b_remaining -= steps * p_swshunt_b[ib];
  }

  p_swshunt_locked = false;
  p_swshunt_adj_count = 0;
  p_swshunt_bprev = p_swshunt_binit;
}

/**
 * Check if bus has a switched shunt
 */
bool gridpack::powerflow::PFBus::hasSwitchedShunt() const
{
  return p_hasSwitchedShunt;
}

/**
 * Get remote bus number for switched shunt control
 */
int gridpack::powerflow::PFBus::getSwitchedShuntRemoteBus() const
{
  return p_swshunt_swrem;
}

/**
 * Get current switched shunt susceptance
 */
double gridpack::powerflow::PFBus::getSwitchedShuntB() const
{
  return p_swshunt_bcurrent;
}


/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto buses
 * @param values array containing voltage magnitude and angle
 */
void gridpack::powerflow::PFBus::setValues(gridpack::ComplexType *values)
{
  double vt = p_v;
  double at = p_a;
  p_a -= real(values[0]);
#ifdef LARGE_MATRIX
  p_v -= real(values[1]);
#else
  if (!p_isPV || p_isIREG_PV) {
    p_v -= real(values[1]);
  }
#endif
  *p_vMag_ptr = p_v;
  double pi = 4.0*atan(1.0);
  if (p_a >= 0.0) {
    *p_vAng_ptr = fmod(p_a+pi,2.0*pi)-pi;
  } else {
    *p_vAng_ptr = fmod(p_a-pi,2.0*pi)+pi;
  }
}

void gridpack::powerflow::PFBus::setValues(gridpack::RealType *values)
{
  double vt = p_v;
  double at = p_a;
  p_a -= values[0];
#ifdef LARGE_MATRIX
  //p_v -= real(values[1]);
  p_v -= values[1];
#else
  if (!p_isPV) {
    p_v -= values[1];
  }
#endif
  *p_vMag_ptr = p_v;
  double pi = 4.0*atan(1.0);
  if (p_a >= 0.0) {
    *p_vAng_ptr = fmod(p_a+pi,2.0*pi)-pi;
  } else {
    *p_vAng_ptr = fmod(p_a-pi,2.0*pi)+pi;
  }
}

/**
 * Return the size of the buffer used in data exchanges on the network.
 * For this problem, the voltage magnitude and phase angle need to be exchanged
 * @return size of buffer
 */
int gridpack::powerflow::PFBus::getXCBufSize(void)
{
  return (2*sizeof(double)+sizeof(bool));
}

/**
 * Assign pointers for voltage magnitude and phase angle
 */
void gridpack::powerflow::PFBus::setXCBuf(void *buf)
{
  p_vAng_ptr = static_cast<double*>(buf);
  p_vMag_ptr = p_vAng_ptr+1;
  void *ptr = static_cast<void*>(p_vMag_ptr+1);
  p_PV_ptr = static_cast<bool*>(ptr);
  // Note: we are assuming that the load function has been called BEFORE
  // the factory setExchange method, so p_a and p_v are set with their initial
  // values.
  *p_vMag_ptr = p_v;
  double pi = 4.0*atan(1.0);
  if (p_a >= 0.0) {
    *p_vAng_ptr = fmod(p_a+pi,2.0*pi)-pi;
  } else {
    *p_vAng_ptr = fmod(p_a-pi,2.0*pi)+pi;
  }
  *p_PV_ptr = p_isPV;
  
}

/**
 * Load values stored in DataCollection object into PFBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::powerflow::PFBus::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  p_data = data.get();
  YMBus::load(data);

  // This routine may be called more than once, so clear all vectors
  p_pg.clear();
  p_qg.clear();
  p_pFac.clear();
  p_pFac_orig.clear();
  p_qFac.clear();
  p_qFac_orig.clear();
  p_gstatus.clear();
  p_gstatus_save.clear();
  p_qmin.clear();
  p_qmax.clear();
  p_qmin_orig.clear();
  p_qmax_orig.clear();
  p_rmpct.clear();
  p_gid.clear();
  p_pt.clear();
  p_pb.clear();
  p_pl.clear();
  p_ql.clear();
  p_ip.clear();
  p_iq.clear();
  p_yp.clear();
  p_yq.clear();
  p_saveIp.clear();
  p_saveIq.clear();
  p_saveYp.clear();
  p_saveYq.clear();
  p_lstatus.clear();
  p_lid.clear();

  bool ok = data->getValue(CASE_SBASE, &p_sbase);
  data->getValue(BUS_VOLTAGE_ANG, &p_angle);
  data->getValue(BUS_VOLTAGE_MAG, &p_voltage);
  // Load bus voltage limits (NVHI/NVLO from RAW) for IREG clamping
  p_vmax = 1.1;  // default
  p_vmin = 0.9;  // default
  data->getValue(BUS_VOLTAGE_MAX, &p_vmax);
  data->getValue(BUS_VOLTAGE_MIN, &p_vmin);
  double pi = 4.0*atan(1.0);

  // Apply initial start mode
  if (p_initStartMode == INIT_START_FLAT) {
    // Flat start: all angles to 0, PQ buses to 1.0 pu
    // (PV/Slack voltage will be set to VS later when processing generators)
    p_v = 1.0;       // Default to 1.0 pu for flat start (will be overridden for PV/Slack)
    p_voltage = 1.0; // Also set p_voltage for resetVoltage()
    p_angle = 0.0;   // Flat angle
    p_a = 0.0;
  } else {
    // Warm start: use values from raw file
    p_v = p_voltage;
    p_angle = p_angle*pi/180.0;
    p_a = p_angle;
  }
  data->getValue(BUS_TYPE, &p_type);
  p_save_type = p_type;  // Save original bus type for restoration after Q limit handling
  if (p_type == 3) {
    setReferenceBus(true);
  }
  if (isIsolated()) {
    p_original_isolated = true;
  } else {
    p_original_isolated = false;
  }
  // Detect synthetic star buses from 3-winding transformers
  std::string busName;
  p_isStarBus = false;
  if (data->getValue(BUS_NAME, &busName)) {
    if (busName.substr(0, 9) == "DUMMY_BUS") {
      p_isStarBus = true;
    }
  }
  p_area = 1;
  data->getValue(BUS_AREA, &p_area);
  p_zone = 1;
  data->getValue(BUS_ZONE, &p_zone);

  // if BUS_TYPE = 2, and gstatus is 1, then bus is a PV bus
  p_isPV = false;

  // added p_pg,p_qg,p_pl,p_ql,p_sbase;

  bool lgen;
  int i, gstatus;
  double pg, qg, vs,qmax,qmin;
  int ngen = 0;
  p_ngen = 0;
  if (data->getValue(GENERATOR_NUMBER, &ngen)) {
    double pcaptot = 0.0;
    double qcaptot = 0.0;   // Total Qmax for Q factor calculation
    for (i=0; i<ngen; i++) {
      lgen = true;
      lgen = lgen && data->getValue(GENERATOR_PG, &pg,i);
      lgen = lgen && data->getValue(GENERATOR_QG, &qg,i);
      lgen = lgen && data->getValue(GENERATOR_VS, &vs,i);
      lgen = lgen && data->getValue(GENERATOR_STAT, &gstatus,i);
      lgen = lgen && data->getValue(GENERATOR_QMAX, &qmax,i);
      lgen = lgen && data->getValue(GENERATOR_QMIN, &qmin,i);
      double pt = 0.0;
      double pb = 0.0;
      double rmpct = 100.0;  // Default RMPCT value per PSS/E
      ok =  data->getValue(GENERATOR_PMAX,&pt,i);
      ok =  data->getValue(GENERATOR_PMIN,&pb,i);
      // Try to get RMPCT, use default 100.0 if not available
      data->getValue(GENERATOR_RMPCT, &rmpct, i);
      if (rmpct < 0.0) rmpct = 0.0;  // Treat negative as zero
      if (lgen) {
        p_gstatus.push_back(gstatus);
        // Store original Q limits from input file BEFORE zeroing for offline generators
        p_qmax_orig.push_back(qmax);
        p_qmin_orig.push_back(qmin);
        if (gstatus == 0) {
          qmax = 0.0;
          qmin = 0.0;
	  pg = 0.0;
	  qg = 0.0;
        }
	p_pg.push_back(pg);
        p_savePg.push_back(pg);
        p_qg.push_back(qg);
        p_saveQg.push_back(qg);
        p_qmax.push_back(qmax);
        p_qmin.push_back(qmin);
        p_pt.push_back(pt);
        p_pb.push_back(pb);
        p_vs.push_back(vs);
        int ireg = 0;
        data->getValue(GENERATOR_IREG, &ireg, i);
        p_ireg.push_back(ireg);
        p_rmpct.push_back(gstatus ? rmpct : 0.0);  // Zero RMPCT for offline generators
	if(gstatus) {
	  p_pFac.push_back(pt - pb);
	  pcaptot += pt - pb;
	  // For p_qFac, use Qmax for generators with reactive capability (Qmax > Qmin)
	  // Generators with Qmax == Qmin (including Qmax=Qmin=0) get zero p_qFac
	  if (qmax > qmin) {
	    p_qFac.push_back(qmax);
	    qcaptot += qmax;
	  } else {
	    p_qFac.push_back(0.0);
	  }
	} else {
	  p_pFac.push_back(0.0);
	  p_qFac.push_back(0.0);
	}
	p_pFac_orig.push_back(pt - pb);
	p_qFac_orig.push_back(qmax > qmin ? qmax : 0.0);

        if (gstatus == 1) {
          // Use VS for PV/Slack buses when:
          // - flat start, OR
          // - warm start with qlim=false
          if (p_initStartMode == INIT_START_FLAT || !p_qlim) {
            p_v = vs;        // Set voltage to generator setpoint
            p_voltage = vs;  // Also update p_voltage so resetVoltage() uses VS
          }
          // Only set PV if generator has reactive power capability (Qmax != Qmin)
          // Generators with Qmax == Qmin cannot regulate voltage and should be
          // treated as PQ injections with fixed Q output
          if (p_type == 2 && qmax != qmin) p_isPV = true;
        }
        std::string id("-1");
        data->getValue(GENERATOR_ID,&id,i);
        p_gid.push_back(id);
        p_ngen++;
      }
    }
    // Detect IREG PV buses: if ALL online generators have remote IREG,
    // use augmented Jacobian (control V at remote bus instead of local bus)
    if (p_isPV && p_ngen > 0) {
      bool all_remote = true;
      int ireg_bus_found = 0;
      double ireg_vs_found = 0.0;
      int orig_idx = getOriginalIndex();
      for (int ig = 0; ig < p_ngen; ig++) {
        if (p_gstatus[ig] == 1) {
          int ireg = p_ireg[ig];
          if (ireg == 0 || ireg == orig_idx) {
            all_remote = false;
            break;
          }
          if (ireg_bus_found == 0) {
            ireg_bus_found = ireg;
            ireg_vs_found = p_vs[ig];
          }
        }
      }
      if (all_remote && ireg_bus_found != 0) {
        // Store IREG info but don't enable augmented Jacobian yet
        // The factory will handle PV swap after all buses are loaded
        p_ireg_remote_bus = ireg_bus_found;
        p_ireg_vs = ireg_vs_found;
        // p_isIREG_PV is set by factory after PV swap setup
      }
    }

    // Normalize p_pFac (P-based factor for P distribution)
    if (pcaptot != 0.0 && p_ngen > 1) {
      for (i=0; i<p_ngen; i++) {
        p_pFac[i] = p_pFac[i]/pcaptot;
        p_pFac_orig[i] = p_pFac[i];
      }
    } else {
      p_pFac[0] = 1.0;
      p_pFac_orig[0] = p_pFac[0];
    }
    // Normalize p_qFac (Q-based factor for Q distribution when qlim=false)
    // Uses Qmax to distribute Q among generators with reactive capability
    // Generators with Qmax == Qmin (zero Q capability) get zero p_qFac
    if (qcaptot != 0.0 && p_ngen > 1) {
      for (i=0; i<p_ngen; i++) {
        p_qFac[i] = p_qFac[i]/qcaptot;
        p_qFac_orig[i] = p_qFac[i];
      }
    } else if (p_ngen > 0) {
      // If no Q capability (qcaptot=0) or single generator, first gen with
      // any capability gets all Q, otherwise distribute equally
      bool found = false;
      for (i=0; i<p_ngen; i++) {
        if (p_qFac_orig[i] > 0.0 && !found) {
          p_qFac[i] = 1.0;
          p_qFac_orig[i] = 1.0;
          found = true;
        } else if (!found && i == p_ngen - 1) {
          // No generators with Q capability - give to first online gen
          // (will be clamped to zero anyway if Qmax=Qmin=0)
          for (int j=0; j<p_ngen; j++) {
            if (p_gstatus[j] == 1) {
              p_qFac[j] = 1.0;
              p_qFac_orig[j] = 1.0;
              break;
            }
          }
        }
      }
    }
  }
  // Warm-start Q-limit pre-saturation: scheduled QG at QMAX/QMIN -> start as PQ.
  if (p_isPV && p_qlim && p_initStartMode != INIT_START_FLAT) {
    double total_qg = 0.0, total_qmax = 0.0, total_qmin = 0.0;
    for (i = 0; i < p_ngen; i++) {
      if (p_gstatus[i] == 1) {
        total_qg   += p_qg[i];
        total_qmax += p_qmax[i];
        total_qmin += p_qmin[i];
      }
    }
    if (total_qg >= total_qmax - p_qlim_deadband
        || total_qg <= total_qmin + p_qlim_deadband) {
      p_isPV = false;
      // Pick the saturated limit from V vs VS (V>VS => QMIN, V<VS => QMAX), not
      // scheduled QG which can be inconsistent. Local regulation only.
      if (p_ireg_remote_bus == 0) {
        double vset = 0.0; int nv = 0;
        for (i = 0; i < p_ngen; i++)
          if (p_gstatus[i] == 1) { vset += p_vs[i]; nv++; }
        if (nv > 0) vset /= nv;
        const double V_init_tol = 1.0e-4;
        if (p_voltage > vset + V_init_tol) {
          for (i = 0; i < p_ngen; i++) if (p_gstatus[i] == 1) p_qg[i] = p_qmin[i];
        } else if (p_voltage < vset - V_init_tol) {
          for (i = 0; i < p_ngen; i++) if (p_gstatus[i] == 1) p_qg[i] = p_qmax[i];
        }
      }
    }
  }

  p_saveisPV = p_isPV;
  p_save2isPV = p_isPV;  // Initialize for Q limit handling - ensures valid value if chkQlim() never called

// Add load
  int lstatus;
  double pl,ql,ip,iq,yp,yq;
  p_load = true;
  p_load = p_load && data->getValue(LOAD_PL, &pl,0);
  p_load = p_load && data->getValue(LOAD_QL, &ql,0);
  p_load = p_load && data->getValue(LOAD_IP, &ip,0);
  p_load = p_load && data->getValue(LOAD_IQ, &iq,0);
  p_load = p_load && data->getValue(LOAD_YP, &yp,0);
  p_load = p_load && data->getValue(LOAD_YQ, &yq,0);
  int nld = 0;
  p_nload = 0;
  if (data->getValue(LOAD_NUMBER, &nld)) {
    for (i=0; i<nld; i++) {
      p_load = true;
      p_load = p_load && data->getValue(LOAD_PL, &pl,i);
      p_load = p_load && data->getValue(LOAD_QL, &ql,i);
      p_load = p_load && data->getValue(LOAD_IP, &ip,i);
      p_load = p_load && data->getValue(LOAD_IQ, &iq,i);
      p_load = p_load && data->getValue(LOAD_YP, &yp,i);
      p_load = p_load && data->getValue(LOAD_YQ, &yq,i);
      p_load = p_load && data->getValue(LOAD_STATUS, &lstatus,i);
      if (p_load) {
	/* Store constant-power load only; ZIP components stored separately */
        p_pl.push_back(pl);
        p_savePl.push_back(pl);
        p_ql.push_back(ql);
        p_saveQl.push_back(ql);
        p_ip.push_back(ip);
        p_iq.push_back(iq);
        p_yp.push_back(yp);
        p_yq.push_back(yq);
        p_saveIp.push_back(ip);
        p_saveIq.push_back(iq);
        p_saveYp.push_back(yp);
        p_saveYq.push_back(yq);
        p_lstatus.push_back(lstatus);
        std::string id("-1");
        data->getValue(LOAD_ID,&id,i);
        p_lid.push_back(id);
        p_nload++;
      }
    }
  }
  // Load switched shunt data
  p_hasSwitchedShunt = false;
  p_swshunt_n.clear();
  p_swshunt_b.clear();
  p_swshunt_steps_on.clear();
  p_swshunt_nblocks = 0;
  p_swshunt_locked = false;
  p_swshunt_adj_count = 0;

  int swshunt_modsw = 0;
  int swshunt_stat = 1;  // Default in-service (v23 format has no STATUS field)
  if (data->getValue(SHUNT_MODSW, &swshunt_modsw)) {
    data->getValue(SHUNT_SWCH_STAT, &swshunt_stat);  // Override if present (v33+)
    if (swshunt_stat == 1 && (swshunt_modsw == 1 || swshunt_modsw == 2)) {
      p_swshunt_modsw = swshunt_modsw;
      p_swshunt_vswhi = 1.0;
      p_swshunt_vswlo = 1.0;
      p_swshunt_swrem = 0;
      p_swshunt_binit = 0.0;
      p_swshunt_adjm = 0;
      data->getValue(SHUNT_VSWHI, &p_swshunt_vswhi);
      data->getValue(SHUNT_VSWLO, &p_swshunt_vswlo);
      if (!data->getValue(SHUNT_SWREM, &p_swshunt_swrem)) {
        data->getValue(SHUNT_SWREG, &p_swshunt_swrem);  // v34+ uses SWREG
      }
      data->getValue(SHUNT_BINIT, &p_swshunt_binit);
      data->getValue(SHUNT_ADJM, &p_swshunt_adjm);

      // Read N1-N8, B1-B8 block data
      const char *n_keys[] = {SHUNT_N1, SHUNT_N2, SHUNT_N3, SHUNT_N4,
                               SHUNT_N5, SHUNT_N6, SHUNT_N7, SHUNT_N8};
      const char *b_keys[] = {SHUNT_B1, SHUNT_B2, SHUNT_B3, SHUNT_B4,
                               SHUNT_B5, SHUNT_B6, SHUNT_B7, SHUNT_B8};
      for (int ib = 0; ib < 8; ib++) {
        int ni = 0;
        double bi = 0.0;
        data->getValue(n_keys[ib], &ni);
        data->getValue(b_keys[ib], &bi);
        if (ni == 0 && bi == 0.0) break;  // End of blocks
        p_swshunt_n.push_back(ni);
        p_swshunt_b.push_back(bi);
        p_swshunt_nblocks++;
      }

      if (p_swshunt_nblocks > 0) {
        // Compute Bmin and Bmax from block data
        p_swshunt_bmin = 0.0;
        p_swshunt_bmax = 0.0;
        for (int ib = 0; ib < p_swshunt_nblocks; ib++) {
          double block_b = p_swshunt_n[ib] * p_swshunt_b[ib];
          if (block_b > 0.0) {
            p_swshunt_bmax += block_b;  // Capacitor block
          } else {
            p_swshunt_bmin += block_b;  // Reactor block
          }
        }

        // Decompose BINIT into initial step counts per block
        p_swshunt_steps_on.resize(p_swshunt_nblocks, 0);
        double b_remaining = p_swshunt_binit;
        for (int ib = 0; ib < p_swshunt_nblocks; ib++) {
          if (p_swshunt_b[ib] == 0.0) continue;
          int steps = static_cast<int>(b_remaining / p_swshunt_b[ib]);
          if (steps < 0) steps = 0;
          if (steps > p_swshunt_n[ib]) steps = p_swshunt_n[ib];
          p_swshunt_steps_on[ib] = steps;
          b_remaining -= steps * p_swshunt_b[ib];
        }

        p_swshunt_bcurrent = p_swshunt_binit;
        p_swshunt_bprev = p_swshunt_binit;
        p_hasSwitchedShunt = true;
      }
    }
  }

  // If this is being called a second time, then update pointers
  if (p_vMag_ptr) *p_vMag_ptr = p_v;
  if (p_vAng_ptr) {
    double pi = 4.0*atan(1.0);
    if (p_a >= 0.0) {
      *p_vAng_ptr = fmod(p_a+pi,2.0*pi)-pi;
    } else {
      *p_vAng_ptr = fmod(p_a-pi,2.0*pi)+pi;
    }
  }
}

/**
 * Set values of YBus matrix. These can then be used in subsequent
 * calculations
 */
void gridpack::powerflow::PFBus::setYBus(void)
{
  YMBus::setYBus();
  gridpack::ComplexType ret;
  ret = YMBus::getYBus();
  p_ybusr = real(ret);
  p_ybusi = imag(ret);
}

/**
 * Get values of YBus matrix. These can then be used in subsequent
 * calculations
 */
gridpack::ComplexType gridpack::powerflow::PFBus::getYBus(void)
{
  return YMBus::getYBus();
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::powerflow::PFBus::setMode(int mode)
{
  if (mode == YBus) {
    YMBus::setMode(gridpack::ymatrix::YBus);
  }
  p_mode = mode;
}

/**
 * Reset voltage and phase angle to initial values
 */
void gridpack::powerflow::PFBus::resetVoltage(void)
{
  p_v = p_voltage;
  p_a = p_angle;
  if (p_vMag_ptr) *p_vMag_ptr = p_v;
  if (p_vAng_ptr) {
    double pi = 4.0*atan(1.0);
    if (p_a >= 0.0) {
      *p_vAng_ptr = fmod(p_a+pi,2.0*pi)-pi;
    } else {
      *p_vAng_ptr = fmod(p_a-pi,2.0*pi)+pi;
    }
  }
}

/**
 * Set voltage limits on bus
 * @param vmin lower value of voltage
 * @param vmax upper value of voltage
 */
void gridpack::powerflow::PFBus::setVoltageLimits(double vmin, double vmax)
{
  p_vmin = vmin;
  p_vmax = vmax;
}

/**
 * Check voltage for violation
 * @return false if there is a voltage violation
 */
bool gridpack::powerflow::PFBus::checkVoltageViolation()
{
  bool ret = true;
  if (*p_vMag_ptr > p_vmax || *p_vMag_ptr < p_vmin) ret = false;
  return ret;
}

/**
 * Return the value of the voltage magnitude on this bus
 * @return: voltage magnitude
 */
double gridpack::powerflow::PFBus::getVoltage()
{
  return *p_vMag_ptr;
}

/**
 * Return whether or not the bus is a PV bus (V held fixed in powerflow
 * equations)
 * @return true if bus is PV bus
 */
bool gridpack::powerflow::PFBus::isPV(void)
{
  return p_isPV;
}

/**
 * Return the value of the phase angle on this bus
 * @return: phase angle
 */
double gridpack::powerflow::PFBus::getPhase()
{
  return *p_vAng_ptr;
}

/**
 * Get generator status
 * @param gen_id generator ID
 * @return vector of generator statuses
 */
bool gridpack::powerflow::PFBus::getGenStatus(std::string gen_id)
{
  int i;
  int gsize = p_gstatus.size();
  for (i=0; i<gsize; i++) {
    if (gen_id == p_gid[i]) {
      return p_gstatus[i];
    }
  }
  return false;
}

/**
 * Check if bus has any online generator
 * @return true if at least one generator is online
 */
bool gridpack::powerflow::PFBus::hasOnlineGenerator()
{
  int gsize = p_gstatus.size();
  for (int i = 0; i < gsize; i++) {
    if (p_gstatus[i] == 1) {
      return true;
    }
  }
  return false;
}

/**
 * Get total capacity (Pmax) of all online generators on this bus
 * @return total Pmax in MW
 */
double gridpack::powerflow::PFBus::getOnlineGenCapacity()
{
  double capacity = 0.0;
  int gsize = p_gstatus.size();
  for (int i = 0; i < gsize; i++) {
    if (p_gstatus[i] == 1) {
      capacity += p_pt[i];
    }
  }
  return capacity;
}

/**
 * Calculate power injection at this bus from connected branches.
 * This updates p_Pinj and p_Qinj for the slack bus.
 * Should be called after power flow solve before checking generator output.
 */
void gridpack::powerflow::PFBus::calculatePowerInjection()
{
  if (!getReferenceBus() && !isIsolated()) return;

  std::vector<boost::shared_ptr<gridpack::component::BaseComponent> > branches;
  getNeighborBranches(branches);
  int size = branches.size();
  double P = 0.0, Q = 0.0, p, q;

  for (int i = 0; i < size; i++) {
    gridpack::powerflow::PFBranch *branch =
      dynamic_cast<gridpack::powerflow::PFBranch*>(branches[i].get());
    branch->getPQ(this, &p, &q);
    P += p;
    Q += q;
  }
  // Also add bus's own Pi, Qi (from shunts)
  P += p_v * p_v * p_ybusr;
  Q += p_v * p_v * (-p_ybusi);

  p_Pinj = P;
  p_Qinj = Q;
}

/**
 * Get total real power output of all generators on this bus (after PF solve)
 * For slack bus, this is calculated from power injection.
 * For PV/PQ buses, this is the sum of scheduled generator outputs.
 * @return total Pgen in MW
 */
double gridpack::powerflow::PFBus::getTotalGenOutput()
{
  double totalPgen = 0.0;
  int gsize = p_gstatus.size();

  if (getReferenceBus()) {
    // For slack bus, calculate from power injection
    // First update p_Pinj if needed
    calculatePowerInjection();
    double pl = 0.0;
    for (int i = 0; i < p_pl.size(); i++) {
      if (p_lstatus[i] == 1) {
        pl += p_pl[i];
      }
    }
    // Add voltage-dependent P load (IP*V + YP*V^2)
    double pzip, qzip;
    getZIPLoadPower(p_v, pzip, qzip);
    pl += pzip;
    totalPgen = p_Pinj * p_sbase + pl;
  } else {
    // For non-slack buses, sum scheduled generator outputs
    for (int i = 0; i < gsize; i++) {
      if (p_gstatus[i] == 1) {
        totalPgen += p_pg[i];
      }
    }
  }
  return totalPgen;
}

/**
 * Check if generator output exceeds capacity on this bus
 * @return true if within limits, false if Pgen > Pmax
 */
bool gridpack::powerflow::PFBus::checkGenCapacity()
{
  double pgen = getTotalGenOutput();
  double pmax = getOnlineGenCapacity();

  // Allow small tolerance for numerical errors
  double tolerance = 0.01 * pmax;  // 1% tolerance
  if (tolerance < 0.1) tolerance = 0.1;  // At least 0.1 MW

  return (pgen <= pmax + tolerance);
}

/**
 * Get list of generator IDs
 * @return vector of generator IDs
 */
std::vector<std::string> gridpack::powerflow::PFBus::getGenerators()
{
  return p_gid;
}

/**
 * Get list of load IDs
 * @return vector of generator IDs
 */
std::vector<std::string> gridpack::powerflow::PFBus::getLoads()
{
  return p_lid;
}

/**
 * Get total reactive power output of all online generators
 * @return total Qgen in MW (system base units)
 */
double gridpack::powerflow::PFBus::getTotalReactiveOutput()
{
  double totalQgen = 0.0;
  for (int i = 0; i < p_ngen; i++) {
    if (p_gstatus[i] == 1) {
      totalQgen += p_qg[i];
    }
  }
  return totalQgen;
}

/**
 * Get real power output for a specific generator
 * @param genId generator ID
 * @return Pgen in MW
 */
double gridpack::powerflow::PFBus::getGenPOutput(const std::string& genId)
{
  for (int i = 0; i < p_ngen; i++) {
    if (p_gid[i] == genId) return p_pg[i];
  }
  return 0.0;
}

/**
 * Get reactive power output for a specific generator
 * @param genId generator ID
 * @return Qgen in MVAr
 */
double gridpack::powerflow::PFBus::getGenQOutput(const std::string& genId)
{
  for (int i = 0; i < p_ngen; i++) {
    if (p_gid[i] == genId) return p_qg[i];
  }
  return 0.0;
}

/**
 * Get maximum reactive power for a specific generator
 * @param genId generator ID
 * @return Qmax in MVAr
 */
double gridpack::powerflow::PFBus::getGenQMax(const std::string& genId)
{
  for (int i = 0; i < p_ngen; i++) {
    if (p_gid[i] == genId) return p_qmax[i];
  }
  return 0.0;
}

/**
 * Get minimum reactive power for a specific generator
 * @param genId generator ID
 * @return Qmin in MVAr
 */
double gridpack::powerflow::PFBus::getGenQMin(const std::string& genId)
{
  for (int i = 0; i < p_ngen; i++) {
    if (p_gid[i] == genId) return p_qmin[i];
  }
  return 0.0;
}

/**
 * Get voltage setpoint for a specific generator
 * @param genId generator ID
 * @return voltage setpoint in pu
 */
double gridpack::powerflow::PFBus::getGenVSetpoint(const std::string& genId)
{
  for (int i = 0; i < p_ngen; i++) {
    if (p_gid[i] == genId) return p_vs[i];
  }
  return 0.0;
}

/**
 * Get IREG (remote regulated bus number) for a generator by index
 * @param idx generator index (0-based)
 * @return IREG bus number (0 = local regulation)
 */
int gridpack::powerflow::PFBus::getIREG(int idx) const
{
  if (idx >= 0 && idx < static_cast<int>(p_ireg.size())) {
    return p_ireg[idx];
  }
  return 0;
}

/**
 * Get generator status by index
 * @param idx generator index (0-based)
 * @return generator status (1 = online, 0 = offline)
 */
int gridpack::powerflow::PFBus::getGenStatusByIdx(int idx) const
{
  if (idx >= 0 && idx < static_cast<int>(p_gstatus.size())) {
    return p_gstatus[idx];
  }
  return 0;
}

/**
 * Get voltage setpoint for a generator by index
 * @param idx generator index (0-based)
 * @return voltage setpoint in pu
 */
double gridpack::powerflow::PFBus::getVSByIdx(int idx) const
{
  if (idx >= 0 && idx < static_cast<int>(p_vs.size())) {
    return p_vs[idx];
  }
  return 1.0;
}

/**
 * Adjust the bus voltage magnitude for remote voltage regulation.
 * Updates both p_v and p_voltage so the adjustment persists across
 * outer iterations.
 * @param dv voltage adjustment in pu
 */
void gridpack::powerflow::PFBus::adjustVoltageForRemoteReg(double dv)
{
  p_v += dv;
  p_voltage += dv;
  // Clamp to bus voltage limits to prevent unbounded accumulation over many
  // outer iterations (each call adds dv; without clamping p_voltage drifts
  // to multi-pu values causing impossible Q requirements)
  if (p_vmax > 0.0 && p_voltage > p_vmax) { p_v = p_vmax; p_voltage = p_vmax; }
  if (p_vmin > 0.0 && p_voltage < p_vmin) { p_v = p_vmin; p_voltage = p_vmin; }
  if (p_vMag_ptr) *p_vMag_ptr = p_v;
}

/**
 * Set generator status
 * @param gen_id generator ID
 * @param status generator status
 * @return true if generator ID found, false otherwise
 */
bool gridpack::powerflow::PFBus::setGenStatus(std::string gen_id, bool status)
{
  int i;
  int gsize = p_gstatus.size();
  for (i=0; i<gsize; i++) {
    if (gen_id == p_gid[i]) {
      // Only modify values if status is actually changing
      // For already-offline generators, calling setGenStatus(false) should be a no-op
      if (p_gstatus[i] == status) {
        return true;  // Status unchanged, but ID was found
      }
      p_gstatus[i] = status;
      if (status == 0) {
        p_pFac[i] = 0.0;
        p_qFac[i] = 0.0;
        p_qmax[i] = 0.0;
        p_qmin[i] = 0.0;
      } else {
        p_pFac[i] = p_pFac_orig[i];
        p_qFac[i] = p_qFac_orig[i];
        p_qmax[i] = p_qmax_orig[i];
        p_qmin[i] = p_qmin_orig[i];
      }

      // Check if any generators are still online on this bus.
      // If not, convert the bus from PV to PQ to avoid inconsistent state
      // where the bus is flagged as PV but has no Q capacity.
      // When turning a generator back on, restore PV status if the bus was originally PV.
      bool hasOnlineGen = false;
      for (int j = 0; j < gsize; j++) {
        if (p_gstatus[j] == 1) {
          hasOnlineGen = true;
          break;
        }
      }

      if (!hasOnlineGen && p_isPV) {
        // Save the original PV status before converting to PQ
        p_saveisPV = p_isPV;
        p_isPV = false;
        if (p_PV_ptr) *p_PV_ptr = false;
      } else if (hasOnlineGen && !p_isPV && p_saveisPV) {
        // Restore PV status if the bus was originally PV and now has an online generator
        p_isPV = true;
        if (p_PV_ptr) *p_PV_ptr = true;
      }

      return true;
    }
  }
  return false;
}

/**
 * Set isPV status
 * @param status isPV status
 */
void gridpack::powerflow::PFBus::setIsPV(int status)
{
  p_saveisPV = p_isPV;
  p_isPV = status;
  if (p_PV_ptr) *p_PV_ptr = status;
  p_v = p_voltage;
}

/**
 * Reset isPV status
 */
void gridpack::powerflow::PFBus::resetIsPV()
{
  p_isPV = p_saveisPV;
  if (p_PV_ptr) *p_PV_ptr = p_saveisPV;
}

/**
 * setSBus
 * BUS = (CG*(GEN(ON,PG) + J*GEN(ON,QG)-(PD+J*QD))/BASEMVA
 */
void gridpack::powerflow::PFBus::setSBus(void)
{
  int i;
  double pg, qg, pl, ql;
  pg = 0.0;
  qg = 0.0;
  pl = 0.0;
  ql = 0.0;
  bool usegen = false;
  for (i=0; i<p_gstatus.size(); i++) {
    if (p_gstatus[i] == 1) {
      pg += p_pg[i];
      qg += p_qg[i];
      usegen = true;
    }
  }
  for (i=0; i<p_lstatus.size(); i++) {
    if (p_lstatus[i] == 1) {
      pl += p_pl[i];
      ql += p_ql[i];
    }
  }
  if (p_gstatus.size() > 0 && usegen) {
    gridpack::ComplexType sBus((pg - pl) / p_sbase, (qg - ql) / p_sbase);
    p_P0 = real(sBus);
    p_Q0 = imag(sBus);
  } else {
    gridpack::ComplexType sBus((- pl) / p_sbase, (- ql) / p_sbase);
    p_P0 = real(sBus);
    p_Q0 = imag(sBus);
  } 
}

/**
 ** Update pg of specified bus element based on their genID
 ** @param busID
 ** @param genID
 ** @param value
 **/
/*
void gridpack::powerflow::PFBus::updatePg(int busID, std::string genID, double value)
{
  if (getOriginalIndex() == busID) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        if (p_gid[i] == genID) {
          p_pg[i] += value;
        }
      }
    }
  }
}
*/

/**
 * Write output from buses to standard out
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool gridpack::powerflow::PFBus::serialWrite(char *string, const int bufsize,
                                             const char *signal)
{
  if (signal == NULL) {
    double pi = 4.0*atan(1.0);
    double angle = p_a*180.0/pi;
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    getNeighborBranches(branches);
    int size = branches.size();
    // Also add bus i's own Pi, Qi
    double P, Q, p, q;
    P = 0.0;
    Q = 0.0;
    int i;
    for (i=0; i<size; i++) {
      gridpack::powerflow::PFBranch *branch
        = dynamic_cast<gridpack::powerflow::PFBranch*>(branches[i].get());
      branch->getPQ(this, &p, &q);
      P += p;
      Q += q;
    }
    P += p_v*p_v*p_ybusr;
    Q += p_v*p_v*(-p_ybusi);
    p_Pinj = P;
    p_Qinj = Q;
    if (!isIsolated() && !p_isStarBus) {
      sprintf(string, "     %6d      %12.6f         %12.6f         %12.6f         %12.6f\n",
            getOriginalIndex(),angle,p_v,p_Pinj,p_Qinj);
    }
/*
  if (signal == NULL) {
    double pi = 4.0*atan(1.0);
    double angle = p_a*180.0/pi;
    if (!isIsolated()) {
      sprintf(string, "     %6d      %12.6f         %12.6f\n",
            getOriginalIndex(),angle,p_v);
    } else {
      return false;
    } */
  } else if (!strcmp(signal,"vr_str")) {
    double pi = 4.0*atan(1.0);
    double angle = p_a*180.0/pi;
    int use_vmag = 1;
    if (p_saveisPV || p_original_isolated) use_vmag = 0;
    int changed = 0;
    if (p_isPV != p_saveisPV) changed = 1;
    sprintf(string, "%6d %20.12e %20.12e %d %d\n",
        getOriginalIndex(),angle,p_v,use_vmag,changed);
  } else if (!strcmp(signal,"vfail_str")) {
    int use_vmag = 1;
    if (p_saveisPV || p_original_isolated) use_vmag = 0;
    int changed = 0;
    sprintf(string, "%6d %20.12e %20.12e %d %d\n",
        getOriginalIndex(),0.0,0.0,use_vmag,changed);
  } else if (!strcmp(signal,"ca")) {
    double pi = 4.0*atan(1.0);
    double angle = p_a*180.0/pi;
    bool found = false;
    if ((p_v > p_vmax) || (p_v < p_vmin) ) {
      found = true;
      if (!isIsolated()) {
        sprintf(string, "     %6d      %12.6f         %12.6f\n",
            getOriginalIndex(),angle,p_v);
      } else {
        sprintf(string, "     %6d      %12.6f         %12.6f\n",
            getOriginalIndex(),0.0,0.0);
      }
    }
    return found;
  } else if (!strcmp(signal,"pq")) {
    gridpack::ComplexType v[2];
    vectorValues(v);
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    getNeighborBranches(branches);
    if (!isIsolated()) {
      sprintf(string, "     %6d      %12.6f      %12.6f      %2d\n",
            getOriginalIndex(),real(v[0]),real(v[1]),
            static_cast<int>(branches.size()));
    } else {
      sprintf(string, "     %6d      %12.6f      %12.6f      %2d\n",
            getOriginalIndex(),0.0, 0.0,
            static_cast<int>(branches.size()));
    }
  } else if (!strcmp(signal,"record")) {
    if (p_isStarBus) return false;
    char sbuf[128];
    char *cptr = string;
    int slen = 0;
    int nld = p_nload;
    int i;
    double pl =0.0;
    double ql =0.0;
    for (i=0; i<nld; i++) {
      if (p_lstatus[i]) pl += p_pl[i];
      if (p_lstatus[i]) ql += p_ql[i];
    }
    // Add voltage-dependent ZIP load at solved voltage
    {
      double pzip, qzip;
      getZIPLoadPower(p_v, pzip, qzip);
      pl += pzip;
      ql += qzip;
    }
    sprintf(sbuf,"%8d, %4d, %16.8f, %16.8f,",getOriginalIndex(),p_type,
           pl/p_sbase,ql/p_sbase);
    int len = strlen(sbuf);
    if (len<=bufsize) {
      sprintf(cptr,"%s",sbuf);
      slen += len;
      cptr += len;
    }
    double pgen = 0.0;
    double qgen = 0.0;
    double qmin = 0.0;
    double qmax = 0.0;
    int ngen = p_ngen;
    for (i=0; i<ngen; i++) {
      if (p_gstatus[i]) pgen += p_pg[i];
      if (p_gstatus[i]) qgen += p_qg[i];
      if (p_gstatus[i]) qmin += p_qmin[i];
      if (p_gstatus[i]) qmax += p_qmax[i];
    }
    sprintf(sbuf," %16.8f, %16.8f, %16.8f, %16.8f,",
            pgen/p_sbase,qgen/p_sbase,qmax/p_sbase,qmin/p_sbase);
    len = strlen(sbuf);
    if (slen+len<=bufsize) {
      sprintf(cptr,"%s",sbuf);
      slen += len;
      cptr += len;
    }
    int gstatus = 0;
    double pt = 0.0;
    double pb = 0.0;
    for (i=0; i<ngen; i++) {
      if (p_gstatus[i]) gstatus = 1;
      if (p_gstatus[i]) pt += p_pt[i];
      if (p_gstatus[i]) pb += p_pb[i];
    }
    sprintf(sbuf," %16.8f, %16.8f, %1d,",pt,pb,gstatus);
    len = strlen(sbuf);
    if (slen+len<=bufsize) {
      sprintf(cptr,"%s",sbuf);
      slen += len;
      cptr += len;
    }
    double gl, bl;
    YMBus::getShuntValues(&bl, &gl);
    int area;
    p_data->getValue(BUS_AREA,&area);
    sprintf(sbuf," %16.8f, %16.8f, %8d,",gl,bl,area);
    len = strlen(sbuf);
    if (slen+len<=bufsize) {
      sprintf(cptr,"%s",sbuf);
      slen += len;
      cptr += len;
    }
    double zero = 0.0;
    int nzone;
    double basekv;
    p_data->getValue(BUS_ZONE,&nzone);
    p_data->getValue(BUS_BASEKV,&basekv);
    double pi = 4.0*atan(1.0);
    double angle = p_a*180.0/pi;
    if (!isIsolated()) {
      sprintf(sbuf," %16.8f, %16.8f, %16.8f, %8d, %4.2f, %4.2f\n",p_v,angle,basekv,nzone,p_vmax,p_vmin);
    } else {
      sprintf(sbuf," %16.8f, %16.8f, %16.8f, %8d, %4.2f, %4.2f\n",0.0,0.0,basekv,nzone,p_vmax,p_vmin);
    }
    len = strlen(sbuf);
    if (slen+len<=bufsize) {
      sprintf(cptr,"%s",sbuf);
      slen += len;
      cptr += len;
    }
  } else if (!strcmp(signal,"power") || !strcmp(signal,"gen_str")) {
    char sbuf[128];
    char *cptr = string;
    int i, len, slen = 0;
    int ngen=p_pFac.size();
    // Evalate p_Pinj and p_Qinj if bus is reference bus. This is skipped when
    // evaluating matrix elements.
#ifndef LARGE_MATRIX
    if (getReferenceBus() || isIsolated()) {
      std::vector<boost::shared_ptr<BaseComponent> > branches;
      getNeighborBranches(branches);
      int size = branches.size();
      double P, Q, p, q;
      P = 0.0;
      Q = 0.0;
      for (i=0; i<size; i++) {
        gridpack::powerflow::PFBranch *branch
          = dynamic_cast<gridpack::powerflow::PFBranch*>(branches[i].get());
        branch->getPQ(this, &p, &q);
        P += p;
        Q += q;
      }
      // Also add bus i's own Pi, Qi
      P += p_v*p_v*p_ybusr;
      Q += p_v*p_v*(-p_ybusi);
      p_Pinj = P;
      p_Qinj = Q;
    }
#endif
    double pl =0.0;
    double ql =0.0;
    for (i=0; i<p_pl.size(); i++) {
      if (p_lstatus[i] == 1) {
        pl += p_pl[i];
        ql += p_ql[i];
      }
    }
    // Add voltage-dependent ZIP load at solved voltage
    {
      double pzip, qzip;
      getZIPLoadPower(p_v, pzip, qzip);
      pl += pzip;
      ql += qzip;
    }
    for (i=0; i<ngen; i++) {
      double pval = 0.0;
      double qval = 0.0;

      if(getReferenceBus()) {
	pval = p_pFac[i]*(p_Pinj*p_sbase+pl);
	// Use p_qFac for Q distribution (based on Qmax capability)
	qval = p_qFac[i]*(p_Qinj*p_sbase+ql);
      } else if (p_isPV) {
	if(p_gstatus[i]) {
	  pval = p_pg[i];
	  if (p_qlim) {
	    // When qlim=true, use p_qg which is set by chkQlim() with RMPCT-based distribution
	    qval = p_qg[i];
	  } else {
	    // When qlim=false, chkQlim() is not called, so calculate Q from p_Qinj
	    // Use p_qFac (Qmax-based) to distribute Q among generators with reactive capability
	    // Generators with Qmax=Qmin=0 get p_qFac=0, so they won't be assigned Q they can't provide
	    qval = p_qFac[i]*(p_Qinj*p_sbase+ql);
	  }
	} else {
	  pval = 0.0;
	  qval = 0.0;
	}
      } else {
	if(p_gstatus[i]) {
	  pval = p_pg[i];
	  qval = p_qg[i];
	} else {
	  pval = 0.0;
	  qval = 0.0;
	}
      }
      if (!strcmp(signal,"power")) {
        sprintf(sbuf, "     %6d      %s   %12.6f      %12.6f\n",
            getOriginalIndex(),p_gid[i].c_str(),pval,qval);
      } else {
        sprintf(sbuf, "%6d %s %20.12e %20.12e\n",
            getOriginalIndex(),p_gid[i].c_str(),pval,qval);
      }
      len = strlen(sbuf);
      if (slen+len<=bufsize) {
        sprintf(cptr,"%s",sbuf);
        slen += len;
        cptr += len;
      }
    }
    if (slen>0) {
      return true;
    } else {
      return false;
    }
  } else if (!strcmp(signal,"pfail_str")) {
    char sbuf[128];
    char *cptr = string;
    int i, len, slen = 0;
    int ngen=p_pFac.size();
    for (i=0; i<ngen; i++) {
      sprintf(sbuf, "%6d %s %20.12e %20.12e\n",
            getOriginalIndex(),p_gid[i].c_str(),0.0,0.0);
      len = strlen(sbuf);
      if (slen+len<=bufsize) {
        sprintf(cptr,"%s",sbuf);
        slen += len;
        cptr += len;
      }
    }
    if (slen>0) {
      return true;
    } else {
      return false;
    }
  } else if (!strcmp(signal,"src_gen")) {
    if (p_source) {
      char sbuf[128];
      char *cptr = string;
      int i, len, slen = 0;
      std::string status;
      for (i=0; i<p_ngen; i++) {
        if (p_gstatus[i]) {
          status = "  active";
        } else {
          status = "inactive";
        }
        sprintf(sbuf,"%8d %s %s %4d %4d %14.4f %14.4f %14.4f %14.4f\n",
              getOriginalIndex(),
              p_gid[i].c_str(),status.c_str(),p_area,p_zone,p_savePg[i],
              p_pg[i],p_pb[i],p_pt[i]);
        len = strlen(sbuf);
        if (slen+len <= bufsize) {
          sprintf(cptr,"%s",sbuf);
          slen += len;
          cptr += len;
        }
      }
      if (slen>0) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  } else if (!strcmp(signal,"sink_load")) {
    if (p_sink) {
      char sbuf[128];
      char *cptr = string;
      int i, len, slen = 0;
      std::string status;
      for (i=0; i<p_nload; i++) {
        if (p_lstatus[i]) {
          status = "  active";
        } else {
          status = "inactive";
        }
        sprintf(sbuf,"%8d %s %s %4d %4d %14.4f %14.4f %14.4f %14.4f\n",
              getOriginalIndex(),
              p_lid[i].c_str(),status.c_str(),p_area,p_zone,p_savePl[i],
              p_pl[i],p_saveQl[i],p_ql[i]);
        len = strlen(sbuf);
        if (slen+len <= bufsize) {
          sprintf(cptr,"%s",sbuf);
          slen += len;
          cptr += len;
        }
      }
      if (slen>0) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  } else if (!strcmp(signal,"swshunt")) {
    if (p_hasSwitchedShunt) {
      sprintf(string, "     %6d  MODSW=%d  Binit=%8.2f  Bcurrent=%8.2f  V=%8.4f  [%6.4f, %6.4f]\n",
              getOriginalIndex(), p_swshunt_modsw,
              p_swshunt_binit, p_swshunt_bcurrent,
              p_v, p_swshunt_vswlo, p_swshunt_vswhi);
      return true;
    }
    return false;
  }
  return true;
}

/**
 * Return the complex voltage on this bus
 * @return the complex voltage
 */
gridpack::ComplexType gridpack::powerflow::PFBus::getComplexVoltage(void)
{
  p_a = *p_vAng_ptr;
  p_v =  *p_vMag_ptr;
  gridpack::ComplexType ret(cos(p_a),sin(p_a));
  ret = ret*p_v;
  return ret;
}

/**
 * Save state variables inside the component to a DataCollection object.
 * This can be used as a way of moving data in a way that is useful for
 * creating output or for copying state data from one network to another.
 * @param data data collection object into which new values are inserted
 */
void gridpack::powerflow::PFBus::saveData(
    boost::shared_ptr<gridpack::component::DataCollection> data)
{
  double rval;
  int i;
  if (!data->setValue("BUS_PF_VMAG",*p_vMag_ptr)) {
    data->addValue("BUS_PF_VMAG",*p_vMag_ptr);
  }
  rval = *p_vAng_ptr;
  double pi = 4.0*atan(1.0);
  rval = 180.0*rval/pi;
  if (!data->setValue("BUS_PF_VANG",rval)) {
    data->addValue("BUS_PF_VANG",rval);
  }
  if (!data->setValue("BUS_TYPE",p_type)) {
    data->addValue("BUS_TYPE",p_type);
  }
  int ngen=p_pFac.size();
  // Evalate p_Pinj and p_Qinj if bus is reference bus. This is skipped when
  // evaluating matrix elements.
#ifndef LARGE_MATRIX
  if (getReferenceBus() || isIsolated()) {
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    getNeighborBranches(branches);
    int size = branches.size();
    double P, Q, p, q;
    P = 0.0;
    Q = 0.0;
    for (i=0; i<size; i++) {
      gridpack::powerflow::PFBranch *branch
        = dynamic_cast<gridpack::powerflow::PFBranch*>(branches[i].get());
      branch->getPQ(this, &p, &q);
      P += p;
      Q += q;
    }
    // Also add bus i's own Pi, Qi
    P += p_v*p_v*p_ybusr;
    Q += p_v*p_v*(-p_ybusi);
    p_Pinj = P;
    p_Qinj = Q;
  }
#endif
  double pl=0.0;
  double ql=0.0;
  for (i=0; i<p_pl.size(); i++) {
    if (p_lstatus[i] == 1) {
      pl += p_pl[i];
      ql += p_ql[i];
    }
  }
  // Add voltage-dependent ZIP load at solved voltage
  {
    double pzip, qzip;
    getZIPLoadPower(p_v, pzip, qzip);
    pl += pzip;
    ql += qzip;
  }
  for (i=0; i<ngen; i++) {
    if(getReferenceBus()) {
      rval = p_pFac[i]*(p_Pinj+pl/p_sbase);
      if (!data->setValue("GENERATOR_PF_PGEN",rval,i)) {
	data->addValue("GENERATOR_PF_PGEN",rval,i);
      }
      // Use p_qFac for Q distribution (based on Qmax capability)
      rval = p_qFac[i]*(p_Qinj+ql/p_sbase);
      if (!data->setValue("GENERATOR_PF_QGEN",rval,i)) {
	data->addValue("GENERATOR_PF_QGEN",rval,i);
      }
    } else if (p_isPV) {
      if (p_gstatus[i]) {
        if (!data->setValue("GENERATOR_PF_PGEN",p_pg[i]/p_sbase,i)) {
	  data->addValue("GENERATOR_PF_PGEN",p_pg[i]/p_sbase,i);
        }
        if (p_qlim) {
          // When qlim=true, use p_qg which is set by chkQlim() with RMPCT-based distribution
          rval = p_qg[i]/p_sbase;
        } else {
          // When qlim=false, chkQlim() is not called, so calculate Q from p_Qinj
          // Use p_qFac (Qmax-based) for generators with reactive capability
          rval = p_qFac[i]*(p_Qinj+ql/p_sbase);
        }
        if (!data->setValue("GENERATOR_PF_QGEN",rval,i)) {
	  data->addValue("GENERATOR_PF_QGEN",rval,i);
        }
      } else {
        if (!data->setValue("GENERATOR_PF_PGEN",0.0,i)) {
	  data->addValue("GENERATOR_PF_PGEN",0.0,i);
        }
        if (!data->setValue("GENERATOR_PF_QGEN",0.0,i)) {
	  data->addValue("GENERATOR_PF_QGEN",0.0,i);
        }
      }
    } else {
      if (p_gstatus[i]) {
        if (!data->setValue("GENERATOR_PF_PGEN",p_pg[i]/p_sbase,i)) {
	  data->addValue("GENERATOR_PF_PGEN",p_pg[i]/p_sbase,i);
        }
      } else {
        if (!data->setValue("GENERATOR_PF_PGEN",0.0,i)) {
	  data->addValue("GENERATOR_PF_PGEN",0.0,i);
        }
      }

      if (p_gstatus[i]) {
        if (!data->setValue("GENERATOR_PF_QGEN",p_qg[i]/p_sbase,i)) {
	  data->addValue("GENERATOR_PF_QGEN",p_qg[i]/p_sbase,i);
        }
      } else {
        if (!data->setValue("GENERATOR_PF_QGEN",0.0,i)) {
	  data->addValue("GENERATOR_PF_QGEN",0.0,i);
        }
      }
    }
  }
}

/**
 * Save state variables inside the component to a DataCollection object.
 * This can be used as a way of moving data in a way that is useful for
 * creating output or for copying state data from one network to another.
 * @param data data collection object into which new values are inserted
 * added by Renke, also modify the original bus mag, ang,
 * and the original generator PG QG in the datacollection
 */
void gridpack::powerflow::PFBus::saveDataAlsotoOrg(
    boost::shared_ptr<gridpack::component::DataCollection> data)
{
  double rval;
  int i;
  if (!data->setValue("BUS_PF_VMAG",*p_vMag_ptr)) {
    data->addValue("BUS_PF_VMAG",*p_vMag_ptr);
  }
  data->setValue(BUS_VOLTAGE_MAG,*p_vMag_ptr);  //also modify the original BUS_VOLTAGE_MAG 
  
  rval = *p_vAng_ptr;
  double pi = 4.0*atan(1.0);
  rval = 180.0*rval/pi;
  if (!data->setValue("BUS_PF_VANG",rval)) {
    data->addValue("BUS_PF_VANG",rval);
  }
  data->setValue(BUS_VOLTAGE_ANG,rval); //also modify the original BUS_VOLTAGE_ANG 
  
  
  if (!data->setValue("BUS_TYPE",p_type)) {
    data->addValue("BUS_TYPE",p_type);
  }
  int ngen=p_pFac.size();
  // Evalate p_Pinj and p_Qinj if bus is reference bus. This is skipped when
  // evaluating matrix elements.
#ifndef LARGE_MATRIX
  if (getReferenceBus() || isIsolated()) {
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    getNeighborBranches(branches);
    int size = branches.size();
    double P, Q, p, q;
    P = 0.0;
    Q = 0.0;
    for (i=0; i<size; i++) {
      gridpack::powerflow::PFBranch *branch
        = dynamic_cast<gridpack::powerflow::PFBranch*>(branches[i].get());
      branch->getPQ(this, &p, &q);
      P += p;
      Q += q;
    }
    // Also add bus i's own Pi, Qi
    P += p_v*p_v*p_ybusr;
    Q += p_v*p_v*(-p_ybusi);
    p_Pinj = P;
    p_Qinj = Q;
  }
#endif
  double pl=0.0;
  double ql=0.0;
  for (i=0; i<p_pl.size(); i++) {
    if (p_lstatus[i] == 1) {
      pl += p_pl[i];
      ql += p_ql[i];
    }
  }
  // Add voltage-dependent ZIP load at solved voltage
  {
    double pzip, qzip;
    getZIPLoadPower(p_v, pzip, qzip);
    pl += pzip;
    ql += qzip;
  }
  for (i=0; i<ngen; i++) {
    rval = p_pFac[i]*(p_Pinj+pl/p_sbase);
    if (!data->setValue("GENERATOR_PF_PGEN",rval,i)) {
      data->addValue("GENERATOR_PF_PGEN",rval,i);
    }
	data->setValue(GENERATOR_PG,rval*100.0,i); //also modify the original GENERATOR_PG

    // Use p_qFac for Q distribution (based on Qmax capability)
    rval = p_qFac[i]*(p_Qinj+ql/p_sbase);
    if (!data->setValue("GENERATOR_PF_QGEN",rval,i)) {
      data->addValue("GENERATOR_PF_QGEN",rval,i);
    }
	data->setValue(GENERATOR_QG,rval*100.0,i); //also modify the original GENERATOR_QG

  }
}

/**
 * Modify parameters inside the bus module. This is designed to be
 * extensible
 * @param name character string describing parameter to be modified
 * @param busID generator bus number
 * @param genID specified genID
 * @param value new value of parameter
 */
void gridpack::powerflow::PFBus::setParam(std::string name, int busID, 
    std::string genID, double value) 
{
  if (getOriginalIndex() == busID) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        if (p_gid[i] == genID) {
          if (name == GENERATOR_PG) {
            p_pg[i] = value;
          } else if (name == GENERATOR_QG) {
            p_qg[i] = value;
          }
        }
      }
    }
  }
}

/**
 * Modify parameters inside the bus module. This is designed to be
 * extensible
 * @param name character string describing parameter to be modified
 * @param busID generator bus number
 * @param genID specified genID
 * @param value new value of parameter
 */
void gridpack::powerflow::PFBus::setParam(int busID, 
    std::string genID, double value) 
{
//  if (name == GENERATOR_PG) {
   if (getOriginalIndex() == busID) {
     if (p_ngen > 0) {
       for (int i = 0; i < p_ngen; i++) {
         if (p_gid[i] == genID) {
           p_pg[i] += value;
         }
       }
     }
   }
   // if (idx >= 0 && idx<p_pg.size()) {
   //   p_pg[idx] = value;
   // }
//  } else if (name == GENERATOR_QG) {
   // if (idx >= 0 && idx<p_qg.size()) {
   //   p_qg[idx] = value;
   // }
//  }
}

/**
 * Access parameters inside the bus module. This is designed to be
 * extensible
 * @param name character string describing parameter to be accessed
 * @param value value of parameter
 * @param idx index (if necessary) of variable to be accessed
 */
void gridpack::powerflow::PFBus::getParam(std::string &name,
    double *value, int idx)
{
  if (name == GENERATOR_PG) {
    if (idx >= 0 && idx<p_pg.size()) {
      *value = p_pg[idx];
    }
  } else if (name == GENERATOR_QG) {
    if (idx >= 0 && idx<p_qg.size()) {
      *value = p_qg[idx];
    }
  }
}

void gridpack::powerflow::PFBus::getParam(std::string &name,
    int *value, int idx)
{
  if (name == GENERATOR_NUMBER) {
    *value = p_pg.size();
  }
}

/**
 * Get index of internal bus element based on character string identifier
 * @param name character string describing element
 * @param tag character string specifying bus element
 * @return index of element
 */
int gridpack::powerflow::PFBus::getElementIndex(std::string &name, std::string &tag)
{
  if (name == "GENERATOR") {
    int i;
    int nsize = static_cast<int>(p_gid.size());
    for (i=0; i<nsize; i++) {
      if (tag == p_gid[i]) {
        return i;
      }
    }
  }
  return -1;
}

/**
 * Set parameter to ignore voltage violations
 * @param flag value of ignore parameter
 */
void gridpack::powerflow::PFBus::setIgnore(bool flag)
{
  p_ignore = flag;
}

/**
 * Get parameter to ignore voltage violations
 * @return value of ignore parameter
 */
bool gridpack::powerflow::PFBus::getIgnore()
{
  return p_ignore;
}

/**
 * Get area parameter for bus
 * @return bus area index
 */
int gridpack::powerflow::PFBus::getArea()
{
  return p_area;
}

/**
 * Get zone parameter for bus
 * @return bus zone index
 */
int gridpack::powerflow::PFBus::getZone()
{
  return p_zone;
}

/**
 * Get owner number for bus
 * @return bus owner number (0 if not set)
 */
int gridpack::powerflow::PFBus::getOwner()
{
  int owner = 0;
  if (p_data) p_data->getValue(BUS_OWNER, &owner);
  return owner;
}

/**
 * Get base voltage for bus in kV
 * @return base kV (0.0 if not set)
 */
double gridpack::powerflow::PFBus::getBaseKV()
{
  double basekv = 0.0;
  if (p_data) p_data->getValue(BUS_BASEKV, &basekv);
  return basekv;
}

/**
 * Get bus name string
 * @return bus name (empty if not set)
 */
std::string gridpack::powerflow::PFBus::getBusName()
{
  std::string name;
  if (p_data) p_data->getValue(BUS_NAME, &name);
  return name;
}

/**
 * Evaluate diagonal block of Jacobian for power flow calculation and
 * return result as an array of real values
 * @param rvals values of Jacobian block
 * @return number of values returned
 */
int gridpack::powerflow::PFBus::diagonalJacobianValues(double *rvals)
{
  if (!isIsolated()) {
    // Compute ZIP Jacobian contributions: d(pzip)/dV and d(qzip)/dV
    double dpzip_dV = 0.0, dqzip_dV = 0.0;
    for (int i = 0; i < p_lstatus.size(); i++) {
      if (p_lstatus[i] == 1) {
        dpzip_dV += p_ip[i] + 2.0 * p_yp[i] * p_v;
        dqzip_dV += p_iq[i] - 2.0 * p_yq[i] * p_v;
      }
    }
#ifdef LARGE_MATRIX
    if (!getReferenceBus()) {
      rvals[0] = -p_Qinj - p_ybusi * p_v *p_v;
      rvals[1] = p_Pinj - p_ybusr * p_v *p_v;
      rvals[2] = p_Pinj / p_v + p_ybusr * p_v;
      rvals[3] = p_Qinj / p_v - p_ybusi * p_v;
      // Fix up matrix elements if bus is PV bus
      if (p_isPV) {
        rvals[1] = 0.0;
        rvals[2] = 0.0;
        rvals[3] = 1.0;
      } else {
        // Add ZIP load derivatives for PQ buses
        rvals[2] += dpzip_dV / p_sbase;
        rvals[3] += dqzip_dV / p_sbase;
      }
      return 4;
    } else {
      rvals[0] = 1.0;
      rvals[1] = 0.0;
      rvals[2] = 0.0;
      rvals[3] = 1.0;
      return 4;
    }
#else
    if (!getReferenceBus() && !p_isPV) {
      // Standard PQ
      rvals[0] = -p_Qinj - p_ybusi * p_v *p_v;
      rvals[1] = p_Pinj - p_ybusr * p_v *p_v;
      rvals[2] = p_Pinj / p_v + p_ybusr * p_v;
      rvals[3] = p_Qinj / p_v - p_ybusi * p_v;
      // Add ZIP load derivatives for PQ buses
      rvals[2] += dpzip_dV / p_sbase;
      rvals[3] += dqzip_dV / p_sbase;
      return 4;
    } else if (!getReferenceBus() && p_isPV && p_isIREG_PV) {
      // IREG PV: equations [ΔP, V_remote-VS], variables [θ, V]
      rvals[0] = -p_Qinj - p_ybusi * p_v * p_v;  // ∂ΔP/∂θ
      rvals[1] = 0.0;                              // ∂(V_remote-VS)/∂θ
      rvals[2] = p_Pinj / p_v + p_ybusr * p_v;    // ∂ΔP/∂V
      rvals[2] += dpzip_dV / p_sbase;
      rvals[3] = 1.0e-10;                            // small pivot helper
      return 4;
    } else if (!getReferenceBus() && p_isPV) {
      // Standard PV
      rvals[0] = -p_Qinj - p_ybusi * p_v *p_v;
      return 1;
    } else {
      return 0;
    }
#endif
  } else {
    return 0;
  }
}

/**
 * Push p_isPV values from exchange buffer to p_isPV variable
 */
void gridpack::powerflow::PFBus::pushIsPV()
{
  if (p_PV_ptr) p_isPV = *p_PV_ptr;
}


/**
 * Evaluate RHS values for powerflow equation and return result as
 * an array of real values
 * @param rvals values of Jacobian block
 * @return number of values returned
 */
int gridpack::powerflow::PFBus::rhsValues(double *rvals)
{
  if (!isIsolated()) {
    if (!getReferenceBus()) {
      int nvals;
      std::vector<boost::shared_ptr<BaseComponent> > branches;
      getNeighborBranches(branches);
      int size = branches.size();
      int i;
      double P, Q, p, q;
      P = 0.0;
      Q = 0.0;
      for (i=0; i<size; i++) {
        gridpack::powerflow::PFBranch *branch
          = dynamic_cast<gridpack::powerflow::PFBranch*>(branches[i].get());
        branch->getPQ(this, &p, &q);
        P += p;
        Q += q;
      }
      // Also add bus i's own Pi, Qi
      P += p_v*p_v*p_ybusr;
      Q += p_v*p_v*(-p_ybusi);
      p_Pinj = P;
      p_Qinj = Q;
      P -= p_P0;
      Q -= p_Q0;
      // Add voltage-dependent ZIP load contribution
      double pzip, qzip;
      getZIPLoadPower(p_v, pzip, qzip);
      P += pzip / p_sbase;
      Q += qzip / p_sbase;
      rvals[0] = P;
#ifdef LARGE_MATRIX
      if (!p_isPV) {
        rvals[1] = Q;
      } else {
        rvals[1] = 0.0;
      }
      nvals = 2;
#else
      nvals = 1;
      if (!p_isPV) {
        rvals[1] = Q;
        nvals = 2;
      } else if (p_isIREG_PV) {
        // V_remote - VS constraint
        double v_remote = (p_ireg_remote_v_ptr) ? *p_ireg_remote_v_ptr : p_ireg_vs;
        rvals[1] = v_remote - p_ireg_vs;
        nvals = 2;
      }
#endif
      return nvals;
    } else {
#ifdef LARGE_MATRIX
      std::vector<boost::shared_ptr<BaseComponent> > branches;
      getNeighborBranches(branches);
      int size = branches.size();
      int i;
      double P, Q, p, q;
      P = 0.0;
      Q = 0.0;
      for (i=0; i<size; i++) {
        gridpack::powerflow::PFBranch *branch
          = dynamic_cast<gridpack::powerflow::PFBranch*>(branches[i].get());
        branch->getPQ(this, &p, &q);
        P += p;
        Q += q;
      }
      // Also add bus i's own Pi, Qi
      P += p_v*p_v*p_ybusr;
      Q += p_v*p_v*(-p_ybusi);
      p_Pinj = P;
      p_Qinj = Q;
      rvals[0] = 0.0;
      rvals[1] = 0.0;
      return 2;
#else
      return 0;
#endif
    }
  } else {
    return false;
  }
}

/**
 * Get vector containing generator participation
 * @return vector of generator participation factors
 */
std::vector<double> gridpack::powerflow::PFBus::getGeneratorParticipation()
{
  return p_pFac;
}

/**
 * Set value of real power on individual generators
 * @param tag generator ID
 * @param value new value of real power
 * @param data data collection object associated with bus
 */
void gridpack::powerflow::PFBus::setGeneratorRealPower(
    std::string tag, double value, gridpack::component::DataCollection *data)
{
  int i, idx;
  idx = -1;
  for (i=0; i<p_ngen; i++) {
    if (p_gid[i] == tag) {
      idx = i;
      break;
    }
  }
  if (idx != -1) {
    p_pg[idx] = value;
    if (!data->setValue(GENERATOR_PG,value,idx)) {
      data->addValue(GENERATOR_PG,value,idx);
    }
  } else {
    printf("setGeneratorRealPower: No generator found on bus %d with id: (%s)\n",getOriginalIndex(),tag.c_str());
  }
}

/**
 * Set generator status
 * @param tag generator ID
 * @param status new value of status
 * @param data data collection object associated with
 bus
 */
void gridpack::powerflow::PFBus::setGeneratorStatus(
    std::string tag, int status, gridpack::component::DataCollection *data)
{
  int i, idx;
  idx = -1;
  for (i=0; i<p_ngen; i++) {
    if (p_gid[i] == tag) {
      idx = i;
      break;
    }
  }
  if (idx != -1) {
    p_gstatus[idx] = status;
    if (!data->setValue(GENERATOR_STAT,status,idx)) {
      data->addValue(GENERATOR_STAT,status,idx);
    }
    if(!status) {
      p_pg[idx] = 0.0;
      p_qg[idx] = 0.0;
      data->setValue(GENERATOR_PG,0.0,idx);
      data->setValue(GENERATOR_QG,0.0,idx);
    }
    // Check PV status and change it if all generators are off
    p_isPV = p_gstatus[0];
    for (i=1; i<p_ngen; i++) {
      p_isPV = p_isPV || p_gstatus[i];
    }
    setIsPV(p_isPV);
  } else {
    printf("setGenertorStatus: No generator found on bus %d with id: (%s)\n",getOriginalIndex(),tag.c_str());
  }
}

/**
 * Scale value of real power on all generators
 * @param character ID for generator
 * @param value scale factor for real power
 */
void gridpack::powerflow::PFBus::scaleGeneratorRealPower(std::string tag,
    double value)
{
  int i;
  for (i=0; i<p_ngen; i++) {
    if (p_gid[i] == tag && p_gstatus[i]) {
      if (value > 0.0) {
        double excess = p_pt[i]-p_pg[i];
        if (excess < 0.0) {
          printf("bus: %d generator: %s excess (pt): %f (pg): %f\n",
              getOriginalIndex(),tag.c_str(),p_pt[i],p_pg[i]);
        }
        p_pg[i] += value*excess;
      } else {
        double slack = p_pg[i]-p_pb[i];
        p_pg[i] += value*slack;
      }
      break;
    }
  }
}

/**
 * Set value of real power on individual generators
 * @param tag generator ID
 * @param value new value of real power
 * @param data data collection object associated with bus
 */
void gridpack::powerflow::PFBus::setLoadRealPower(
    std::string tag, double value, gridpack::component::DataCollection *data)
{
  int i, idx;
  idx = -1;
  for (i=0; i<p_nload; i++) {
    if (p_lid[i] == tag) {
      idx = i;
      break;
    }
  }
  if (idx != -1) {
    if (!data->setValue(LOAD_PL,value,idx)) {
      data->addValue(LOAD_PL,value,idx);
    }
    p_pl[idx] = value;
  } else {
    printf("No load found for tag: (%s)\n",tag.c_str());
  }
}

/**
 * Scale value of real and reactive power on loads
 * @param character ID for load
 * @param value scale factor for real power
 */
void gridpack::powerflow::PFBus::scaleLoadPower(std::string tag, double value)
{
  int i;
  for (i=0; i<p_nload; i++) {
    if (p_lid[i] == tag && p_lstatus[i] == 1) {
      p_pl[i] = value*p_pl[i];
      p_ql[i] = value*p_ql[i];
      p_ip[i] = value*p_ip[i];
      p_iq[i] = value*p_iq[i];
      p_yp[i] = value*p_yp[i];
      p_yq[i] = value*p_yq[i];
      break;
    }
  }
}

/**
 * Reset power for generators and loads back to original values
 */
void gridpack::powerflow::PFBus::resetPower()
{
  resetVoltage();
  int i;
  for (i=0; i<p_ngen; i++) {
    p_pg[i] = p_savePg[i];
  }
  for (i=0; i<p_nload; i++) {
    p_pl[i] = p_savePl[i];
    p_ql[i] = p_saveQl[i];
    p_ip[i] = p_saveIp[i];
    p_iq[i] = p_saveIq[i];
    p_yp[i] = p_saveYp[i];
    p_yq[i] = p_saveYq[i];
  }
}

/**
 * Get available margin for generator
 * @param tag character ID for generator
 * @param current initial generation
 * @param pmin minimum allowable generation
 * @param pmax maximum allowable generation
 * @param status current status of generator
 */
void gridpack::powerflow::PFBus::getGeneratorMargins(
    std::vector<std::string> &tag, std::vector<double> &current,
    std::vector<double> &pmin, std::vector<double> &pmax,
    std::vector<int> &status)
{
  tag.clear();
  current.clear();
  pmin.clear();
  pmax.clear();
  status.clear();
  int i;
  for (i=0; i<p_ngen; i++) {
    tag.push_back(p_gid[i]);
    current.push_back(p_savePg[i]);
    pmin.push_back(p_pb[i]);
    pmax.push_back(p_pt[i]);
    status.push_back(p_gstatus[i]);
  }
}

/**
 * Get current value of loads
 * @param tag character ID for load
 * @param pl initial value of load real power
 * @param ql initial value of load reactive power
 * @param status current status of load
 */
void gridpack::powerflow::PFBus::getLoadPower(
    std::vector<std::string> &tag, std::vector<double> &pl,
    std::vector<double> &ql, std::vector<int> &status)
{
  tag.clear();
  pl.clear();
  ql.clear();
  status.clear();
  int i;
  for (i=0; i<p_nload; i++) {
    tag.push_back(p_lid[i]);
    pl.push_back(p_savePl[i]);
    ql.push_back(p_saveQl[i]);
    status.push_back(p_lstatus[i]);
  }
}

/**
 * Label bus as a source for real time path rating
 * @param flag identify bus as source
 */
void gridpack::powerflow::PFBus::setSource(bool flag)
{
  p_source = flag;
}

/**
 * Label bus as a sink for real time path rating
 * @param flag identify bus as sink
 */
void gridpack::powerflow::PFBus::setSink(bool flag)
{
  p_sink = flag;
}

/**
 * Store scale factor
 * @param scale factor for scaling generation or loads
 */
void gridpack::powerflow::PFBus::setScale(double scale)
{
  p_rtpr_scale = scale;
}

/**
 * Compute total voltage-dependent ZIP load power at given voltage
 * PSS/E convention:
 *   P_load = PL + IP*V + YP*V^2  (all terms add to real load)
 *   Q_load = QL + IQ*V - YQ*V^2  (IQ=inductive adds, YQ=capacitive subtracts)
 * This returns only the voltage-dependent parts:
 *   pzip = sum(IP*V + YP*V^2)
 *   qzip = sum(IQ*V - YQ*V^2)
 */
void gridpack::powerflow::PFBus::getZIPLoadPower(double V,
    double &pzip, double &qzip) const
{
  pzip = 0.0;
  qzip = 0.0;
  for (int i = 0; i < p_lstatus.size(); i++) {
    if (p_lstatus[i] == 1) {
      pzip += p_ip[i] * V + p_yp[i] * V * V;
      qzip += p_iq[i] * V - p_yq[i] * V * V;
    }
  }
}

/**
 *  Simple constructor
 */
gridpack::powerflow::PFBranch::PFBranch(void)
{
  p_reactance.clear();
  p_resistance.clear();
  p_tap_ratio.clear();
  p_phase_shift.clear();
  p_charging.clear();
  p_shunt_admt_g1.clear();
  p_shunt_admt_b1.clear();
  p_shunt_admt_g2.clear();
  p_shunt_admt_b2.clear();
  p_xform.clear();
  p_shunt.clear();
  p_elems = 0;
  p_theta = 0.0;
  p_sbase = 0.0;
  p_mode = YBus;
  // LTC defaults
  p_hasLTC = false;
  p_ltc_elem = -1;
  p_ltc_code = 0;
  p_ltc_cont = 0;
  p_ltc_cont_is_to = false;
  p_ltc_rma = 1.1;
  p_ltc_rmi = 0.9;
  p_ltc_vma = 1.1;
  p_ltc_vmi = 0.9;
  p_ltc_ntp = 33;
  p_ltc_step = 0.0;
  p_ltc_tap_init = 0.0;
  p_ltc_tap_prev = 0.0;
  p_ltc_locked = false;
  p_ltc_adj_count = 0;
}

/**
 *  Simple destructor
 */
gridpack::powerflow::PFBranch::~PFBranch(void)
{
}

/**
 * Return size of off-diagonal matrix block contributed by the component
 * for the forward/reverse directions
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::powerflow::PFBranch::matrixForwardSize(int *isize, int *jsize) const
{
  if (p_mode == Jacobian) {
    gridpack::powerflow::PFBus *bus1
      = dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
    gridpack::powerflow::PFBus *bus2
      = dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
    bool ok = !bus1->getReferenceBus();
    ok = ok && !bus2->getReferenceBus();
    ok = ok && !bus1->isIsolated();
    ok = ok && !bus2->isIsolated();
    ok = ok && (p_active);
    if (ok) {
#ifdef LARGE_MATRIX
      *isize = 2;
      *jsize = 2;
      return true;
#else
      // IREG PV buses have 2 vars/eqs (like PQ), not 1 (like standard PV)
      bool bus1PV = bus1->isPV() && !bus1->isIREG_PV();
      bool bus2PV = bus2->isPV() && !bus2->isIREG_PV();
      if (bus1PV && bus2PV) {
        *isize = 1;
        *jsize = 1;
        return true;
      } else if (bus1PV) {
        *isize = 1;
        *jsize = 2;
        return true;
      } else if (bus2PV) {
        *isize = 2;
        *jsize = 1;
        return true;
      } else {
        *isize = 2;
        *jsize = 2;
        return true;
      }
#endif
    } else {
      return false;
    }
  } else if (p_mode == YBus) {
    return YMBranch::matrixForwardSize(isize,jsize);
  }
  return false;
}
bool gridpack::powerflow::PFBranch::matrixReverseSize(int *isize, int *jsize) const
{
  if (p_mode == Jacobian) {
    gridpack::powerflow::PFBus *bus1
      = dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
    gridpack::powerflow::PFBus *bus2
      = dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
    bool ok = !bus1->getReferenceBus();
    ok = ok && !bus2->getReferenceBus();
    ok = ok && !bus1->isIsolated();
    ok = ok && !bus2->isIsolated();
    ok = ok && (p_active);
    if (ok) {
#ifdef LARGE_MATRIX
      *isize = 2;
      *jsize = 2;
      return true;
#else
      // IREG PV buses have 2 vars/eqs (like PQ), not 1 (like standard PV)
      bool bus1PV = bus1->isPV() && !bus1->isIREG_PV();
      bool bus2PV = bus2->isPV() && !bus2->isIREG_PV();
      if (bus1PV && bus2PV) {
        *isize = 1;
        *jsize = 1;
        return true;
      } else if (bus1PV) {
        *isize = 2;
        *jsize = 1;
        return true;
      } else if (bus2PV) {
        *isize = 1;
        *jsize = 2;
        return true;
      } else {
        *isize = 2;
        *jsize = 2;
        return true;
      }
#endif
    } else {
      return false;
    }
  } else if (p_mode == YBus) {
    return YMBranch::matrixReverseSize(isize,jsize);
  }
  return false;
}

/**
 * Return the values of the off-diagonal matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::powerflow::PFBranch::matrixForwardValues(ComplexType *values)
{
  if (p_mode == Jacobian) {
    double rvals[4];
    int nvals = forwardJacobianValues(rvals);
    for (int i=0; i<nvals; i++) values[i] = rvals[i];
    if (nvals == 0) {
      return false;
    } else {
      return true;
    }
  } else if (p_mode == YBus) {
    return YMBranch::matrixForwardValues(values);
  }
  return false;
}

bool gridpack::powerflow::PFBranch::matrixForwardValues(RealType *values)
{
  if (p_mode == Jacobian) {
    int nvals = forwardJacobianValues(values);
    if (nvals == 0) {
      return false;
    } else {
      return true;
    }
  }
  return false;
}

bool gridpack::powerflow::PFBranch::matrixReverseValues(ComplexType *values)
{
  if (p_mode == Jacobian) {
    double rvals[4];
    int nvals = reverseJacobianValues(rvals);
    for (int i=0; i<nvals; i++) values[i] = rvals[i];
    if (nvals == 0) {
      return false;
    } else {
      return true;
    }
  } else if (p_mode == YBus) {
    return YMBranch::matrixReverseValues(values);
  }
  return false;
}

bool gridpack::powerflow::PFBranch::matrixReverseValues(RealType *values)
{
  if (p_mode == Jacobian) {
    int nvals = reverseJacobianValues(values);
    if (nvals == 0) {
      return false;
    } else {
      return true;
    }
  }
  return false;
}

// Calculate contributions to the admittance matrix from the branches
void gridpack::powerflow::PFBranch::setYBus(void)
{
  YMBranch::setYBus();
  gridpack::ComplexType ret;
  ret = YMBranch::getForwardYBus();
  p_ybusr_frwd = real(ret);
  p_ybusi_frwd = imag(ret);
  ret = YMBranch::getReverseYBus();
  p_ybusr_rvrs = real(ret);
  p_ybusi_rvrs = imag(ret);
  // Not really a contribution to the admittance matrix but might as well
  // calculate phase angle difference between buses at each end of branch
  gridpack::powerflow::PFBus *bus1 =
    dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
  gridpack::powerflow::PFBus *bus2 =
    dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
  double pi = 4.0*atan(1.0);
  p_theta = (bus1->getPhase() - bus2->getPhase());
}

/**
 * Load values stored in DataCollection object into PFBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 */
void gridpack::powerflow::PFBranch::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  YMBranch::load(data);

  // This routine may be called more than once so clear all vectors
  p_reactance.clear();
  p_resistance.clear();
  p_tap_ratio.clear();
  p_phase_shift.clear();
  p_charging.clear();
  p_shunt_admt_g1.clear();
  p_shunt_admt_b1.clear();
  p_shunt_admt_g2.clear();
  p_shunt_admt_b2.clear();
  p_xform.clear();
  p_shunt.clear();
  p_rateA.clear();
  p_rateB.clear();
  p_rateC.clear();
  p_branch_status.clear();
  p_ckt.clear();
  p_ignore.clear();

  bool ok = true;
  data->getValue(BRANCH_NUM_ELEMENTS, &p_elems);
  double rvar;
  int ivar;
  double pi = 4.0*atan(1.0);
  p_active = false;
  ok = data->getValue(CASE_SBASE, &p_sbase);
  int idx;
  for (idx = 0; idx<p_elems; idx++) {
    bool xform = true;
    xform = xform && data->getValue(BRANCH_X, &rvar, idx);
    if (rvar <1.0e-5 && rvar >=0.0) rvar = 1.0e-5;
    if (rvar >-1.0e-5 && rvar <0.0) rvar =-1.0e-5;
    p_reactance.push_back(rvar);
    xform = xform && data->getValue(BRANCH_R, &rvar, idx);
    p_resistance.push_back(rvar);
    ok = data->getValue(BRANCH_SHIFT, &rvar, idx);
    rvar = -rvar*pi/180.0; 
    p_phase_shift.push_back(rvar);
    ok = data->getValue(BRANCH_TAP, &rvar, idx);
    p_tap_ratio.push_back(rvar); 
    if (rvar != 0.0) {
      p_xform.push_back(xform);
    } else {
      p_xform.push_back(false);
    }
    ivar = 1;
    ok = data->getValue(BRANCH_STATUS, &ivar, idx);
    p_branch_status.push_back(static_cast<bool>(ivar));
    if (ivar == 1 && ok) p_active = true;
    std::string tag;
    ok = data->getValue(BRANCH_CKT, &tag, idx);
    p_ckt.push_back(tag);
    bool shunt = true;
    shunt = shunt && data->getValue(BRANCH_B, &rvar, idx);
    p_charging.push_back(rvar);
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_G1, &rvar, idx);
    p_shunt_admt_g1.push_back(rvar);
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_B1, &rvar, idx);
    p_shunt_admt_b1.push_back(rvar);
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_G2, &rvar, idx);
    p_shunt_admt_g2.push_back(rvar);
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_B2, &rvar, idx);
    p_shunt_admt_b2.push_back(rvar);
    p_shunt.push_back(shunt);
    bool rate = true;
    rate = rate && data->getValue(BRANCH_RATING_A,&rvar,idx);
    p_rateA.push_back(rvar);
    rate = rate && data->getValue(BRANCH_RATING_B,&rvar,idx);
    p_rateB.push_back(rvar);
    rate = rate && data->getValue(BRANCH_RATING_C,&rvar,idx);
    p_rateC.push_back(rvar);
    p_ignore.push_back(false);
  }

  // Load LTC control data — scan elements for COD1=1 (voltage control)
  p_hasLTC = false;
  p_ltc_elem = -1;
  p_ltc_locked = false;
  p_ltc_adj_count = 0;
  for (int idx = 0; idx < p_elems; idx++) {
    int code1 = 0;
    int cont1 = 0;
    // v33+: TRANSFORMER_CODE1 and TRANSFORMER_CONT1
    // v23: TRANSFORMER_CONTROL stores controlled bus (nonzero implies voltage control)
    bool has_ltc_data = false;
    if (data->getValue(TRANSFORMER_CODE1, &code1, idx) && code1 == 1) {
      if (!data->getValue(TRANSFORMER_CONT1, &cont1, idx)) {
        data->getValue(TRANSFORMER_CONTROL, &cont1, idx);
      }
      has_ltc_data = (cont1 != 0);
    } else if (data->getValue(TRANSFORMER_CONTROL, &cont1, idx) && cont1 != 0) {
      // v23 format: TRANSFORMER_CONTROL present with nonzero bus implies voltage control
      code1 = 1;
      has_ltc_data = true;
    }
    if (has_ltc_data) {

      double rma = 1.1, rmi = 0.9, vma = 1.1, vmi = 0.9;
      int ntp = 33;
      data->getValue(TRANSFORMER_RMA, &rma, idx);
      data->getValue(TRANSFORMER_RMI, &rmi, idx);
      data->getValue(TRANSFORMER_VMA, &vma, idx);
      data->getValue(TRANSFORMER_VMI, &vmi, idx);
      data->getValue(TRANSFORMER_NTP, &ntp, idx);

      // Compute step size: use TRANSFORMER_STEP if available (v23), else from tap range
      double step = 0.0;
      if (!data->getValue(TRANSFORMER_STEP, &step, idx) || step < 1.0e-6) {
        if (ntp > 1) {
          step = (rma - rmi) / (ntp - 1);
        }
      }
      if (step < 1.0e-6) step = 0.00625;  // Default 0.625% step

      p_hasLTC = true;
      p_ltc_elem = idx;
      p_ltc_code = code1;
      p_ltc_cont = cont1;
      // Determine if controlled bus is the to-bus (tap direction reversal needed)
      int frombus = 0, tobus = 0;
      data->getValue(BRANCH_FROMBUS, &frombus);
      data->getValue(BRANCH_TOBUS, &tobus);
      p_ltc_cont_is_to = (cont1 == tobus);
      p_ltc_rma = rma;
      p_ltc_rmi = rmi;
      p_ltc_vma = vma;
      p_ltc_vmi = vmi;
      p_ltc_ntp = ntp;
      p_ltc_step = step;
      p_ltc_tap_init = p_tap_ratio[idx];
      p_ltc_tap_prev = p_tap_ratio[idx];
      break;  // Only support one LTC per branch
    }
  }
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
/**
 * Adjust LTC tap ratio to bring controlled voltage within deadband.
 * @param v_controlled voltage at controlled bus
 * @return true if a tap adjustment was made
 */
bool gridpack::powerflow::PFBranch::adjustLTC(double v_controlled)
{
  if (!p_hasLTC || p_ltc_locked) return false;
  int idx = p_ltc_elem;

  // Check if voltage is within deadband
  double vmi = p_ltc_vmi;
  double vma = p_ltc_vma;
  if (vma - vmi < 0.002) {
    double vmid = (vma + vmi) / 2.0;
    vmi = vmid - 0.001;
    vma = vmid + 0.001;
  }
  if (v_controlled >= vmi && v_controlled <= vma) {
    return false;  // No action needed
  }

  double tap_current = p_tap_ratio[idx];
  double tap_new = tap_current;

  // Tap direction depends on which bus is controlled:
  // - Tap on from-bus side: increasing tap INCREASES from-bus V, DECREASES to-bus V
  // - If controlling to-bus: low V → decrease tap; high V → increase tap
  // - If controlling from-bus: low V → increase tap; high V → decrease tap
  double step = p_ltc_step;
  if (p_ltc_cont_is_to) step = -step;  // Reverse direction for to-bus control

  if (v_controlled < p_ltc_vmi) {
    // Voltage too low
    tap_new = tap_current + step;
  } else {
    // Voltage too high
    tap_new = tap_current - step;
  }

  // Clamp to [RMI, RMA]
  if (tap_new > p_ltc_rma) tap_new = p_ltc_rma;
  if (tap_new < p_ltc_rmi) tap_new = p_ltc_rmi;

  // No change possible (already at limit)
  if (fabs(tap_new - tap_current) < 1.0e-8) {
    p_ltc_locked = true;
    return false;
  }

  // Cycle detection: if this adjustment would return to previous tap, lock
  if (p_ltc_adj_count > 0 && fabs(tap_new - p_ltc_tap_prev) < 1.0e-8) {
    p_ltc_locked = true;
    return false;
  }

  // Apply the tap change to both PFBranch and YMBranch tap ratio vectors
  p_ltc_tap_prev = tap_current;
  p_tap_ratio[idx] = tap_new;
  YMBranch::setParam(BRANCH_TAP, tap_new, idx);

  p_ltc_adj_count++;
  if (p_ltc_adj_count >= 10) {
    p_ltc_locked = true;
  }

  return true;
}

/**
 * Reset LTC to initial tap ratio
 */
void gridpack::powerflow::PFBranch::resetLTC()
{
  if (!p_hasLTC) return;
  int idx = p_ltc_elem;
  p_tap_ratio[idx] = p_ltc_tap_init;
  YMBranch::setParam(BRANCH_TAP, p_ltc_tap_init, idx);
  p_ltc_tap_prev = p_ltc_tap_init;
  p_ltc_locked = false;
  p_ltc_adj_count = 0;
}

/**
 * Check if branch has LTC control
 */
bool gridpack::powerflow::PFBranch::hasLTC() const
{
  return p_hasLTC;
}

/**
 * Get controlled bus number for LTC
 */
int gridpack::powerflow::PFBranch::getLTCControlledBus() const
{
  return p_ltc_cont;
}

/**
 * Get index of the LTC-controlled element within this branch
 */
int gridpack::powerflow::PFBranch::getLTCElementIndex() const
{
  return p_ltc_elem;
}

void gridpack::powerflow::PFBranch::setMode(int mode)
{
  if (mode == YBus) {
    YMBranch::setMode(gridpack::ymatrix::YBus);
  }
  p_mode = mode;
}

/**
 * Return the contribution to the Jacobian for the powerflow equations from
 * a branch
 * @param bus: pointer to the bus making the call
 * @param values: an array of 4 doubles that holds return metrix elements
 */
void gridpack::powerflow::PFBranch::getJacobian(gridpack::powerflow::PFBus *bus, double *values)
{
  double v;
  double cs, sn;
  double ybusr, ybusi;
  if (bus == getBus1().get()) {
    gridpack::powerflow::PFBus *bus2 =
      dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
    v = bus2->getVoltage();
    cs = cos(p_theta);
    sn = sin(p_theta);
    ybusr = p_ybusr_frwd;
    ybusi = p_ybusi_frwd;
  } else if (bus == getBus2().get()) {
    gridpack::powerflow::PFBus *bus1 =
      dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
    v = bus1->getVoltage();
    cs = cos(-p_theta);
    sn = sin(-p_theta);
    ybusr = p_ybusr_rvrs;
    ybusi = p_ybusi_rvrs;
  } else {
    // TODO: Some kind of error
    return;
  }
  values[0] = v*(ybusr*sn - ybusi*cs);
  values[1] = -v*(ybusr*cs + ybusi*sn);
  values[2] = (ybusr*cs + ybusi*sn);
  values[3] = (ybusr*sn - ybusi*cs);
}

/**
 * Return contribution to constraints
 * @param p: real part of constraint
 * @param q: imaginary part of constraint
 */
void gridpack::powerflow::PFBranch::getPQ(gridpack::powerflow::PFBus *bus, double *p, double *q)
{
  gridpack::powerflow::PFBus *bus1 = 
    dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
  double v1 = bus1->getVoltage();
  gridpack::powerflow::PFBus *bus2 =
    dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
  double v2 = bus2->getVoltage();
  double cs, sn;
  double ybusr, ybusi;
  p_theta = bus1->getPhase() - bus2->getPhase();
  if (bus == bus1) {
    cs = cos(p_theta);
    sn = sin(p_theta);
    ybusr = p_ybusr_frwd;
    ybusi = p_ybusi_frwd;
  } else if (bus == bus2) {
    cs = cos(-p_theta);
    sn = sin(-p_theta);
    ybusr = p_ybusr_rvrs;
    ybusi = p_ybusi_rvrs;
  } else {
    // TODO: Some kind of error
    return;
  }
  *p = v1*v2*(ybusr*cs+ybusi*sn);
  *q = v1*v2*(ybusr*sn-ybusi*cs);
}

/**
 * Return complex power for line element
 * @param tag describing line element on branch
 * @return complex power
 */
gridpack::ComplexType gridpack::powerflow::PFBranch::getComplexPower(
        std::string tag)
{
  gridpack::ComplexType vi, vj, Yii, Yij, s;
  s = ComplexType(0.0,0.0);
  gridpack::powerflow::PFBus *bus1 = 
    dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
  vi = bus1->getComplexVoltage();
  gridpack::powerflow::PFBus *bus2 =
    dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
  vj = bus2->getComplexVoltage();
  getLineElements(tag,&Yii,&Yij);
  s = vi*conj(Yii*vi+Yij*vj)*p_sbase;
  return s;
}

/**
 * Return complex power at the "to" end of line element
 * @param tag describing line element on branch
 * @return complex power at receiving end
 */
gridpack::ComplexType gridpack::powerflow::PFBranch::getReversePower(
        std::string tag)
{
  gridpack::ComplexType vi, vj, Yjj, Yji, s;
  s = ComplexType(0.0,0.0);
  gridpack::powerflow::PFBus *bus1 =
    dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
  vi = bus1->getComplexVoltage();
  gridpack::powerflow::PFBus *bus2 =
    dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
  vj = bus2->getComplexVoltage();
  getRvrsLineElements(tag,&Yjj,&Yji);
  s = vj*conj(Yjj*vj+Yji*vi)*p_sbase;
  return s;
}

/**
 * Write output from branches to standard out
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if branch is contributing string to output, false otherwise
 */
bool gridpack::powerflow::PFBranch::serialWrite(char *string, const int bufsize,
                                                const char *signal)
{
  char buf[128];
  gridpack::powerflow::PFBus *bus1
    = dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
  gridpack::powerflow::PFBus *bus2
    = dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
  bool ok = p_active;
  if (ok) {


  if (signal == NULL || !strcmp(signal,"flow_str")) {
    bool rating = false;
    if (signal != NULL) rating = !strcmp(signal,"flow_str");
    gridpack::ComplexType s;
    std::vector<std::string> tags = getLineTags();
    int i;
    int ilen = 0;
    for (i=0; i<p_elems; i++) {
      s = getComplexPower(tags[i]);
      double p = real(s);
      double q = imag(s);
      if (!p_branch_status[i]) p = 0.0;
      if (!p_branch_status[i]) q = 0.0;
      if (bus1->isIsolated() || bus2->isIsolated()) p=0.0;
      if (bus1->isIsolated() || bus2->isIsolated()) q=0.0;
      double perf = 0.0;
      int viol = 0;
      if (p_rateA[i] > 0.0) {
        perf = abs(s)/p_rateA[i];
        if (perf > 1.0) viol = 1;
        perf = perf*perf;
      }
      if (rating) {
        sprintf(buf, "%6d %6d %s %20.12e %20.12e %20.12e %20.12e %1d\n",
            getBus1OriginalIndex(),getBus2OriginalIndex(),tags[i].c_str(),
            p,q,perf,p_rateA[i],viol);
      } else {
        sprintf(buf, "     %6d      %6d     %s   %12.6f         %12.6f\n",
            getBus1OriginalIndex(),getBus2OriginalIndex(),tags[i].c_str(),
            p,q);
      }
      int len = strlen(buf);
      if (ilen + len < bufsize) {
        sprintf(string,"%s",buf);
        string += len;
        ilen += len;
      }
    }
    return true;
  } else if (!strcmp(signal,"fail_str")) {
    std::vector<std::string> tags = getLineTags();
    int i;
    int ilen = 0;
    for (i=0; i<p_elems; i++) {
        sprintf(buf, "     %6d      %6d     %s   %12.6f         %12.6f %12.6f %12.6f %1d\n",
            getBus1OriginalIndex(),getBus2OriginalIndex(),tags[i].c_str(),
            0.0,0.0,0.0,0.0,0);
      int len = strlen(buf);
      if (ilen + len < bufsize) {
        sprintf(string,"%s",buf);
        string += len;
        ilen += len;
      }
    }
    return true;
  } else if (!strcmp(signal,"flow")) {
    gridpack::ComplexType s;
    std::vector<std::string> tags = getLineTags();
    int i;
    bool found = false;
    int ilen = 0;
    for (i=0; i<p_elems; i++) {
      s = getComplexPower(tags[i]);
      double p = real(s);
      double q = imag(s);
      if (!p_branch_status[i]) p = 0.0;
      if (!p_branch_status[i]) q = 0.0;
      if (bus1->isIsolated() || bus2->isIsolated()) p=0.0;
      if (bus1->isIsolated() || bus2->isIsolated()) q=0.0;
      double S = sqrt(p*p+q*q);
      if (S > p_rateA[i] && p_rateA[i] != 0.0){
        sprintf(buf, "     %6d      %6d        %s  %12.6f         %12.6f     %8.2f     %8.2f%s\n",
    	  getBus1OriginalIndex(),getBus2OriginalIndex(),tags[i].c_str(),
          p,q,p_rateA[i],S/p_rateA[i]*100,"%");
        int len = strlen(buf);
        if (ilen + len < bufsize) {
          sprintf(string,"%s",buf);
          string += len;
          ilen += len;
        }
        found = true;
      }
    }
    return found;
  } else if (!strcmp(signal,"record")) {
    char *cptr = string;
    double pi = 4.0*atan(1.0);
    int slen = 0;
    int i, idx, jdx;
    for (i = 0; i<p_elems; i++) {
      idx = getBus1OriginalIndex();
      jdx = getBus2OriginalIndex();
      sprintf(buf,"%8d, %8d, %2s,",idx,jdx,p_ckt[i].c_str());
      int len = strlen(buf);
      if (len<=bufsize) {
        sprintf(cptr,"%s",buf);
        slen += len;
        cptr += len;
      }
//      double yi = 0.0;
//      double yj = 0.0;
      sprintf(buf," %12.6f, %12.6f, %12.6f, 0.0, 0.0, %8.2f, %8.2f, %8.2f,",
         p_resistance[i],p_reactance[i],p_charging[i],p_rateA[i],p_rateB[i],p_rateC[i]);
      len = strlen(buf);
      if (slen+len<=bufsize) {
        sprintf(cptr,"%s",buf);
        slen += len;
        cptr += len;
      }
      double rval = -180.0*p_phase_shift[i]/pi;
      idx = 0;
      if (p_branch_status[i]) idx = 1;
      sprintf(buf," %12.4f, %12.4f, %1d\n",p_tap_ratio[i],rval,idx);
      len = strlen(buf);
      if (slen+len<=bufsize) {
        sprintf(cptr,"%s",buf);
        slen += len;
        cptr += len;
      }
    }
    return true;
  }
  return false;
  } // ok loop
  return false;
}

/**
 * Get the status of the branch element
 * @param tag character string identifying branch element
 * @return status of branch element
 */
bool gridpack::powerflow::PFBranch::getBranchStatus(std::string tag)
{
  int i;
  int bsize = p_branch_status.size();
  for (i=0; i<bsize; i++) {
    if (tag == p_ckt[i]) {
      return p_branch_status[i];
    }
  }
  return false;
}

/**
 * Set the status of the branch element
 * @param tag character string identifying branch element
 * @param status status of branch element
 * @return true if circuit ID found, false otherwise
 */
bool gridpack::powerflow::PFBranch::setBranchStatus(std::string tag, bool status)
{
  int i;
  int bsize = p_branch_status.size();
  for (i=0; i<bsize; i++) {
    if (tag == p_ckt[i]) {
      p_branch_status[i] = status;
      YMBranch::setLineStatus(tag,status);
      return true;
    }
  }
  return false;
}

/**
 * get branch rating value
 * @param tag transmission element ID
 * @return branch rating value
 */
double gridpack::powerflow::PFBranch::getBranchRatingA(std::string tag)
{
  int i;
  int bsize = p_ckt.size();
  double ret = 0.0;
  for (i=0; i<bsize; i++) {
    if (tag == p_ckt[i]) {
      return p_rateA[i];
    }
  }
  return ret;
}

/**
 * get branch rating B value
 * @param tag transmission element ID
 * @return branch rating value
 */
double gridpack::powerflow::PFBranch::getBranchRatingB(std::string tag)
{
  int i;
  int bsize = p_ckt.size();
  double ret = 0.0;
  for (i=0; i<bsize; i++) {
    if (tag == p_ckt[i]) {
      return p_rateB[i];
    }
  }
  return ret;
}


/**
 * get branch rating C value
 * @param tag transmission element ID
 * @return branch rating value
 */
double gridpack::powerflow::PFBranch::getBranchRatingC(std::string tag)
{
  int i;
  int bsize = p_ckt.size();
  double ret = 0.0;
  for (i=0; i<bsize; i++) {
    if (tag == p_ckt[i]) {
      return p_rateC[i];
    }
  }
  return ret;
}

/**
 * Get list of line IDs
 * @return list of line identifiers
 */
std::vector<std::string> gridpack::powerflow::PFBranch::getLineIDs()
{
  return p_ckt;
}

/**
 * Set parameter to ignore voltage violations
 * @param tag identifier of line element
 * @param flag value of ignore parameter
 */
void gridpack::powerflow::PFBranch::setIgnore(std::string tag, bool flag)
{
  int i;
  int bsize = p_ckt.size();
  for (i=0; i<bsize; i++) {
    if (tag == p_ckt[i]) {
      p_ignore[i] = flag;
    }
  }
}

/**
 * Get parameter to ignore voltage violations
 * @param tag identifier of line element
 * @return value of ignore parameter
 */
bool gridpack::powerflow::PFBranch::getIgnore(std::string tag)
{
  int i;
  int bsize = p_ckt.size();
  bool ret = false;
  for (i=0; i<bsize; i++) {
    if (tag == p_ckt[i]) {
      return p_ignore[i];
    }
  }
  return ret;
}

/**
 * Evaluate off-diagonal block of Jacobian for power flow calculation
 * and return result as an array of real values
 * @param rvals values of Jacobian block
 * @return number of values returned
 */
int gridpack::powerflow::PFBranch::forwardJacobianValues(double *rvals)
{
  gridpack::powerflow::PFBus *bus1
    = dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
  gridpack::powerflow::PFBus *bus2
    = dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
  bool ok = !bus1->getReferenceBus();
  ok = ok && !bus2->getReferenceBus();
  ok = ok && !bus1->isIsolated();
  ok = ok && !bus2->isIsolated();
  ok = ok && (p_active);
  int nvals;
  if (ok) {
    double t11, t12, t21, t22;
    double cs = cos(p_theta);
    double sn = sin(p_theta);
    // IREG PV buses have 2 vars/eqs (like PQ), not 1 (like standard PV)
    bool bus1PV = bus1->isPV() && !bus1->isIREG_PV();
    bool bus2PV = bus2->isPV() && !bus2->isIREG_PV();
#ifdef LARGE_MATRIX
    rvals[0] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
    rvals[1] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
    rvals[2] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
    rvals[3] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
    rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
    rvals[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
    rvals[2] *= bus1->getVoltage();
    rvals[3] *= bus1->getVoltage();
    // fix up matrix if one or both buses at the end of the branch is a PV bus
    if (bus1PV && bus2PV) {
      rvals[1] = 0.0;
      rvals[2] = 0.0;
      rvals[3] = 0.0;
    } else if (bus1PV) {
      rvals[1] = 0.0;
      rvals[3] = 0.0;
    } else if (bus2PV) {
      rvals[2] = 0.0;
      rvals[3] = 0.0;
    }
    nvals = 4;
#else
    if (bus1PV && bus2PV) {
      rvals[0] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      nvals = 1;
    } else if (bus1PV) {
      rvals[0] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
      rvals[1] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[1] *= bus1->getVoltage();
      nvals = 2;
    } else if (bus2PV) {
      rvals[0] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
      rvals[1] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
      nvals = 2;
    } else {
      rvals[0] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
      rvals[1] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
      rvals[2] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
      rvals[3] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[2] *= bus1->getVoltage();
      rvals[3] *= bus1->getVoltage();
      nvals = 4;
    }
#endif
    // For IREG PV bus1: replace Q-row with V_remote constraint
    if (bus1->isIREG_PV()) {
      int remote = bus1->getIREGRemoteBus();
      int bus2_idx = bus2->getOriginalIndex();
      if (nvals >= 2) rvals[1] = 0.0;  // ∂(V_remote)/∂θ_bus2 = 0
      if (nvals >= 4) {
        rvals[3] = (bus2_idx == remote) ? 1.0 : 0.0;
      }
    }
    return nvals;
  } else {
    return 0;
  }
}

int gridpack::powerflow::PFBranch::reverseJacobianValues(double *rvals)
{
  gridpack::powerflow::PFBus *bus1
    = dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
  gridpack::powerflow::PFBus *bus2
    = dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
  bool ok = !bus1->getReferenceBus();
  ok = ok && !bus2->getReferenceBus();
  ok = ok && !bus1->isIsolated();
  ok = ok && !bus2->isIsolated();
  ok = ok && (p_active);
  int nvals;
  if (ok) {
    double t11, t12, t21, t22;
    double cs = cos(-p_theta);
    double sn = sin(-p_theta);
    // IREG PV buses have 2 vars/eqs (like PQ), not 1 (like standard PV)
    bool bus1PV = bus1->isPV() && !bus1->isIREG_PV();
    bool bus2PV = bus2->isPV() && !bus2->isIREG_PV();
#ifdef LARGE_MATRIX
    rvals[0] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
    rvals[1] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
    rvals[2] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
    rvals[3] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
    rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
    rvals[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
    rvals[2] *= bus2->getVoltage();
    rvals[3] *= bus2->getVoltage();
    // fix up matrix if one or both buses at the end of the branch is a PV bus
    if (bus1PV && bus2PV) {
      rvals[1] = 0.0;
      rvals[2] = 0.0;
      rvals[3] = 0.0;
    } else if (bus1PV) {
      rvals[2] = 0.0;
      rvals[3] = 0.0;
    } else if (bus2PV) {
      rvals[1] = 0.0;
      rvals[3] = 0.0;
    }
    nvals = 4;
#else
    if (bus1PV && bus2PV) {
      rvals[0] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      nvals = 1;
    } else if (bus1PV) {
      rvals[0] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
      rvals[1] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
      nvals = 2;
    } else if (bus2PV) {
      rvals[0] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
      rvals[1] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[1] *= bus2->getVoltage();
      nvals = 2;
    } else {
      rvals[0] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
      rvals[1] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
      rvals[2] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
      rvals[3] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[2] *= bus2->getVoltage();
      rvals[3] *= bus2->getVoltage();
      nvals = 4;
    } 
#endif
    return nvals;
  } else {
    return 0;
  }
}
