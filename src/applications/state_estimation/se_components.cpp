/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_components.cpp
 * @author Yousu Chen
 * @date   2/24/2014 
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include <iostream>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "se_components.hpp"

//#define LARGE_MATRIX

/**
 *  Simple constructor
 */
gridpack::state_estimation::SEBus::SEBus(void)
{
  p_shunt_gs = 0.0;
  p_shunt_bs = 0.0;
  p_v = 0.0;
  p_a = 0.0;
  p_theta = 0.0;
  p_angle = 0.0;
  p_voltage = 0.0;
  p_pl = 0.0;
  p_ql = 0.0;
  p_sbase = 0.0;
  p_mode = YBus;
  p_rowJidx.clear();
  p_rowRidx.clear();
  p_colJidx.clear();
  p_colRidx.clear();
  setReferenceBus(false);
}

/**
 *  Simple destructor
 */
gridpack::state_estimation::SEBus::~SEBus(void)
{
}

/**
 *  Return size of matrix block contributed by the component
 *  @param isize, jsize: number of rows and columns of matrix block
 *  @return: false if network component does not contribute matrix element
 */
bool gridpack::state_estimation::SEBus::matrixDiagSize(int *isize, int *jsize) const
{
  if (p_mode == YBus) {
    return YMBus::matrixDiagSize(isize,jsize);
  } else if (p_mode == Jacobian_H) {
    *isize = 1;
    *jsize = 1;
    return true;
  }
  return true;
}

/**
 * Return the values of the matrix block. The values are
 * returned in row-major order.
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::state_estimation::SEBus::matrixDiagValues(ComplexType *values)
{
  if (p_mode == YBus) {
    return YMBus::matrixDiagValues(values);
/*  } else if (p_mode == Jacobian_H) {
    std::vector<gridpack::state_estimation::Measurement>
    meas = p_meas; //p_meas supposed to be all measurements on this bus
    int nmeas = meas.size();
    int i;
    for (i=0; i<nmeas; i++ ) {
       if (meas[i].p_type == "VM") {
          values[0] = 0.0; 
          values[1] = 1.0; 
       } else if (meas[i].p_type == "PI") {
         std::vector<boost::shared_ptr<BaseComponent> > branches;
         getNeighborBranches(branches);
         int size = branches.size();
         int j;
         double ret1 = 0.0;
         double ret2 = 0.0;
         for (j=0; j<size; j++) {
           gridpack::state_estimation::SEBranch *branch
             = dynamic_cast<gridpack::state_estimation::SEBranch*>(branches[i].get());
          branch->getVTheta(this,&v,&theta);
          ret1 += p_v * v * (p_ybusr_frwd*sin(theta) + p_ybusi_frwd*cos(theta)) - p_v * p_v * p_ybusi;
          ret2 +=  v * (p_ybusr_frwd*cos(theta) + p_ybusi_frwd*sin(theta)) + p_v * p_ybusr;
         }
          values[0] = ret1;
          values[1] = ret2;
       //} // to add other bus measurements
       }
*/
  }
}

/**
 * Return the size of the block that this component contributes to the
 * vector
 * @param size: size of vector block
 * @return: false if component does not contribute to vector
 */
bool gridpack::state_estimation::SEBus::vectorSize(int *size) const
{
  return true;
}

/**
 * Return the values of the vector block
 * @param values: pointer to vector values
 * @return: false if network component does not contribute
 *        vector element
 */
bool gridpack::state_estimation::SEBus::vectorValues(ComplexType *values)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto buses 
 * @param values array containing voltage magnitude and angle
 */
void gridpack::state_estimation::SEBus::setValues(gridpack::ComplexType *values)
{
  double vt = p_v;
  double at = p_a;
  p_a -= real(values[0]);
#ifdef LARGE_MATRIX
  p_v -= real(values[1]);
#else
  if (!p_isPV) {
    p_v -= real(values[1]);
  }
#endif
  *p_vAng_ptr = p_a;
  *p_vMag_ptr = p_v;
//  printf("at: %12.6f vt: %12.6f da: %12.6f dv: %12.6f  p_a: %12.6f p_v: %12.6f\n",
//      at,vt,real(values[0]),real(values[1]),p_a,p_v);
}

/**
 * Return the size of the buffer used in data exchanges on the network.
 * For this problem, the voltage magnitude and phase angle need to be exchanged
 * @return size of buffer
 */
int gridpack::state_estimation::SEBus::getXCBufSize(void)
{
  return 2*sizeof(double);
}

/**
 * Assign pointers for voltage magnitude and phase angle
 */
void gridpack::state_estimation::SEBus::setXCBuf(void *buf)
{
  p_vAng_ptr = static_cast<double*>(buf);
  p_vMag_ptr = p_vAng_ptr+1;
  // Note: we are assuming that the load function has been called BEFORE
  // the factory setExchange method, so p_a and p_v are set with their initial
  // values.
  *p_vAng_ptr = p_a;
  *p_vMag_ptr = p_v;
}

/**
 * Load values stored in DataCollection object into SEBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::state_estimation::SEBus::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  YMBus::load(data);

  bool ok = data->getValue(CASE_SBASE, &p_sbase);
  data->getValue(BUS_VOLTAGE_ANG, &p_angle);
  data->getValue(BUS_VOLTAGE_MAG, &p_voltage); 
  p_v = p_voltage;
  double pi = 4.0*atan(1.0);
  p_angle = p_angle*pi/180.0;
  p_a = p_angle;
  int itype;
  data->getValue(BUS_TYPE, &itype);
  if (itype == 3) {
    setReferenceBus(true);
  }

  // if BUS_TYPE = 2 then bus is a PV bus
  p_isPV = false;
  // if (itype == 2) p_isPV = true;

  // added p_pg,p_qg,p_pl,p_ql,p_sbase;
  p_load = true;
  p_load = p_load && data->getValue(LOAD_PL, &p_pl);
  p_load = p_load && data->getValue(LOAD_QL, &p_ql);
  //printf("p_pl=%f,p_ql=%f\n",p_pl,p_ql);
  bool lgen;
  int i, ngen, gstatus;
  double pg, qg, vs;
  ngen = 0;
  if (data->getValue(GENERATOR_NUMBER, &ngen)) {
    for (i=0; i<ngen; i++) {
      lgen = true;
      lgen = lgen && data->getValue(GENERATOR_PG, &pg,i);
      lgen = lgen && data->getValue(GENERATOR_QG, &qg,i);
      lgen = lgen && data->getValue(GENERATOR_VS, &vs,i);
      lgen = lgen && data->getValue(GENERATOR_STAT, &gstatus,i);
      if (lgen) {
        p_pg.push_back(pg);
        p_qg.push_back(qg);
        p_gstatus.push_back(gstatus);
        if (gstatus == 1) {
          p_v = vs; //reset initial PV voltage to set voltage
          if (itype == 2) p_isPV = true;
        }
      }
    }
  }

}

/**
 * Set values of YBus matrix. These can then be used in subsequent
 * calculations
 */
void gridpack::state_estimation::SEBus::setYBus(void)
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
gridpack::ComplexType gridpack::state_estimation::SEBus::getYBus(void)
{
  return YMBus::getYBus();
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::state_estimation::SEBus::setMode(int mode)
{
  if (mode == YBus) {
    YMBus::setMode(gridpack::ymatrix::YBus);
  }
  p_mode = mode;
}

/**
 * Return the value of the voltage magnitude on this bus
 * @return: voltage magnitude
 */
double gridpack::state_estimation::SEBus::getVoltage()
{
  return *p_vMag_ptr;
}

/**
 * Return whether or not the bus is a PV bus (V held fixed in powerflow
 * equations)
 * @return true if bus is PV bus
 */
bool gridpack::state_estimation::SEBus::isPV(void)
{
  return p_isPV;
}

/**
 * Return whether or not a bus is isolated
 * @return true if bus is isolated
 */
bool gridpack::state_estimation::SEBus::isIsolated(void) const
{
  return YMBus::isIsolated();
}

/**
 * Return the value of the phase angle on this bus
 * @return: phase angle
 */
double gridpack::state_estimation::SEBus::getPhase()
{
  return *p_vAng_ptr;
}

/**
 * Write output from buses to standard out
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool gridpack::state_estimation::SEBus::serialWrite(char *string,
    const int bufsize, const char *signal)
{
  if (signal == NULL) {
    double pi = 4.0*atan(1.0);
    double angle = p_a*180.0/pi;
    sprintf(string, "     %6d      %12.6f         %12.6f\n",
        getOriginalIndex(),angle,p_v);
  } else if (!strcmp(signal,"pq")) {
    gridpack::ComplexType v[2];
    vectorValues(v);
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    getNeighborBranches(branches);
    sprintf(string, "     %6d      %12.6f         %12.6f      %2d\n",
        getOriginalIndex(),real(v[0]),real(v[1]),branches.size());
  }
  return true;
}

/**
 * Add a measurement to the bus
 * @param measurement a measurement struct that will be used to
 * assign
 * internal paramters
 */
void gridpack::state_estimation::SEBus::addMeasurement(
    gridpack::state_estimation::Measurement measurement)
{
  p_meas.push_back(measurement);
  //TODO: Implement this method
}

/**
 * Return the complex voltage on this bus
 * @return the complex voltage
 */
gridpack::ComplexType gridpack::state_estimation::SEBus::getComplexVoltage(void)
{
  gridpack::ComplexType ret(cos(p_a),sin(p_a));
  ret = ret*p_v;
  return ret;
}

/**
 * Configure buses with state estimation parameters. These can be
 * used in other methods
 */
void gridpack::state_estimation::SEBus::configureSE(void)
{
  // Calculate the number of matrix values associated with this bus
  int nmeas = p_meas.size(); // Suppose p_meas is the vector of all the measurements on this bus
  int ncnt = 0;
  int i, j, nsize;
  int busid = getOriginalIndex();
  for (i=0; i<nmeas; i++) {
    std::string type = p_meas[i].p_type;
    if (type == "VM" || type == "VA") {
      if (!getReferenceBus()) { 
        ncnt += 2;
      } else {
        ncnt++;
      }
    } else if (type == "PI" || type == "QI") {
      std::vector<boost::shared_ptr<BaseComponent> > branch_nghbrs;
      getNeighborBranches(branch_nghbrs);
      nsize = branch_nghbrs.size();
      double ret1, ret2;
      for (j=0; j<nsize; j++) {
        SEBranch *branch
          = dynamic_cast<SEBranch*>(branch_nghbrs[j].get());
        SEBus *bus = dynamic_cast<SEBus*>(branch->getBus1().get());
        if (bus == this) bus = dynamic_cast<SEBus*>(branch->getBus2().get());
        if (!bus->getReferenceBus()) {
          ncnt += 2;
        } else {
          ncnt++;
        }
      }
      if (!getReferenceBus()) { 
        ncnt += 2;
      } else {
        ncnt++;
      }
    }
  } 
  p_numElements = ncnt;
}

/**
 * Return number of rows in matrix from component
 * @return number of rows from component
 */
int gridpack::state_estimation::SEBus::matrixNumRows() const
{
  return p_meas.size();
}

/**
 * Return number of cols in matrix from component
 * @return number of cols from component
 */
int gridpack::state_estimation::SEBus::matrixNumCols() const
{
  if (p_mode == Jacobian_H) {
    // Check to see if this bus has measurements or is attached to anything that
    // has measurements
    bool meas = false;
    if (p_meas.size() > 0) meas = true;
    if (!meas) {
      std::vector<boost::shared_ptr<BaseComponent> > branch_nghbrs;
      getNeighborBranches(branch_nghbrs);
      std::vector<boost::shared_ptr<BaseComponent> > bus_nghbrs;
      getNeighborBuses(bus_nghbrs);
      int nsize = branch_nghbrs.size();
      int i;
      gridpack::state_estimation::SEBus *bus;
      gridpack::state_estimation::SEBranch *branch;
      for (i=0; i<nsize && !meas; i++) {
        bus = dynamic_cast<SEBus*>(bus_nghbrs[i].get());
        branch = dynamic_cast<SEBranch*>(branch_nghbrs[i].get());
        if (bus->matrixNumRows() > 0) meas = true;
        if (branch->matrixNumRows() > 0) meas = true;
      }
    }
    if (!meas) return 0;
    // Bus has measurements associated with it.
    if (!getReferenceBus()) {
      return 2;
    } else {
      return 1;
    }
  } else if (p_mode == R_inv) {
    return p_meas.size();
  }
}

/**
 * Set row indices corresponding to the rows contributed by this component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @param idx matrix index of row irow
 */
void gridpack::state_estimation::SEBus::matrixSetRowIndex(int irow, int idx)
{
  if (p_mode == Jacobian_H) {
    if (irow < p_rowJidx.size()) {
      p_rowJidx[irow] = idx;
    } else {
      p_rowJidx.push_back(idx);
    }
  } else if (p_mode == R_inv) {
    if (irow < p_rowRidx.size()) {
      p_rowRidx[irow] = idx;
    } else {
      p_rowRidx.push_back(idx);
    }
  }
}

/**
 * Set column indices corresponding to the columns contributed by this component
 * @param icol index of column contributed by this component (e.g. if component
 * contributes 3 columns then icol is between 0 and 2)
 * @param idx matrix index of column icol
 */
void gridpack::state_estimation::SEBus::matrixSetColIndex(int icol, int idx)
{
  if (p_mode == Jacobian_H) {
    if (icol < p_colJidx.size()) {
      p_colJidx[icol] = idx;
    } else {
      p_colJidx.push_back(idx);
    }
  } else if (p_mode == R_inv) {
    if (icol < p_colRidx.size()) {
      p_colRidx[icol] = idx;
    } else {
      p_colRidx.push_back(idx);
    }
  }
}

/**
 * Get the row index corresponding to the rows contributed by this component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @return matrix index of row irow
 */
int gridpack::state_estimation::SEBus::matrixGetRowIndex(int idx)
{
  if (p_mode == Jacobian_H) {
  return p_rowJidx[idx];
  } else if (p_mode == R_inv) {
  return p_rowRidx[idx];
  }
}

/**
 * Get the column index corresponding to the columns contributed by this component
 * @param icol index of column contributed by this component (e.g. if component
 * contributes 3 columns then icol is between 0 and 2)
 * @return matrix index of column icol
 */
int gridpack::state_estimation::SEBus::matrixGetColIndex(int idx)
{
  if (p_mode == Jacobian_H) {
    if (idx >= p_colJidx.size())
      printf("violation in bus:matrixGetColIndex bus: %d size: %d idx: %d\n",
          getOriginalIndex(),idx,p_colJidx.size());
    return p_colJidx[idx];
  } else if (p_mode == R_inv) {
    return p_colRidx[idx];
  }
}

/**
 * Return the number of matrix values contributed by this component
 * @return number of matrix values
 */
int gridpack::state_estimation::SEBus::matrixNumValues() const
{
  if (p_mode == Jacobian_H) {
    return p_numElements;
  } else if (p_mode == R_inv) {
    return p_meas.size();
  }
}

/**
 * Return values from a matrix block
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
*/
void gridpack::state_estimation::SEBus::matrixGetValues(ComplexType *values, int *rows, int *cols)
{
  if (p_mode == Jacobian_H) {
    int nmeas = p_meas.size(); // Suppose p_meas is the vector of all the measurements on this bus
    int ncnt = 0;
    int i, j, im, jm, nsize;
    double v, theta, yfbusr,yfbusi;
    std::string ctk, type;
    for (i=0; i<nmeas; i++) {
      im = matrixGetRowIndex(i);
      ctk = p_meas[i].p_ckt;
      type = p_meas[i].p_type;
      if (type == "VM") {
#if 0
        std::vector<boost::shared_ptr<BaseComponent> > branch_nghbrs;
        getNeighborBranches(branch_nghbrs);
        nsize = branch_nghbrs.size();
        for (j=0; j<nsize; j++) {
          jm = matrixGetColIndex(0);
          values[ncnt] = gridpack::ComplexType(0.0,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
          jm = matrixGetColIndex(1);
          values[ncnt] = gridpack::ComplexType(0.0,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
        }
#endif
        if (!getReferenceBus()) { 
          jm = matrixGetColIndex(0);
          values[ncnt] = gridpack::ComplexType(0.0,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
          jm = matrixGetColIndex(1);
          values[ncnt] = gridpack::ComplexType(1.0,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
        } else {
          jm = matrixGetColIndex(0);
          values[ncnt] = gridpack::ComplexType(1.0,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
        }
      } else if (type == "PI") {
        std::vector<boost::shared_ptr<BaseComponent> > branch_nghbrs;
        getNeighborBranches(branch_nghbrs);
        nsize = branch_nghbrs.size();
        double ret1, ret2;
        for (j=0; j<nsize; j++) {
          SEBranch *branch
            = dynamic_cast<SEBranch*>(branch_nghbrs[j].get());
          SEBus *bus = dynamic_cast<SEBus*>(branch->getBus1().get());
          if (bus == this) bus = dynamic_cast<SEBus*>(branch->getBus2().get());
          branch->getVTheta(this, &v, &theta);
          ComplexType yfbus=branch->getForwardYBus();
          yfbusr = real (yfbus);
          yfbusi = imag (yfbus);
          // to discuss, how to use YBus branch data in bus 
          ret1 += p_v * v * (-yfbusr*sin(theta) + yfbusi*cos(theta));
          if (!bus->getReferenceBus()) {
            values[ncnt] = gridpack::ComplexType(p_v*v*(yfbusr*sin(theta)-yfbusi*cos(theta)),0.0);
            jm = bus->matrixGetColIndex(0);
            rows[ncnt] = im;
            cols[ncnt] = jm;
            ncnt++;
            jm = bus->matrixGetColIndex(1);
          } else {
            jm = bus->matrixGetColIndex(0);
          }
          ret2 +=  v * (yfbusr*cos(theta) + yfbusi*sin(theta));
          values[ncnt] = gridpack::ComplexType(p_v*(yfbusr*cos(theta)+yfbusi*sin(theta)),0.0);
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
        }
        if (!getReferenceBus()) {
          ret1 -= p_v * p_v * p_ybusi;
          jm = matrixGetColIndex(0);
          values[ncnt] = gridpack::ComplexType(ret1,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
          jm = matrixGetColIndex(1);
        } else {
          jm = matrixGetColIndex(0);
        }
        ret2 += p_v * p_ybusr;
        values[ncnt] = gridpack::ComplexType(ret2,0.0); 
        rows[ncnt] = im;
        cols[ncnt] = jm;
        ncnt++;
      } else if (type == "QI") {
        std::vector<boost::shared_ptr<BaseComponent> > branch_nghbrs;
        getNeighborBranches(branch_nghbrs);
        nsize = branch_nghbrs.size();
        double ret1, ret2;
        for (j=0; j<nsize; j++) {
          SEBranch *branch
            = dynamic_cast<SEBranch*>(branch_nghbrs[j].get());
          SEBus *bus = dynamic_cast<SEBus*>(branch->getBus1().get());
          if (bus == this) bus = dynamic_cast<SEBus*>(branch->getBus2().get());
          branch->getVTheta(this, &v, &theta);
          ComplexType yfbus=branch->getForwardYBus();
          yfbusr = real (yfbus);
          yfbusi = imag (yfbus);
          ret1 += p_v * v * (yfbusr*cos(theta) + yfbusi*sin(theta));
          if (!bus->getReferenceBus()) {
            values[ncnt] = gridpack::ComplexType(p_v*v*(-yfbusr*cos(theta)-yfbusi*sin(theta)),0.0);
            jm = bus->matrixGetColIndex(0);
            rows[ncnt] = im;
            cols[ncnt] = jm;
            ncnt++;
            jm = bus->matrixGetColIndex(1);
          } else {
            jm = bus->matrixGetColIndex(0);
          }
          ret2 +=  v * (yfbusr*cos(theta) + yfbusi*sin(theta));
          values[ncnt] = gridpack::ComplexType(p_v*(yfbusr*sin(theta)-yfbusi*cos(theta)),0.0);
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
        }
        ret1 -= p_v * p_v * p_ybusr;
        if (!getReferenceBus()) {
          jm = matrixGetColIndex(0);
          values[ncnt] = gridpack::ComplexType(ret1,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
          jm = matrixGetColIndex(1);
        } else {
          jm = matrixGetColIndex(0);
        }
        ret2 -= p_v * p_ybusi;
        values[ncnt] = gridpack::ComplexType(ret2,0.0); 
        rows[ncnt] = im;
        cols[ncnt] = jm;
        ncnt++;
      } else if (type == "VA") {
#if 0
        std::vector<boost::shared_ptr<BaseComponent> > branch_nghbrs;
        getNeighborBranches(branch_nghbrs);
        nsize = branch_nghbrs.size();
        for (j=0; j<nsize; j++) {
          jm = matrixGetColIndex(0);
          values[ncnt] = gridpack::ComplexType(0.0,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
          jm = matrixGetColIndex(1);
          values[ncnt] = gridpack::ComplexType(0.0,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
        }
#endif
        if (!getReferenceBus()) {
          jm = matrixGetColIndex(0);
          values[ncnt] = gridpack::ComplexType(1.0,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
          jm = matrixGetColIndex(1);
        } else {
          jm = matrixGetColIndex(0);
        }
        values[ncnt] = gridpack::ComplexType(0.0,0.0); 
        rows[ncnt] = im;
        cols[ncnt] = jm;
        ncnt++;
      }
    } 
  } else if (p_mode == R_inv) {
    int nsize = p_meas.size();
    int i;
    for (i=0; i<nsize; i++) {
      if (p_meas[i].p_deviation != 0.0) {
        values[i] = 1.0/(p_meas[i].p_deviation*p_meas[i].p_deviation);
      } else {
        values[i] = 0.0;
      }
      rows[i] = matrixGetRowIndex(i);
      cols[i] = matrixGetColIndex(i);
    }
  }
}

/**
 * Return values from a vector
 * @param values: pointer to vector values (z-h(x))
 * @param idx: pointer to vector index 
*/
void gridpack::state_estimation::SEBus:: VectorGetElementValues(ComplexType *values, int *idx)
{
  if (p_mode == Jacobian_H) {
    int nmeas = p_meas.size(); // Suppose p_meas is the vector of all the measurements on this bus
    int ncnt = 0;
    int i, j, im, jm, nsize;
    double v, theta,yfbusr,yfbusi;
    for (i=0; i<nmeas; i++) {
       vectorGetElementIndices(idx);
       if (p_meas[i].p_type == "VM") {
         int index = getGlobalIndex();
//         values[ncnt] = static_cast<double>(index),meas[i].p_value-p_v; 
//         values[ncnt] = p_meas[i].p_value-p_v; 
         values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-p_v),0.0);
         ncnt++;
       } else if (p_meas[i].p_type == "PI") {
         std::vector<boost::shared_ptr<BaseComponent> > branch_nghbrs;
         getNeighborBranches(branch_nghbrs);
         nsize = branch_nghbrs.size();
         gridpack::state_estimation::SEBranch *branch
           = dynamic_cast<gridpack::state_estimation::SEBranch*>(branch_nghbrs[i].get());
         ComplexType yfbus=branch->getForwardYBus();
         yfbusr = real (yfbus);
         yfbusi = imag (yfbus);
         double ret;
         for (j=0; j<nsize; j++) {
           branch->getVTheta(this, &v, &theta);
           ret +=  v * (yfbusr*cos(theta) + yfbusi*sin(theta));
           ncnt++;
         }
         ret *= p_v; 
         int index = getGlobalIndex();
         vectorGetElementIndices(idx);
//         values[ncnt] = static_cast<double>(index), meas[i].p_value-ret;
//         values[ncnt] = p_meas[i].p_value-ret;
         values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-ret),0.0);
         ncnt++;
       } else if (p_meas[i].p_type == "QI") {
         std::vector<boost::shared_ptr<BaseComponent> > branch_nghbrs;
         getNeighborBranches(branch_nghbrs);
         nsize = branch_nghbrs.size();
         gridpack::state_estimation::SEBranch *branch
           = dynamic_cast<gridpack::state_estimation::SEBranch*>(branch_nghbrs[i].get());
         ComplexType yfbus=branch->getForwardYBus();
         yfbusr = real (yfbus);
         yfbusi = imag (yfbus);
         double ret;
         for (j=0; j<nsize; j++) {
           branch->getVTheta(this,&v,&theta);
           ret +=  v * (yfbusr*sin(theta) - yfbusi*cos(theta));
           ncnt++;
         }
         ret *= p_v; 
         int index = getGlobalIndex();
         vectorGetElementIndices(idx);
//         values[ncnt] = static_cast<double>(index), p_meas[i].p_value-ret;
//         values[ncnt] = p_meas[i].p_value-ret;
         values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-ret),0.0);
         ncnt++;
      } else if (p_meas[i].p_type == "VA") {
         int index = getGlobalIndex();
//         values[ncnt] = static_cast<double>(index),p_meas[i].p_value-p_a; 
//         values[ncnt] = p_meas[i].p_value-p_a; 
         values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-p_a),0.0);
         ncnt++;
      }
    } 
  } else if (p_mode == R_inv) {
  }
}

/**
 *  Simple constructor
 */
gridpack::state_estimation::SEBranch::SEBranch(void)
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
  p_branch_status.clear();
  p_elems = 0;
  p_theta = 0.0;
  p_sbase = 0.0;
  p_rowJidx.clear();
  p_rowRidx.clear();
  p_colJidx.clear();
  p_colRidx.clear();
  p_mode = YBus;
}

/**
 *  Simple destructor
 */
gridpack::state_estimation::SEBranch::~SEBranch(void)
{
}

/**
 * Return size of off-diagonal matrix block contributed by the component
 * for the forward/reverse directions
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::state_estimation::SEBranch::matrixForwardSize(int *isize, int *jsize) const
{
  if (p_mode == YBus) {
    return YMBranch::matrixForwardSize(isize,jsize);
  }
}
bool gridpack::state_estimation::SEBranch::matrixReverseSize(int *isize, int *jsize) const
{
  if (p_mode == YBus) {
    return YMBranch::matrixReverseSize(isize,jsize);
  }
}

/**
 * Return the values of the off-diagonal matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::state_estimation::SEBranch::matrixForwardValues(ComplexType *values)
{
  if (p_mode == YBus) {
//    values[0] = p_ybusr_frwd;
//    values[1] = p_ybusi_frwd;
//    values[2] = -p_ybusi_frwd;
//    values[3] = p_ybusr_frwd;
    return YMBranch::matrixForwardValues(values);
  }
}

bool gridpack::state_estimation::SEBranch::matrixReverseValues(ComplexType *values)
{
  if (p_mode == YBus) {
  //  values[0] = p_ybusr_rvrs;
  //  values[1] = p_ybusi_rvrs;
  //  values[2] = -p_ybusi_rvrs;
  //  values[3] = p_ybusr_rvrs;
    return YMBranch::matrixForwardValues(values);
  }
}

// Calculate contributions to the admittance matrix from the branches
void gridpack::state_estimation::SEBranch::setYBus(void)
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
  gridpack::state_estimation::SEBus *bus1 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());
  gridpack::state_estimation::SEBus *bus2 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());
//  if (p_xform) {
//    printf ("from %d-> to %d: p_phase_shift = %f, a = %f+%fi\n", bus1->getOriginalIndex(), bus2->getOriginalIndex(), p_phase_shift, real(a), imag(a) );
//  }
  //p_theta = bus1->getPhase() - bus2->getPhase();
  double pi = 4.0*atan(1.0);
  p_theta = (bus1->getPhase() - bus2->getPhase());
  //printf("p_phase_shift: %12.6f\n",p_phase_shift);
  //printf("p_theta: %12.6f\n",p_theta);
  //printf("p_tap_ratio: %12.6f\n",p_tap_ratio);

}

/**
 * Load values stored in DataCollection object into SEBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 */
void gridpack::state_estimation::SEBranch::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  YMBranch::load(data);

  bool ok = true;
  data->getValue(BRANCH_NUM_ELEMENTS, &p_elems);
  double rvar;
  int ivar;
  std::string svar;
  double pi = 4.0*atan(1.0);
  p_active = false;
  ok = data->getValue(CASE_SBASE, &p_sbase);
  int idx;
  for (idx = 0; idx<p_elems; idx++) {
    bool xform = true;
    xform = xform && data->getValue(BRANCH_X, &rvar, idx);
    p_reactance.push_back(rvar);
    xform = xform && data->getValue(BRANCH_R, &rvar, idx);
    p_resistance.push_back(rvar);
    ok = ok && data->getValue(BRANCH_SHIFT, &rvar, idx);
    rvar = -rvar*pi/180.0; 
    p_phase_shift.push_back(rvar);
    ok = ok && data->getValue(BRANCH_TAP, &rvar, idx);
    p_tap_ratio.push_back(rvar); 
    ok = ok && data->getValue(BRANCH_CKT, &svar, idx);
    p_tag.push_back(svar);
    if (rvar != 0.0) {
      p_xform.push_back(xform);
    } else {
      p_xform.push_back(false);
    }
    ivar = 1;
    data->getValue(BRANCH_STATUS, &ivar, idx);
    p_branch_status.push_back(static_cast<bool>(ivar));
    if (ivar == 1) p_active = true;
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
  }
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::state_estimation::SEBranch::setMode(int mode)
{
  if (mode == YBus) {
    YMBranch::setMode(gridpack::ymatrix::YBus);
  }
  p_mode = mode;
}

/**
 * Return the complex admittance of the branch
 * @return: complex addmittance of branch
 */
gridpack::ComplexType gridpack::state_estimation::SEBranch::getAdmittance(void)
{
  int i;
  gridpack::ComplexType ret(0.0,0.0);
  for (i=0; i<p_elems; i++) {
    gridpack::ComplexType tmp(p_resistance[i], p_reactance[i]);
    if (!p_xform[i] && p_branch_status[i]) {
      tmp = -1.0/tmp;
    } else {
      tmp = gridpack::ComplexType(0.0,0.0);
    }
    ret += tmp;
  }
  return ret;
}

/**
 * Return transformer contribution from the branch to the calling
 * bus
 * @param bus: pointer to the bus making the call
 * @return: contribution to Y matrix from branch
 */
gridpack::ComplexType
gridpack::state_estimation::SEBranch::getTransformer(gridpack::state_estimation::SEBus *bus)
{
  int i;
  gridpack::ComplexType ret(0.0,0.0);
  for (i=0; i<p_elems; i++) {
    gridpack::ComplexType tmp(p_resistance[i],p_reactance[i]);
    gridpack::ComplexType tmpB(0.0,0.5*p_charging[i]);
    if (p_xform[i] && p_branch_status[i]) {
      tmp = -1.0/tmp;
      tmp = tmp - tmpB;
      gridpack::ComplexType a(cos(p_phase_shift[i]),sin(p_phase_shift[i]));
      a = p_tap_ratio[i]*a;
      if (bus == getBus1().get()) {
        tmp = tmp/(conj(a)*a);
      } else if (bus == getBus2().get()) {
        // tmp is unchanged
      }
    } else {
      tmp = gridpack::ComplexType(0.0,0.0);
    }
    ret += tmp;
  }
  return ret;
}

/**
 * Return the contribution to a bus from shunts
 * @param bus: pointer to the bus making the call
 * @return: contribution to Y matrix from shunts associated with branches
 */
gridpack::ComplexType
gridpack::state_estimation::SEBranch::getShunt(gridpack::state_estimation::SEBus *bus)
{
  double retr, reti;
  retr = 0.0;
  reti = 0.0;
  int i;
  for (i=0; i<p_elems; i++) {
    double tmpr, tmpi;
    if (p_shunt[i] && p_branch_status[i]) {
      tmpr = 0.0;
      tmpi = 0.0;
      if (!p_xform[i]) {
        tmpi = 0.5*p_charging[i];
        tmpr = 0.0;
      }
      // HACK: pointer comparison, maybe could handle this better
      if (bus == getBus1().get()) {
        tmpr += p_shunt_admt_g1[i];
        tmpi += p_shunt_admt_b1[i];
      } else if (bus == getBus2().get()) {
        tmpr += p_shunt_admt_g2[i];
        tmpi += p_shunt_admt_b2[i];
      } else {
        // TODO: Some kind of error
      }
    } else {
      tmpr = 0.0;
      tmpi = 0.0;
    }
    retr += tmpr;
    reti += tmpi;
  }
  return gridpack::ComplexType(retr,reti);
}

/**
 * Write output from branches to standard out
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if branch is contributing string to output, false otherwise
 */
bool gridpack::state_estimation::SEBranch::serialWrite(char *string,
    const int bufsize, const char *signal)
{
  gridpack::ComplexType v1, v2, y, s;
  gridpack::state_estimation::SEBus *bus1 = 
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());
  v1 = bus1->getComplexVoltage();
  gridpack::state_estimation::SEBus *bus2 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());
  v2 = bus2->getComplexVoltage();
  y = gridpack::ComplexType(p_ybusr_frwd,p_ybusi_frwd);
  s = v1*conj(y*(v1-v2));
  double p = real(s)*p_sbase;
  double q = imag(s)*p_sbase;
  sprintf(string, "     %6d      %6d      %12.6f         %12.6f\n",
      bus1->getOriginalIndex(),bus2->getOriginalIndex(),p,q);
  return true;
}
/**
 * Add a measurement to the branch
 * @param measurement a measurement struct that will be used to
 * assign
 * internal paramters
 */
void gridpack::state_estimation::SEBranch::addMeasurement(
    gridpack::state_estimation::Measurement measurement)
{
  p_meas.push_back(measurement);
  //TODO: Implement this method
}

/**
  * Return contribution to constraints
  * @param v: voltage at the other bus
  * @param theta: angle difference between two buses
  */
void gridpack::state_estimation::SEBranch::getVTheta(gridpack::state_estimation::SEBus *bus, double *v, double *theta)
{
  gridpack::state_estimation::SEBus *bus1 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());
  double v1 = bus1->getVoltage();
  gridpack::state_estimation::SEBus *bus2 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());
  double v2 = bus2->getVoltage();
  if (bus == bus1) {
     *v = v2;
     *theta = bus1->getPhase() - bus2->getPhase();  
  }  else if (bus == bus2) {
     *v = v1;
     *theta = bus2->getPhase() - bus1->getPhase();  
  }
}
 
/**
  * Return contribution to constraints
  * @param v1, v2: voltages at buses
  * @param theta: angle difference between two buses
  */
void gridpack::state_estimation::SEBranch::getV1V2Theta(gridpack::state_estimation::SEBranch *branch, double *v1, double *v2, double *theta)
{
  gridpack::state_estimation::SEBus *bus1 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());
  *v1 = bus1->getVoltage();
  gridpack::state_estimation::SEBus *bus2 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());
  *v2 = bus2->getVoltage();
  *theta = bus1->getPhase() - bus2->getPhase();  
}
 
/**
 * Configure branches with state estimation parameters. These can be
 * used in other methods
 */
void gridpack::state_estimation::SEBranch::configureSE(void)
{
  // Calculate the number of matrix elements associated witht this branch
  int reference = 1; // TBD: to be read from XML
  gridpack::state_estimation::SEBus *bus1 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());
  gridpack::state_estimation::SEBus *bus2 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());
  int nmeas = p_meas.size(); // Suppose p_meas is the vector of all the measurements on this branch
  int ncnt = 0;
  int i, j, im, jm, nsize;
  for (i=0; i<nmeas; i++) {
    std::string type = p_meas[i].p_type;
    std::string ckt = p_meas[i].p_ckt;
    if (type == "PIJ" || type == "QIJ" || type == "IIJ") {
      int nsize = p_tag.size();
      for (j=0; j<nsize; j++) {
        if (p_tag[j] == ckt) {
          if (!bus1->getReferenceBus()) {
            ncnt += 2;
          } else {  // reference bus, only for dPIJ/DVI
            ncnt++;
          }
          if (!bus2->getReferenceBus()) {
            ncnt += 2;
          } else {  // reference bus, only for dPIJ/DVJ
            ncnt++;
          }
        } 
      }
    }
  }
  p_numElements = ncnt;
}

/**
 * Return number of rows in matrix from component
 * @return number of rows from component
 */
int gridpack::state_estimation::SEBranch::matrixNumRows() const
{
  return p_meas.size();
}

/**
 * Return number of cols in matrix from component
 * @return number of cols from component
 */
int gridpack::state_estimation::SEBranch::matrixNumCols() const
{
  if (p_mode == Jacobian_H) {
    return 0;
  } else if (p_mode == R_inv) {
    return p_meas.size();
  }
}

/**
 * Set row indices corresponding to the rows contributed by this component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @param idx matrix index of row irow
 */
void gridpack::state_estimation::SEBranch::matrixSetRowIndex(int irow, int idx)
{
  if (p_mode == Jacobian_H) {
    if (irow < p_rowJidx.size()) {
      p_rowJidx[irow] = idx;
    } else {
      p_rowJidx.push_back(idx);
    }
  } else if (p_mode == R_inv) {
    if (irow < p_rowRidx.size()) {
      p_rowRidx[irow] = idx;
    } else {
      p_rowRidx.push_back(idx);
    }
  }
}

/**
 * Set column indices corresponding to the columns contributed by this component
 * @param icol index of column contributed by this component (e.g. if component
 * contributes 3 columns then icol is between 0 and 2)
 * @param idx matrix index of column icol
 */
void gridpack::state_estimation::SEBranch::matrixSetColIndex(int icol, int idx)
{
  if (p_mode == Jacobian_H) {
    if (icol < p_colJidx.size()) {
      p_colJidx[icol] = idx;
    } else {
      p_colJidx.push_back(idx);
    }
  } else if (p_mode == R_inv) {
    if (icol < p_colRidx.size()) {
      p_colRidx[icol] = idx;
    } else {
      p_colRidx.push_back(idx);
    }
  }
}

/**
 * Get the row index corresponding to the rows contributed by this component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @return matrix index of row irow
 */
int gridpack::state_estimation::SEBranch::matrixGetRowIndex(int idx)
{
  if (p_mode == Jacobian_H) {
    return p_rowJidx[idx];
  } else if (p_mode == R_inv) {
    return p_rowRidx[idx];
  }
}

/**
 * Get the column index corresponding to the columns contributed by this component
 * @param icol index of column contributed by this component (e.g. if component
 * contributes 3 columns then icol is between 0 and 2)
 * @return matrix index of column icol
 */
int gridpack::state_estimation::SEBranch::matrixGetColIndex(int idx)
{
  if (p_mode == Jacobian_H) {
    if (idx >= p_colJidx.size())
      printf("violation in branch:matrixGetColIndex branch: %d %d size: %d idx: %d\n",
          getBus1OriginalIndex(),getBus2OriginalIndex(),idx,p_colJidx.size());
    return p_colJidx[idx];
  } else if (p_mode == R_inv) {
    return p_colRidx[idx];
  }
}

/**
 * Return the number of matrix values contributed by this component
 * @return number of matrix values
 */
int gridpack::state_estimation::SEBranch::matrixNumValues() const
{
  if (p_mode == Jacobian_H) {
    return p_numElements;
  } else if (p_mode == R_inv) {
    return p_meas.size();
  }
}

/**
 * Return values from a matrix block
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
*/
void gridpack::state_estimation::SEBranch:: matrixGetValues(ComplexType *values, int *rows, int *cols)
{
  if (p_mode == Jacobian_H) {
    SEBus *bus1 = dynamic_cast<SEBus*>(getBus1().get());
    SEBus *bus2 = dynamic_cast<SEBus*>(getBus2().get());
    int nmeas = p_meas.size(); // Suppose p_meas is the vector of all the measurements on this branch
    int ncnt = 0;
    int i, j, im, jm, nsize;
    double v1, v2, theta;
    std::string ckt, type;
    v1 = bus1->getVoltage();
    v2 = bus2->getVoltage();
    theta = bus1->getPhase() - bus2->getPhase();  
    //    int ref = getRef(this);
    for (i=0; i<nmeas; i++) {
      im = matrixGetRowIndex(i);
      ckt = p_meas[i].p_ckt;
      type = p_meas[i].p_type;
      if (type == "PIJ") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (p_tag[j] == ckt) {
            if (!bus1->getReferenceBus()) {
              jm = bus1->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(v1*v2*(p_resistance[j]*sin(theta)
                    -p_reactance[j]*cos(theta)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus1->matrixGetColIndex(1);
              values[ncnt] = gridpack::ComplexType(-v2*(p_resistance[j]*cos(theta)
                    +p_reactance[j]*sin(theta))
                  +2*(p_resistance[j]+p_shunt_admt_g1[j])*v1,0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            } else {  // reference bus, only for dPIJ/DVI
              jm = bus1->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-v2*(p_resistance[j]*cos(theta)
                    +p_reactance[j]*sin(theta))+2*(p_resistance[j]
                    +p_shunt_admt_g1[j])*v1,0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }
            if (!bus2->getReferenceBus()) {
              jm = bus2->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-v1*v2*(p_resistance[j]*sin(theta)
                    -p_reactance[j]*cos(theta)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus2->matrixGetColIndex(1);
              values[ncnt] = gridpack::ComplexType(-v1*(p_resistance[j]*cos(theta)
                    +p_reactance[j]*sin(theta)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            } else {  // reference bus, only for dPIJ/DVJ
              jm = bus2->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-v1*(p_resistance[j]*cos(theta)
                    +p_reactance[j]*sin(theta)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }
          } 
        }
      } else if (type == "QIJ") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (p_tag[j] == ckt) {
            if (!bus1->getReferenceBus()) {
              jm = bus1->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-v1*v2*(p_resistance[j]*cos(theta)
                    +p_reactance[j]*sin(theta)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus1->matrixGetColIndex(1);
              values[ncnt] = gridpack::ComplexType(-v2 * (p_resistance[j] * sin(theta)
                    - p_reactance[j] * cos(theta))-2*(p_reactance[j]+p_shunt_admt_b1[j])*v1,0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            } else {  // reference bus, only for dQIJ/DVI
              jm = bus1->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-v2*(p_resistance[j]*sin(theta)
                    -p_reactance[j]*cos(theta))-2*(p_reactance[j]
                    +p_shunt_admt_b1[j])*v1,0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }
            if (!bus2->getReferenceBus()) {
              jm = bus2->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(v1 * v2 * (p_resistance[j] * cos(theta)
                    + p_reactance[j] * sin(theta)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus2->matrixGetColIndex(1);
              values[ncnt] = gridpack::ComplexType(-v2 * (p_resistance[j] * sin(theta)
                    - p_reactance[j] * cos(theta)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            } else {  // reference bus, only for dQIJ/DVJ
              jm = bus2->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-v2 * (p_resistance[j] * sin(theta)
                    - p_reactance[j] * cos(theta)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }
          } 
        }
      } else if (type == "IIJ") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (p_tag[j] == ckt) {
            double Iij = sqrt((p_resistance[j]*p_resistance[j]+p_reactance[j]*p_reactance[j])
                *(v1*v1+v2*v2-2*v1*v2*cos(theta))); 
            if (!bus1->getReferenceBus()) {
              jm = bus1->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType((p_resistance[j]*p_resistance[j]
                    +p_reactance[j]*p_reactance[j])*v1*v2*sin(theta)/Iij,0.0);  
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus1->matrixGetColIndex(1);
              values[ncnt] = gridpack::ComplexType((p_resistance[j]*p_resistance[j]
                    +p_reactance[j]*p_reactance[j])*(v1-v2*cos(theta))/Iij,0.0);  
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            } else {  // reference bus, only for dIIJ/DVI
              jm = bus1->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType((p_resistance[j]*p_resistance[j]
                    +p_reactance[j]*p_reactance[j])*(v1-v2*cos(theta))/Iij,0.0);  
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }
            if (!bus2->getReferenceBus()) {
              jm = bus2->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-(p_resistance[j]*p_resistance[j]
                    +p_reactance[j]*p_reactance[j])*v1*v2*sin(theta)/Iij,0.0);  
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus2->matrixGetColIndex(1);
              values[ncnt] = gridpack::ComplexType((p_resistance[j]*p_resistance[j]
                    +p_reactance[j]*p_reactance[j])*(v2-v1*cos(theta))/Iij,0.0);  
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            } else {  // reference bus, only for dIIJ/DVJ
              jm = bus2->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType((p_resistance[j]*p_resistance[j]
                    +p_reactance[j]*p_reactance[j])*(v2-v1*cos(theta))/Iij,0.0);  
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }
          } 
        }
      }
    }
  } else if (p_mode == R_inv) {
    int nsize = p_meas.size();
    int i;
    for (i=0; i<nsize; i++) {
      if (p_meas[i].p_deviation != 0.0) {
        values[i] = 1.0/(p_meas[i].p_deviation*p_meas[i].p_deviation);
      } else {
        values[i] = 0.0;
      }
      rows[i] = matrixGetRowIndex(i);
      cols[i] = matrixGetColIndex(i);
    }
  }
}

/**
 * Return values from a vector
 * @param values: pointer to vector values (z-h(x))
 * @param idx: pointer to vector index 
*/
void gridpack::state_estimation::SEBranch:: VectorGetElementValues(ComplexType *values, int *idx)
{
  if (p_mode == Jacobian_H) {
    gridpack::state_estimation::SEBus *bus1 =
      dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());
    gridpack::state_estimation::SEBus *bus2 =
      dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());
    int nmeas = p_meas.size(); // Suppose p_meas is the vector of all the measurements on this branch
    int ncnt = 0;
    int i, j, im, jm, nsize;
    double ret, ret1, ret2;
    double v1,v2,theta;
    v1 = bus1->getVoltage();
    v2 = bus2->getVoltage();
    theta = bus1->getPhase() - bus2->getPhase();  
    for (i=0; i<nmeas; i++) {
      vectorGetElementIndices(idx);
      if (p_meas[i].p_type == "PIJ") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (p_tag[j] == p_meas[i].p_ckt) {
            ret1 =  v1*v1* (p_resistance[j] + p_shunt_admt_g1[j])
              - v1*v2*(p_resistance[j]*cos(theta) + p_reactance[j]*sin(theta));
          }
        }
        vectorGetElementIndices(idx);
        values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-ret1),0.0);
        //         values] = gridpack::ComplexType(static_cast<double>(idx1+idx2),0.0);
        //         values[ncnt] = p_meas[i].p_value-ret1;
        ncnt++;
      } else if (p_meas[i].p_type == "QIJ") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (p_tag[j] == p_meas[i].p_ckt) {
            ret2 = - v1*v1* (p_reactance[j] + p_shunt_admt_b1[j])
              - v1*v2*(p_resistance[j]*sin(theta) - p_reactance[j]*cos(theta));
          }
        }
        vectorGetElementIndices(idx);
        values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-ret2),0.0);
        //         values[ncnt] = p_meas[i].p_value-ret2;
        ncnt++;
      } else if (p_meas[i].p_type == "IIJ") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (p_tag[j] == p_meas[i].p_ckt) {
            ret1 =  v1*v1* (p_resistance[j] + p_shunt_admt_g1[j])
              - v1*v2*(p_resistance[j]*cos(theta) + p_reactance[j]*sin(theta));
            ret2 = -v1*v1* (p_reactance[j] + p_shunt_admt_b1[j])
              - v1*v2*(p_resistance[j]*sin(theta) - p_reactance[j]*cos(theta));
          }
        }
        ret = sqrt(ret1*ret1+ret2*ret2)/v1;
        vectorGetElementIndices(idx);
        //         values[ncnt] = p_meas[i].p_value-ret;
        values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-ret),0.0);
        ncnt++;
      }
    } 
  } else if (p_mode == R_inv) {
  }
}
