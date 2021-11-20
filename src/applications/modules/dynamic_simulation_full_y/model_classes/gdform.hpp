/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   gdform.hpp
 * @author Renke Huang
 * @Last modified:   Jan. 16, 2021
 * 
 * @brief  
 * 
 * 
 */


#ifndef _gdform_h_
#define _gdform_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_generator_model.hpp"

namespace gridpack {
namespace dynamic_simulation {
class GridFormingGenerator : public BaseGeneratorModel
{
  public:
    /**
     * Basic constructor
     */
    GridFormingGenerator();

    /**
     * Basic destructor
     */
    virtual ~GridFormingGenerator();

    /**
     * Load parameters from DataCollection object into generator model
     * @param data collection of generator parameters from input files
     * @param index of generator on bus
     * TODO: might want to move this functionality to BaseGeneratorModel
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);
   
    /**
     * Initialize generator model before calculation
     * @param mag voltage magnitude
     * @param ang voltage angle
     * @param ts time step 
     */
    void init(double mag, double ang, double ts);

    /**
     * Return contribution to Norton current
     * @return contribution to Norton vector
     */
    gridpack::ComplexType INorton();

    /**
     * Return Norton impedence
     * @return value of Norton impedence
     */
    gridpack::ComplexType NortonImpedence();

    /**
     * Predict new state variables for time step
     * @param t_inc time step increment
     * @param flag initial step if true
     */
    void predictor(double t_inc, bool flag);

    /**
     * Correct state variables for time step
     * @param t_inc time step increment
     * @param flag initial step if true
     */
    void corrector(double t_inc, bool flag);

    /**
     * Predict part calculate current injections
     * @param flag initial step if true
     */
    void predictor_currentInjection(bool flag);

    /**
     * Corrector part calculate current injections
     * @param flag initial step if true
     */
    void corrector_currentInjection(bool flag);

    /**
     * Set voltage on each generator
     */
    void setVoltage(gridpack::ComplexType voltage);

    /**
	* return true if trip generator successfully
	* 
	*/
    bool tripGenerator();
	
	/**
	* return true if modify the generator parameters successfully
	* input controlTyp: 0: GFI mp adjust; 1: GFI mq adjust; 2: GFI Pset adjust; 3: GFI Qset adjust; others: invalid 
    * input newParValScaletoOrg:  GFI new parameter scale factor to the very initial parameter value at the begining of dynamic simulation
	* 
	*/
    bool applyGeneratorParAdjustment(int controlType, double newParValScaletoOrg);

    /**
     * Write out generator state
     * @param string buffer that contains output
     * @param bufsize size of string buffer
     * @param signal character string used to determine behavior
     */
    bool serialWrite(char* string, const int bufsize, const char* signal);

    /**
     * return a vector containing any generator values that are being
     * watched
     * @param vals vector of watched values
     */
    void getWatchValues(std::vector<double> &vals);

  private:

    double p_sbase;
    double p_pg, p_qg;
    int p_status;
	bool p_tripped;
	bool bmodel_debug;
	bool bCurrentLimitFlag;
    
    double MVABase;
    double XL, Ts, Vset, mq, kpv, kiv, Emax, Emin, mp, kppmax, kipmax, Pset, Pmax, Pmin, Imax;
    double fset, Ra;
	double delta_omega_lim;
	double mp_org, mq_org, Vset_org, Pset_org;
    
    //double Vterm, Theta, Ir, Ii;
    
    double x1E_0, x2d_0, x3_0, x4_0, xvterm_0, xp_0, xq_0;
    double x1E_1, x2d_1, x3_1, x4_1, xvterm_1, xp_1, xq_1;
    double dx1E_0, dx2d_0, dx3_0, dx4_0, dxvterm_0, dxp_0, dxq_0;
    double dx1E_1, dx2d_1, dx3_1, dx4_1, dxvterm_1, dxp_1, dxq_1;

    double E_term, IrNorton, IiNorton, omega, genP, genQ;
	double Poutctrl;
	double Vterm_delay, genP_delay, genQ_delay;
    
    gridpack::ComplexType p_INorton;
	gridpack::ComplexType p_Norton_Ya;

    double presentMag, presentAng;
	double Ir, Ii, Iinjr, Iinji, Vterm, Theta, Igen_mag;
	double B, G;

    std::string p_ckt;
    int p_bus_id;

    friend class boost::serialization::access;

    template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
        ar &
          & p_sbase
          & p_pg & p_qg
          & p_status
          & MVABase 
          & p_INorton
          & p_bus_id;
      }

};
}  // dynamic_simulation
}  // gridpack
#endif
