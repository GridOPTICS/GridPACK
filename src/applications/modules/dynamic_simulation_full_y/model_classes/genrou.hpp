/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   genrou.hpp
 * @author Shuangshuang Jin 
 * @Last modified: Oct 9, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef _genrou_h_
#define _genrou_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_generator_model.hpp"

namespace gridpack {
namespace dynamic_simulation {
class GenrouGenerator : public BaseGeneratorModel
{
  public:
    /**
     * Basic constructor
     */
    GenrouGenerator();

    /**
     * Basic destructor
     */
    virtual ~GenrouGenerator();

    /**
     * Load parameters from DataCollection object into generator model
     * @param data collection of generator parameters from input files
     * @param index of generator on bus
     * TODO: might want to move this functionality to BaseGeneratorModel
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Saturation function
     * @ param x
     */
    double Sat(double x);
    
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
	* return true if trip generator successfully
	* 
	*/
    bool tripGenerator();
	
    /**
     * Set voltage on each generator
     */
    void setVoltage(gridpack::ComplexType voltage);

    /** 
     * Get the value of the field voltage parameter
     * @return value of field voltage
     */
    double getFieldVoltage();

    /**
     * Write out generator state
     * @param signal character string used to determine behavior
     * @param string buffer that contains output
     */
	 // Yuan commented below 20201011
    //void write(const char* signal, char* string);
	// Yuan commend end
	
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

    /**
     * Set internal state parameter in generator
     * @param name character string corresponding to state variable
     * @param value new value for state parameter
     * @return false if no variable corresponding to name is found
     */
    bool setState(std::string name, double value);

    /**
     * Get internal state parameter in generator
     * @param name character string corresponding to state variable
     * @param value current value for state parameter
     * @return false if no variable corresponding to name is found
     */
    bool getState(std::string name, double *value);

  private:

    double p_sbase;
    double p_pg, p_qg;
    int p_status;
	bool p_tripped;
    
    double MBase;
    double H, D, Ra, Xd, Xq, Xdp, Xdpp, Xl, Xqp, Xqpp;
    double Tdop, Tdopp, Tqopp, Tqop, S10, S12;
	double genP, genQ;
    
    double Vterm, Theta, Ir, Ii;
    
    double x1d, x2w, x3Eqp, x4Psidp, x5Psiqp, x6Edp;
    double x1d_1, x2w_1, x3Eqp_1, x4Psidp_1, x5Psiqp_1, x6Edp_1;
    double dx1d, dx2w, dx3Eqp, dx4Psidp, dx5Psiqp, dx6Edp;
    double dx1d_1, dx2w_1, dx3Eqp_1, dx4Psidp_1, dx5Psiqp_1, dx6Edp_1;
    double Id, Iq;
    
    double Efd, LadIfd, Pmech;
    
    double B, G;

    double IrNorton, IiNorton;
    
    gridpack::ComplexType p_INorton;
	gridpack::ComplexType p_Norton_Ya;

    boost::shared_ptr<BaseGovernorModel> p_governor;
    boost::shared_ptr<BaseExciterModel> p_exciter;

    double presentMag, presentAng;

    double Efdinit, Pmechinit;
	
	bool enableSat, printFlag;
            
    ////////////////////////////////////////////////////

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
          & MBase & Ra
          & p_INorton
          & p_bus_id;
      }

};
}  // dynamic_simulation
}  // gridpack
#endif
