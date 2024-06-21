#include <genplb.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Genplb::Genplb(void)
{
  nxgen   = 0; // Number of variables for this model when integration type is implicit or implicit-explicit
  has_playback_file = false;
}

void Genplb::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxgen = 0;
  *nvar = nxgen;
}

Genplb::~Genplb(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseGeneratorModel
 */
void Genplb::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTGenModel::load(data,idx); // load parameters in base generator model
  gridpack::ComplexType Zsource;
  char line[1000];
  char *out;

  // load parameters for the model type
  data->getValue(BUS_NUMBER, &bid);
  data->getValue(GENERATOR_PLAYBACK_FILE,&playback_file,idx);
  data->getValue(GENERATOR_PLAYBACK_ISCALE,&Iscale,idx);

  // Open file
  fp = fopen(playback_file.c_str(), "r");
  
  if(!fp) {
    has_playback_file = false;
    printf("Could not open generator playback file %s\n",playback_file.c_str());
    exit(1);
  } else {
    has_playback_file = true;
  
    // Ignore first line (header info)
    out = fgets(line,1000,fp);
    // Read second line (Time = 0)
    out = fgets(line,1000,fp);
    
    sscanf(line, "%lf, %lf, %lf, %lf",&file_time1,&file_I1[0],&file_I1[1],&file_I1[2]);
    
    // Read third line (first time instant
    out = fgets(line,1000,fp);

    sscanf(line, "%lf, %lf, %lf, %lf",&file_time2,&file_I2[0],&file_I2[1],&file_I2[2]);
  }

  
  data->getValue(GENERATOR_ZSOURCE,&Zsource,idx);
  p_Rs = real(Zsource);
  p_Xdp = imag(Zsource);

  p_L = p_Xdp/OMEGA_S;

}

/**
 * Initialize generator model before calculation
 * @param [output] xin - array where initialized generator variables should be set
 */
void Genplb::init(gridpack::RealType* xin)
{
  double IGD,IGQ; // Machine currents in cartesian coordinates
  double Pg, Qg;  // Generator real and reactive power
  double dw=0.0;  // Initial machine speed deviation
  gridpack::RealType *x = xin+offsetb; // generator array starts from this location
  double VD,VQ;

  Pg = pg/sbase;
  Qg = qg/sbase;

  VD = p_Vm0*cos(p_Va0);
  VQ = p_Vm0*sin(p_Va0);

  gridpack::ComplexType V = gridpack::ComplexType(VD,VQ);
  gridpack::ComplexType S = gridpack::ComplexType(Pg,Qg);
  gridpack::ComplexType I;
  gridpack::ComplexType Z = gridpack::ComplexType(p_Rs,p_Xdp);

  I = conj(S/V);

  Im = abs(I);
  Ia = arg(I);

  double ia,ib,ic;

  ia = Im*sin(Ia);
  ib = Im*sin(Ia - 2*PI/3.0);
  ic = Im*sin(Ia + 2*PI/3.0);

  p_iabc[0] = ia;
  p_iabc[1] = ib;
  p_iabc[2] = ic;
}

/**
 * Write output from generators to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Genplb::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Genplb::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto generators
 * @param values array containing generator state variables
*/
void Genplb::setValues(gridpack::RealType *values)
{
  gridpack::RealType *x = values+offsetb; // generator array starts from this location

}

/**
   Prestep function
*/
void Genplb::preStep(double time ,double timestep)
{

}

/**
   Poststep function
*/
void Genplb::postStep(double time)
{
}



/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Genplb::vectorGetValues(gridpack::RealType *values)
{

}

  /**
   * Return the generator current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
void Genplb::getCurrent(double *ia, double *ib, double *ic)
{
  char line[100];
  char *out;

  if(!has_playback_file) {
    p_iabc[0] = Im*sin(OMEGA_S*p_time + Ia);
    p_iabc[1] = Im*sin(OMEGA_S*p_time + Ia - TWOPI_OVER_THREE);
    p_iabc[2] = Im*sin(OMEGA_S*p_time + Ia + TWOPI_OVER_THREE);
  } else {
    if(file_time2 < p_time) {
      file_time1 = file_time2;
      file_I1[0] = file_I2[0];
      file_I1[1] = file_I2[1];
      file_I1[2] = file_I2[2];

      out = fgets(line,1000,fp);
      if(out == NULL) {
	printf("End of file");
	exit(1);
      } else {
	sscanf(line, "%lf, %lf, %lf, %lf",&file_time2,&file_I2[0],&file_I2[1],&file_I2[2]);
      }
    }

    /* Do interpolation */
    double I_slope[3];
    I_slope[0] = (file_I2[0] - file_I1[0])/(file_time2 - file_time1);
    I_slope[1] = (file_I2[1] - file_I1[1])/(file_time2 - file_time1);
    I_slope[2] = (file_I2[2] - file_I1[2])/(file_time2 - file_time1);

    double dt = p_time - file_time1;

    // Calculate current. Negative sign because it is a load
    p_iabc[0] = -(file_I1[0] + dt*I_slope[0])/Iscale; 
    p_iabc[1] = -(file_I1[1] + dt*I_slope[1])/Iscale; 
    p_iabc[2] = -(file_I1[2] + dt*I_slope[2])/Iscale; 
    
  }     

  *ia = p_iabc[0];
  *ib = p_iabc[1];
  *ic = p_iabc[2];
}

/**
 * Return the global location for the generator current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void Genplb::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = -1;
}


/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values
 */
int Genplb::matrixNumValues()
{
  int numVals  = 0;

  return numVals;
}

/**
 * Return values from Jacobian matrix
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void Genplb::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;
  
  *nvals = ctr;
}
