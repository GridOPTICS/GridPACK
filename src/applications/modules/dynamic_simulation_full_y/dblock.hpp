/**
 * @file dblock.hpp
 * @brief Header file for discontinuous blocks
 *
 */

#ifndef DBLOCK_HPP
#define DBLOCK_HPP

#include <math.h>
#include <algorithm>

/*
  GAINLIMITER:

  input  : u 
  output : y
      
                 
                
                              dymax,ymax
                             -----
        -------------       /
        |           |      /
  u ----|     K     |----------- y
        |           |    /
        -------------   /
                    ---
                    dymin, ymin
        
        

   Output:
    y = K*u,  ymin <= y <= ymax,dymin <= dydt <= dymax 

    Notes:
     This block can be used in various ways
     i) Simple Gain - Set limits = infty
    ii) Simple limiter - Set the gain K to 1.0, set the limits
   iii) Gain with limiter - Set the gain and limits (ymin,ymax) appropriately
    iv) Rate limiter - Set the gain and rate limits (dymin,dymax) appropriately

*/
class GainLimiter
{
 public:
  GainLimiter();

  /**
     GETOUTPUT - Returns the block output

     INPUTS:
       u       - block input
       y       - block output
  **/
  double getoutput(double u);

  /**
     GETOUTPUT - Returns the block output, applies dynamic limits

     INPUTS:
       u       - block input
       ymin    - Min. limit on output
       ymax    - Max. limit on output

     OUPUTS:
       y       - block output
  **/
  double getoutput(double u, double ymin, double ymax);

  
  /**
     SETPARAMS - Set the gain and limits

     INPUTS:
       K          Gain
       ymin       Min. limit for output y
       ymax       Max. limit for output y
  **/
  void setparams(double K, double ymin,double ymax);

  /**
     SETPARAMS - Set the gain and limits
     
     INPUTS:
       K          Gain
       ymin       Min. limit for output y
       ymax       Max. limit for output y
       dymin      Min. limit for rate of change of y
       dymax      Max. limit for rate of change of y
  **/
  void setparams(double K, double ymin,double ymax,double dymin, double dymax);

  /**
     GETOUTPUT - Returns the block output

     INPUTS:
       u       - block input
      dt       - time-step
    isupdate   - move to next step?
  **/
  double getoutput(double u, double dt, bool isupdate);

  private:

  double p_K; // Gain
  double p_ymin; // Lower limit on output
  double p_ymax; // Upper limit on output
  double p_dymin; // Lower limit on rate of change of output
  double p_dymax; // Upper limit on rate of change of output
  double p_dt;    // time-step
  double p_yprev; // Output at the previous step
};

/*
  SLOPE:

  input  : u 
  output : y
        ----------------           
        |        ymax  |
        |       ----   |
        |      /       |      
  u ----|     /dy_du   |----- y
        | ___/         |
        | ymin         |    
        ---------------   
                    

   Output:
    if u <= u0,     y = ya
    if u >= u1,     y = yb
    else:
      y = y0 + dy_du*(u - u0)

    where, the slope dy_du is calculated from the points (u0,y0), (u1,y1) as

     dy_du = (y1 - y0)/(u1 - u0)
*/
class Slope
{
 public:
  Slope();

  /**
     GETOUTPUT - Returns the block output

     INPUTS:
       u       - block input
       y       - block output
  **/
  double getoutput(double u);
  
  /**
     SETPARAMS - Set the gain and the slop points

     INPUTS:
       u       size 2 array of input pair [u0,u1] for slope calculation
       y       size 2 array of output pair [y0,y1] for slope calculation
       ya      Output y when u <= u0
       yb      Output y when u >= u1
  **/
  void setparams(double *, double *,double, double);

  private:

  double p_u[2]; // Input points [u0,u1] for slope calculation
  double p_y[2]; // Output points [y0,y1] for slope calculatimn
  double p_ya;   // Output y when u <= u0
  double p_yb;   // Output y when u >= u1
  double p_dydu; // slope
};

/**
  Deadband:

  input  : u 
  output : y
                   
        ---------------       
        |           /  |      
  u ----|   _______/   |
        |  /u1    u2   |----------- y
        | /            |    
        ---------------   
                    

   Output:
    if u <= u1,     y = u
    if u >= u2,     y = u
    else:
      y = 0
**/
class Deadband
{
  public:
  Deadband();

    /**
     GETOUTPUT - Returns the block output

     INPUTS:
       u       - block input
       y       - block output
  **/
  double getoutput(double u);
  
  /**
     SETPARAMS - Set the gain and the slop points

     INPUTS:
       u1      lower limit for deadband
       u2      upper limit for deadband
  **/
  void setparams(double, double);
  private:
  double p_u1,p_u2; // Lower and upper limits for deadband
};

/*
  PIECEWISESLOPE: Implements a piecewise linear slope function

  input  : u 
  output : y
        ----------------           
        |               |
        |   PIECEWISE   |
        |     SLOPE     |      
  u ----|               |----- y
        |               |    
        ----------------   
                    

   Output:
   The output y is decided on which slope segment u lies on
   See the description for Slope class, i.e., if u lies on
   segment i with boundary points u_i and u_i{i+1} output y
   is
      y = y_i + (dy_du)_{i}*(u - u_i)

    where, the slope dy_du is calculated from the points (u_i,y_i), (u_{i+1},y_{i+1}) as

     (dy_du)_{i+1} = (y_{i+1} - y_i)/(u_{i+1} - u_i)
     
*/
class PiecewiseSlope
{
 public:
  PiecewiseSlope();

  /**
     GETOUTPUT - Returns the block output

     INPUTS:
       u       - block input
       y       - block output
  **/
  double getoutput(double u);
  
  /**
     SETPARAMS - Set the gain and the slope points

     INPUTS:
       n       number of slope segments (max 4.)
       u       array of size n+1 input points [u0,u1,...] for the piecewise function
       y       array of size n+1 output points [y0,y1,...]] for slope calculation

       Notes:
         u <= u0 => y = y0
	 u >= un => y = yn

	 Max 4. segments (n = 4) supported currently.
	 
	 Segments with [u_i,y_i] = [0,0] are ignored
  **/
  void setparams(int, double *, double *);

  double init_given_y(double yout);

  private:
    Slope p_slopes[4];
    double p_u[5]; // Input points for slope calculations
    double p_y[5]; // Output points for slope calculations
    int   p_n;     // Number of segments 
    bool  p_increasing; // Increasing (true) or decreasing (false) function
};


#endif
