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
  GAIN:

  input  : u 
  output : y
      
                 
                
                              ymax
                             -----
        -------------       /
        |           |      /
  u ----|     K     |----------- y
        |           |    /
        -------------   /
                    ---
                    ymin
        
        

   Output:
    y = K*u,  ymin <= y <= ymax

    Notes:
     This block can be used in various ways
     i) Simple Gain - Set limits = infty
    ii) Simple limiter - Set the gain K to 1.0, set the limits
   iii) Gain with limiter - Set the gain and limits appropriately

*/
class Gain
{
 public:
  Gain();

  /**
     GETOUTPUT - Returns the block output

     INPUTS:
       u       - block input
       y       - block output
  **/
  double getoutput(double u);
  
  /**
     SETPARAMS - Set the gain and limits

     INPUTS:
       K          Gain
       ymin       Min. limit for output y
       ymax       Max. limit for output y
  **/
  void setparams(double K, double ymin,double ymax);

  private:

  double p_K; // Gain
  double p_ymin; // Lower limit on output
  double p_ymax; // Upper limit on output
};

/*
  SLOPE:

  input  : u 
  output : y

      
                 
                             
                             
        -------------       
        |           |      
  u ----|   Slope   |
        |  Control  |----------- y
        |           |    
        -------------   
                    
                    
        
        

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


#endif
