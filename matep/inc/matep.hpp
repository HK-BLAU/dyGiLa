/*
 * 
 * This is the *.hpp header file of the Strong Coupling Correction Object (SCCO) 
 * or Material Parameters Object (Matep).
 * 
 * Member functions are declared at here, and then defined at Matep.cpp.
 *
 * author: Quang. Zhang (timohyva@github)
 *
 */ 

#ifndef MATEP_HPP
#define MATEP_HPP

//#include <string>
#include <iostream>
#include <cstddef>
#include <cmath>
#include <vector>


using real_t = float; // or double
namespace matep {
class Matep {
public:
        Matep() {};                                                  // default constructor

  // ************************************************************************** //
  // >>>>>>>>>>>        interfaces of dimensional qualities        <<<<<<<<<<<< //
  // ************************************************************************** //
		#pragma hila loop_function
        real_t Fa0p(real_t p);                 //dimensionless Landau Coefficient
		
		#pragma hila loop_function
        real_t Tcp(real_t p);                  // in unit of Kelvin
		
		#pragma hila loop_function
        real_t Tcp_mK(real_t p);               // in unit of mK
		
		#pragma hila loop_function
        real_t tauQP(real_t p, real_t T);      // Fitted QP releaxing time in \mus
		
		#pragma hila loop_function
        real_t mEffp(real_t p);                // quisiparticle effective mass 
		
		#pragma hila loop_function
        real_t vFp(real_t p);                  // Fermi velocity
		
		#pragma hila loop_function
        real_t xi0p(real_t p);                 // zero Temperature coherent length
		
		#pragma hila loop_function
        real_t xi0GLp(real_t p);               // zero temperature GL choherent length
		
		#pragma hila loop_function
        real_t tGL(real_t p);                  // GL time, time unit
		
		#pragma hila loop_function
        double N0p(real_t p);                  // deisty of state on Fermi surface
        

  // ************************************************************************* //
  // >>>>  interfaces of dimensionless coeficients; SC-correction parts: <<<<< //
  // ************************************************************************* //
		#pragma hila loop_function
        real_t alpha_td(real_t p, real_t T);
		
		#pragma hila loop_function
        real_t beta1_td(real_t p, real_t T);
		
		#pragma hila loop_function
        real_t beta2_td(real_t p, real_t T);
		
		#pragma hila loop_function
        real_t beta3_td(real_t p, real_t T);
		
		#pragma hila loop_function
        real_t beta4_td(real_t p, real_t T);
		
		#pragma hila loop_function
        real_t beta5_td(real_t p, real_t T);

        // damping term coefficient gamma, in unit of tGL^-1
		#pragma hila loop_function
        real_t gamma_td(real_t p, real_t T);        

        // Quadratic H-term coefficient gz_td
		#pragma hila loop_function
        real_t gz_td(real_t p);        

  // >>>>>>>>>    interfaces for beta_A, beta_B, gaps and tAB_RWS    <<<<<<<<< //
		#pragma hila loop_function
        real_t beta_A_td(real_t p, real_t T);
		
		#pragma hila loop_function
        real_t beta_B_td(real_t p, real_t T);
		
		#pragma hila loop_function
        real_t gap_A_td(real_t p, real_t T);
		
		#pragma hila loop_function
        real_t gap_B_td(real_t p, real_t T);
		
		#pragma hila loop_function
        real_t gap_td(real_t p, real_t T);      // gap for given p,T, and showing message

		#pragma hila loop_function
        real_t tAB_RWS(real_t p);               

  // >>>>>>> interfaces for f_{A}, f_{B}, in unit of (1/3)(Kb Tc)^2 N(0) <<<<< //
		#pragma hila loop_function
        real_t f_A_td(real_t p, real_t T);
		
		#pragma hila loop_function
        real_t f_B_td(real_t p, real_t T);
        

  // >>>>>>>>>>>>>>> menmber funcition Levi_Civita symbol <<<<<<<<<<<<<<<<<<<< //
		#pragma hila loop_function
        real_t epsilon(int al, /*alpha*/
         	       int be, /*beta*/
		       int ga /*gamma*/);
  
private:		
		#pragma hila loop_function
        // linear interpolation function:
        real_t lininterp(const real_t *cX_arr, real_t p);
};
}
#endif
