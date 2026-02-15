#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>

#include "plumbing/hila.h"
#include "plumbing/fft.h"

#include "glsol.hpp"
#include "matep.hpp"


void glsol::next_UniT_Hfield_PID() {

  static hila::timer next_timer("timestep");
  Field<phi_t> deltaPi;
  Field<Vector<3,Complex<real_t>>> djAaj;

  int bc = config.boundaryConditions;

  next_timer.start();

  std::initializer_list<int> coordsList {0,0,0};
  const CoordinateVector originpoints(coordsList);

  onsites(ALL) {
    matep::Matep MPonsites;
    
    real_t gapa = MPonsites.gap_A_td(config.Inip, T[X]);
    real_t gapb = MPonsites.gap_B_td(config.Inip, T[X]);

    A[X] += config.dt * pi[X];

    if (bc == 1) {
      if ((X.coordinate(e_x) == 0) || (X.coordinate(e_x) == 1)) {
        foralldir(d1) foralldir(d2) {
          if (d1 == d2) {
            A[X].e(d1,d2).re = 1.0;
            A[X].e(d1,d2).im = 0.0;
          }
          else {
            A[X].e(d1,d2).re = 0.0;
            A[X].e(d1,d2).im = 0.0;
          }	
        }
        A[X] = gapb * A[X]/sqrt(3.0);
      }
      else if ((X.coordinate(e_x) == (config.lx - 1)) || (X.coordinate(e_x) == (config.lx - 2))) {
        foralldir(d1) foralldir(d2) {
          if (d1 == 0 && d2 == 0) {
            A[X].e(d1,d2).re = 1.0;
            A[X].e(d1,d2).im = 0.0;
          }
          else if (d1 == 0 && d2 == 1) {
            A[X].e(d1,d2).re = 0.0;
            A[X].e(d1,d2).im = 1.0;
          }
          else {
            A[X].e(d1,d2).re = 0.0;
            A[X].e(d1,d2).im = 0.0;
          }
        }
        A[X] = gapa * A[X]/sqrt(2.0);
      }
    }
    else if (bc == 2) {
      if (X.coordinate(e_x) == 0 || X.coordinate(e_x) == (config.lx - 1) ||
          X.coordinate(e_x) == 1 || X.coordinate(e_x) == (config.lx - 2) ||
          X.coordinate(e_y) == 0 || X.coordinate(e_y) == (config.ly - 1) ||
          X.coordinate(e_y) == 1 || X.coordinate(e_y) == (config.ly - 2) ||
          X.coordinate(e_z) == 0 || X.coordinate(e_z) == (config.lz - 1) ||
          X.coordinate(e_z) == 1 || X.coordinate(e_z) == (config.lz - 2)) {
        A[X] = 0.0;	    
      }
    }
  }

  onsites(ALL) {
    matep::Matep MPonsites;
    
    real_t beta0 = MPonsites.alpha_td(config.Inip, T[X]);
    real_t beta1 = MPonsites.beta1_td(config.Inip, T[X]);
    real_t beta2 = MPonsites.beta2_td(config.Inip, T[X]);
    real_t beta3 = MPonsites.beta3_td(config.Inip, T[X]);
    real_t beta4 = MPonsites.beta4_td(config.Inip, T[X]);
    real_t beta5 = MPonsites.beta5_td(config.Inip, T[X]);

    auto AxAt = A[X]*A[X].transpose();
    auto AxAd = A[X]*A[X].dagger();

    deltaPi[X] = - beta0*A[X]
      - 2.0*beta1*A[X].conj()*AxAt.trace()
      - 2.0*beta2*A[X]*AxAd.trace()
      - 2.0*beta3*AxAt*A[X].conj()
      - 2.0*beta4*AxAd*A[X]
      - 2.0*beta5*A[X].conj()*A[X].transpose()*A[X]
      - MPonsites.gz_td(config.Inip)*H[X]*(H[X].transpose()*A[X]);
  }

  onsites(ALL) {
    djAaj[X] = 0;
    foralldir(j) {
      djAaj[X] += A[X + j].column(j) - A[X - j].column(j);
    }
  }

  onsites(ALL) {
    phi_t mat;
    foralldir(d) {
      auto col = djAaj[X+d] - djAaj[X-d];
      for (int i=0; i<NDIM; i++) mat.e(i,d) = col[i];
    }

    deltaPi[X] += (1.0/(2.0*(config.dx*config.dx)))*mat;
  }

  onsites(ALL) {
    deltaPi[X] += (1.0/(config.dx*config.dx)) * (A[X + e_x] + A[X - e_x]
                                                     + A[X + e_y] + A[X - e_y]
                                                     + A[X + e_z] + A[X - e_z]
                                                     - 6.0*A[X]);
  }

  if (t < config.tdif) {
    // initial diffusive phase
    pi[ALL] = deltaPi[X]/(config.difFac);
    t += config.dt/config.difFac;
  }
  else if (t < config.tdis && config.gamma.squarenorm() > 0 )
    {
      pi[ALL] = pi[X] + (deltaPi[X] - 2.0 * config.gamma * pi[X])*(config.dt);
      t += config.dt;
    }
  else
    {
      pi[ALL] = pi[X] + deltaPi[X]*config.dt;
      t += config.dt;
    }

  if (config.Tstabilization == 1 && t >= config.Tstabilization_start) {
    real_t currentVolumeFraction = phaseCounting();
    updatePhaseControl(currentVolumeFraction);
  }

  next_timer.stop();

} // next_UniT_Hfield_PID() ends here