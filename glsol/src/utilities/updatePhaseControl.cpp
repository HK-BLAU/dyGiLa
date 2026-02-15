#define USE_MPI 
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>
#include <algorithm>
#include <cmath>

#include "plumbing/hila.h"
#include "glsol.hpp"
#include "matep.hpp"

void glsol::updatePhaseControl(real_t currentVolumeFraction) {
    
    // calculate error
    real_t volumeError = currentVolumeFraction - config.targetVolumeFraction;
    
    // add to history (obsolete currently)
    volumeHistory.push_back(currentVolumeFraction);
    if (volumeHistory.size() > config.historySize) {
        volumeHistory.erase(volumeHistory.begin());
    }
    
    // proportional term
    real_t P_term = config.Kp * volumeError;
    
    // integral term
    errorIntegral += volumeError * config.dt;
    real_t I_term = config.Ki * errorIntegral;
    
    // derivative term
    real_t D_term = 0.0;
    if (lastVolumeError != 0.0) {
        real_t errorDerivative = (volumeError - lastVolumeError) / config.dt;
        D_term = config.Kd * errorDerivative;
    }
    
    // total temperature change from PID
    real_t deltaT = P_term + I_term + D_term;
    
    // limit to maximum allowed change
    deltaT = std::max(-config.maxTempChange, std::min(config.maxTempChange, deltaT));
    
    // apply temperature change uniformly
    //const real_t Tcp_mK = MP.Tcp_mK(config.Inip);
    onsites(ALL) {
        T[X] += deltaT;
    }
    
    // update last error for next iteration
    lastVolumeError = volumeError;
    
    // debug output
    if (hila::myrank() == 0) {
        hila::out0 << " VolumeError=" << volumeError
                   << " P=" << P_term << " I=" << I_term << " D=" << D_term
                   << " deltaT=" << deltaT << std::endl;
    }
}