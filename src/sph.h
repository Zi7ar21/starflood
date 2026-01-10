#pragma once

#include "simulation.h"

#ifdef ENABLE_SPH
int solve_sph(const sim_t* restrict sim_ptr);
#endif
