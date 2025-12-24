#pragma once

#include "types.h"

int solver_run(real* restrict acc, real* restrict pot, real* restrict rho, real* restrict prs, const real* restrict mas, const real* restrict rad, const real* restrict pos, const real* restrict vel, real scale_factor, unsigned int N, unsigned int step_number);
