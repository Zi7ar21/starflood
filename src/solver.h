#pragma once

#include "types.h"

real probe_potential(const real* restrict mas, const real* restrict pos, const real* probe_pos, unsigned int N, unsigned int step_number);

int solver_run(real* restrict acc, real* restrict pot, real* restrict rho, real* restrict prs, const real* restrict mas, const real* restrict rad, const real* restrict pos, const real* restrict vel, unsigned int N, unsigned int step_number);
