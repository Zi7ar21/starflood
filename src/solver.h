#pragma once

#include "types.h"

int solver_run(real* volatile acc, real* volatile pot, const real* volatile mas, const real* volatile rad, const real* volatile pos, unsigned int N, unsigned int step_number);
