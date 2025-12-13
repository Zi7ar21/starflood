#pragma once

#include "types.h"

int solver_run(real* volatile pot, real* volatile acc, const real* volatile mas, const real* volatile pos, unsigned int N, unsigned int step_number);
