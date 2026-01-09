#include "solvers.h"

#include "common.h"
#include "config.h"
#include "simulation.h"
#include "gravity.h"
#include "sph.h"

int solvers_run(sim_t* restrict sim_ptr) {
	sim_t sim = *sim_ptr;

	solve_gravity(&sim);

	#ifdef ENABLE_SPH
	solve_sph(&sim);
	#endif

	return STARFLOOD_SUCCESS;
}
