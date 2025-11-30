#pragma once

#include <stddef.h>

#include "types.h"

typedef struct {
	unsigned int N; // Number of bodies in the simulation

	void* mem; // Shared Memory

	real* pot; // Potential Energy (Scalar)
	real* kin; // Kinetic Energy (Scalar)
	real* mas; // Mass (Scalar)
	real* pos; // Position (3-Vector)
	real* vel; // Velocity (3-Vector)
	real* acc; // Acceleration (3-Vector)
} simulation_t;

int simulation_init(simulation_t* simulation, unsigned int N);

int simulation_step(simulation_t* simulation);

int simulation_free(simulation_t* simulation);
