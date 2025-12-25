#pragma once

#include <stddef.h>

#include "config.h"
#include "types.h"

typedef struct {
	unsigned int step_number; // Number of steps elapsed

	unsigned int N; // Number of bodies in the simulation

	void* mem; // Shared memory for the simulation

	real* mas; //             Mass (  Scalar)
	real* rad; // Smoothing Radius (  Scalar)
	real* pos; //         Position (3-Vector)
	real* vel; //         Velocity (3-Vector)

	real* acc; // Acceleration (3-Vector)

	real* pot; // Gravitational Potential (Scalar)

	real* rho; //  Density (Scalar)
	real* prs; // Pressure (Scalar)

	real* ken; //   Kinetic Energy (Scalar)
	real* pen; // Potential Energy (Scalar)
} sim_t;

// Initialize a simulation
int simulation_init(sim_t* restrict simulation, unsigned int N);

// Free a simulation
int simulation_free(sim_t* restrict simulation);

// Read a simulation state from a file
int simulation_read(sim_t* restrict simulation, const char* restrict filename);

// Save a simulation state to a file
int simulation_save(const sim_t* restrict simulation, const char* restrict filename);

// Updates a simulation (runs a timestep)
int simulation_step(sim_t* restrict simulation);
