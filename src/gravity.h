#pragma once

#include "simulation.h"
#include "types.h"

#ifdef ENABLE_GRID
// Solve gravity using the mesh (particle-particle/particle-mesh method)
int solve_gravity_part_mesh(sim_t* restrict sim_ptr);
#endif

#ifdef ENABLE_TREE
// Solve gravity using the tree (particle-particle/particle-tree method)
int solve_gravity_part_tree(sim_t* restrict sim_ptr);
#endif

// Solve gravity using O(N^2) pairwise summation/brute-force (particle-particle method)
int solve_gravity_part_part(sim_t* restrict sim_ptr);
