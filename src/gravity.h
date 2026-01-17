#pragma once

#include "config.h"
#include "simulation.h"
#include "types.h"

#ifdef ENABLE_GRID
// Solve gravity using the mesh (particle-particle/particle-mesh method)
int solve_gravity_part_mesh(const sim_t* restrict sim_ptr);
#endif

#ifdef ENABLE_TREE
// Solve gravity using the tree (particle-particle/particle-tree method)
#ifdef SOLVER_VIS
int solve_gravity_part_tree(const sim_t* restrict sim_ptr, f32* render_buffer, unsigned int sizex, unsigned int sizey);
#else
int solve_gravity_part_tree(const sim_t* restrict sim_ptr);
#endif
#endif

// Solve gravity using O(N^2) pairwise summation/brute-force (particle-particle method)
int solve_gravity_part_part(const sim_t* restrict sim_ptr);
