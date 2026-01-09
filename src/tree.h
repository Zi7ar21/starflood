#pragma once

#include "simulation.h"
#include "types.h"

enum tree_conf {
	TREE_MAS_X,
	TREE_MAS_Y,
	TREE_MAS_Z,
	TREE_MAS_W,
	TREE_DOF
};

typedef struct {
	int child[8];

	real param[TREE_DOF];
} node_t;

typedef struct {
	unsigned int max_depth;
	unsigned int max_nodes;
	unsigned int num_nodes;

	real bounds_min[3];
	real bounds_max[3];

	node_t* node;
} tree_t;

int tree_init(tree_t* restrict tree_ptr, unsigned int max_nodes, const real* restrict bounds_min, const real* restrict bounds_max);

int tree_free(tree_t* restrict tree_ptr);

int tree_build(tree_t* restrict tree_ptr, const sim_t* restrict sim_ptr);
