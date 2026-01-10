#pragma once

#include "types.h"

enum tree_conf {
	TREE_MAS_X,
	TREE_MAS_Y,
	TREE_MAS_Z,
	TREE_MAS_W,
	TREE_DOF
};

typedef struct {
	i32 parent; // parent of this node (-1: this node is, -2: this node is the tree root)
	i32 bodies; // number of bodies in this node
	i32 child[8]; // children of this node (8 >= bodies: children are particles, 8 < bodies: children are nodes)
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

int tree_init(tree_t* restrict tree_ptr, unsigned int max_nodes);

int tree_free(tree_t* restrict tree_ptr);

int tree_build(tree_t* restrict tree_ptr, unsigned int N, const real* restrict pos, const real* restrict mas);
