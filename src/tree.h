#pragma once

#include "types.h"

// Enumeration of tree parameters (configuration space)
// for Barnes-Hut, this is the mass and center of mass
enum tree_conf {
	TREE_MAS_X, // monopole offset (X-component of center of mass, relative to node center)
	TREE_MAS_Y, // monopole offset (Y-component of center of mass, relative to node center)
	TREE_MAS_Z, // monopole offset (Z-component of center of pass, relative to node center)
	TREE_MAS_W, // monopole weight (mass within the node)
	TREE_DOF // degrees of freedom
};

// Structure containing data of a tree node
typedef struct {
	uint parent; // parent of this node
	uint bodies; // number of particles in this node and children

	uint child[8]; // children of this node (8 >= bodies: children are particles, 8 < bodies: children are nodes)

	real bounds_min[3]; // node bounds min
	real bounds_max[3]; // node bounds max

	real param[TREE_DOF]; // tree parameters
} node_t;

// Structure containing data of a tree
typedef struct {
	uint max_nodes; // max number of nodes in the tree
	uint num_nodes; //     number of nodes in the tree

	real bounds_min[3]; // tree bounds min
	real bounds_max[3]; // tree bounds max

	node_t* node; // dynamically allocated array containing the tree nodes
} tree_t;

// Initializes a tree
int tree_init(tree_t* restrict tree_ptr, unsigned int max_nodes);

// Frees resources allocated by a tree
int tree_free(tree_t* restrict tree_ptr);

// Builds a tree for calculations
int tree_build(tree_t* restrict tree_ptr, unsigned int num_parts, const real* restrict pos, const real* restrict mas);
