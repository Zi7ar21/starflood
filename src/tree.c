#include "tree.h"

#include <stdlib.h>

#include "common.h"
#include "simulation.h"

int tree_init(tree_t* restrict tree_ptr, unsigned int max_nodes, const real* restrict bounds_min, const real* restrict bounds_max) {
	tree_t tree = *tree_ptr;

	tree.max_depth = 24u;
	tree.max_nodes = max_nodes;
	tree.num_nodes = 1u;

	for(int i = 0; i < 3; i++) {
		tree.bounds_min[i] = bounds_min[i];
		tree.bounds_max[i] = bounds_max[i];
	}

	tree.node = (node_t*)malloc(sizeof(node_t) * tree.max_nodes);

	if(NULL == (void*)tree.node) {
		return STARFLOOD_FAILURE;
	}

	node_t null_node;

	for(unsigned int j = 0u; j < 8u; j++) {
		null_node.child[j] = -1;
	}

	for(unsigned int j = 0u; j < (unsigned int)TREE_DOF; j++) {
		null_node.param[j] = (real)0.0;
	}

	for(unsigned int i = 0u; i < tree.max_nodes; i++) {
		tree.node[i] = null_node;
	}

	*tree_ptr = tree;

	return STARFLOOD_SUCCESS;
}

int tree_free(tree_t* restrict tree_ptr) {
	tree_t tree = *tree_ptr;

	{
		free(tree.node);
	}

	*tree_ptr = tree;

	return STARFLOOD_SUCCESS;
}

int tree_build(tree_t* restrict tree_ptr, const sim_t* restrict sim_ptr) {
	return STARFLOOD_FAILURE;
}
