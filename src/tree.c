// Needed for timing.h
#define _POSIX_C_SOURCE 199309L

#include "tree.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "timing.h"
#include "types.h"

int tree_init(tree_t* restrict tree_ptr, unsigned int max_nodes) {
	tree_t tree = *tree_ptr;

	tree.max_depth = 24u;
	tree.max_nodes = max_nodes;
	tree.num_nodes = 1u;

	for(int i = 0; i < 3; i++) {
		tree.bounds_min[i] = (real)0.0;
		tree.bounds_max[i] = (real)0.0;
	}

	size_t tree_size = sizeof(node_t) * tree.max_nodes;

	void* mem = malloc(tree_size);

	if(NULL == mem) {
		return STARFLOOD_FAILURE;
	}

	memset(mem, 0, tree_size);

	tree.node = (node_t*)mem;

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

int tree_build(tree_t* restrict tree_ptr, unsigned int N, const real* restrict pos) {
	TIMING_INIT();

	tree_t tree = *tree_ptr;

	tree.num_nodes = 1u;

	node_t* node = tree.node;

	TIMING_START();

	// Purge the tree
	{
		node_t node_init;

		node_init.parent = 0;
		node_init.bodies = 0;

		for(int i = 0; i < 8; i++) {
			node_init.child[i] = -1;
		}

		for(int i = 0; i < (int)TREE_DOF; i++) {
			node_init.param[i] = (real)0.0;
		}

		for(unsigned int i = 0u; i < N; i++) {
			node[i] = node_init;
		}
	}
	
	TIMING_STOP();
	TIMING_PRINT("tree_build()", "purge");
	TIMING_START();

	for(unsigned int i = 0u; i < N; i++) {
		real body_pos[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		unsigned int cur_depth = 0u;

		i32 curr_node = (i32)( 0);
		i32 prev_node = (i32)(-1);

		real bounds_min[3];
		real bounds_max[3];

		for(unsigned int j = 0u; j < 3u; j++) {
			bounds_min[j] = tree.bounds_min[j];
			bounds_max[j] = tree.bounds_max[j];
		}

		for(unsigned int depth = 0u; depth < tree.max_depth; depth++) {
			int bodies = node[curr_node].bodies;

			real bounds_mid[3];

			for(unsigned int j = 0u; j < 3u; j++) {
				bounds_mid[j] = (real)0.5 * (bounds_min[j] + bounds_max[j]);
			}

			// This node has room
			if(8 > bodies) {
				node[curr_node].child[bodies] = (int)i;
				node[curr_node].bodies += 1;
				break;
			}

			// This node is full, split it
			if(8 == bodies) {
				if( tree.max_nodes < (tree.num_nodes + 8u) ) {
					fprintf(stderr, "%s error: Exceeded maximum number of tree nodes!\n", "tree_build()");
					return STARFLOOD_FAILURE;
				}

				// Children bodies of this current node
				int child[8];

				for(int child_num = 0; child_num < 8; child_num++) {
					child[child_num] = node[curr_node].child[child_num]; // Save the body ID
					node[curr_node].child[child_num] = (int)tree.num_nodes + child_num; // Set the node's child ID to the new node's ID
				}

				for(int child_num = 0; child_num < 8; child_num++) {
					unsigned int id = (unsigned int)child[child_num];

					real child_body_pos[3] = {
						pos[3u*id+0u],
						pos[3u*id+1u],
						pos[3u*id+2u]
					};

					int child_octant =
						(bounds_mid[0] < child_body_pos[0] ? 1 : 0) +
						(bounds_mid[1] < child_body_pos[1] ? 2 : 0) +
						(bounds_mid[2] < child_body_pos[2] ? 4 : 0);

					// Get the ID of the new node
					int new_node_id = (int)tree.num_nodes + child_octant;

					// Get the number of bodies in the new node
					int new_node_bodies = node[new_node_id].bodies;

					// Add the child body ID from the current node to the new node
					node[new_node_id].child[new_node_bodies] = child[child_num];
					node[new_node_id].bodies = new_node_bodies;
				}

				node[curr_node].bodies = 8 + 1; // we will change this later

				tree.num_nodes += 8u;
			}

			int octant =
				(bounds_mid[0] < body_pos[0] ? 1 : 0) +
				(bounds_mid[1] < body_pos[1] ? 2 : 0) +
				(bounds_mid[2] < body_pos[2] ? 4 : 0);

			for(int j = 0; j < 3; j++) {
				bounds_min[j] = bounds_mid[j] <  body_pos[j] ? bounds_mid[j] : bounds_min[j];
				bounds_max[j] = bounds_mid[j] >= body_pos[j] ? bounds_mid[j] : bounds_max[j];
			}

			curr_node = node[curr_node].child[octant];

			cur_depth = depth + 1;
		}

		if(tree.max_depth == cur_depth) {
			fprintf(stderr, "%s error: Exceeded maximum tree depth!\n", "tree_build()");
			return STARFLOOD_FAILURE;
		}
	}

	TIMING_STOP();
	TIMING_PRINT("tree_build()", "insertion");

	printf("  %u\n", tree.num_nodes);

	return STARFLOOD_SUCCESS;
}
