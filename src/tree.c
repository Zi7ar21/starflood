// Needed for posix_memalign() and timing.h
#define _POSIX_C_SOURCE 200112L

#include "tree.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "timing.h"
#include "types.h"

int tree_init(tree_t* restrict tree_ptr, uint max_nodes) {
	{
		tree_ptr->max_nodes = (uint)0u;
		tree_ptr->num_nodes = (uint)0u;

		for(int i = 0; i < 6; i++) {
			tree_ptr->bounds[i] = (real)0;
		}

		tree_ptr->node = (node_t*)NULL;
	}

	const size_t mem_size = sizeof(node_t) * max_nodes;

	void* mem = NULL;

	#ifdef STARFLOOD_ALIGNMENT
	posix_memalign(&mem, (size_t)STARFLOOD_ALIGNMENT, mem_size);
	#else
	mem = malloc(mem_size);
	#endif

	if(NULL == mem) {
		#ifdef STARFLOOD_ALIGNMENT
		fprintf(stderr, "%s error: mem is NULL after posix_memalign(&mem, %zu, %zu", "tree_init()", (size_t)STARFLOOD_ALIGNMENT, mem_size);
		#else
		fprintf(stderr, "%s error: mem is NULL after malloc(%zu", "tree_init()", mem_size);
		#endif
		perror(")");
		return STARFLOOD_FAILURE;
	}

	memset(mem, 0, sizeof(node_t) * tree_ptr->max_nodes);
	
	{
		tree_ptr->max_nodes = max_nodes;

		tree_ptr->node = (node_t*)mem;
	}

	return STARFLOOD_SUCCESS;
}

int tree_free(tree_t* restrict tree_ptr) {
	free(tree_ptr->node);

	{
		tree_ptr->max_nodes = (uint)0u;
		tree_ptr->num_nodes = (uint)0u;

		for(int i = 0; i < 6; i++) {
			tree_ptr->bounds[i] = (real)0;
		}

		tree_ptr->node = (node_t*)NULL;
	}

	return STARFLOOD_SUCCESS;
}

int tree_build(tree_t* restrict tree_ptr, unsigned int num_parts, const real* restrict pos, const real* restrict mas) {
	TIMING_INIT();

	tree_t tree = *tree_ptr;

	tree.num_nodes = 1u;

	node_t* node = tree.node;

	TIMING_START();

	// Purge the tree
	{
		node_t node_init;

		node_init.parent = (uint)0u;
		node_init.bodies = (uint)0u;

		for(int i = 0; i < 8; i++) {
			node_init.branch[i] = (uint)0;
		}

		for(int i = 0; i < 3; i++) {
			node_init.bounds[i] = (real)0;
		}

		for(int i = 0; i < (int)TREE_DOF; i++) {
			node_init.params[i] = (real)0;
		}

		for(unsigned int i = 0u; i < (unsigned int)TREE_NODES_MAX; i++) {
			node[i] = node_init;
		}
	}

	TIMING_STOP();
	TIMING_PRINT("tree_build()", "purge");

	// set the root node's bounds
	for(int i = 0; i < 6; i++) {
		node[0].bounds[i] = tree.bounds[i];
	}

	TIMING_START();

	// Insert all of the particles into the tree
	for(unsigned int idx = 0u; idx < num_parts; idx++) {
		// particle being inserted's position
		const real part_pos[3] = {
			pos[3u*idx+0u],
			pos[3u*idx+1u],
			pos[3u*idx+2u]
		};

		// particle being inserted's mass
		const real part_mas = mas[idx];

		// Current node bounds
		real bounds_min[3];
		real bounds_max[3];

		// Initialize the bounds tracking variables
		for(int i = 0; i < 3; i++) {
			bounds_min[i] = tree.bounds[3*0+i];
			bounds_max[i] = tree.bounds[3*1+i];
		}

		// Check if the particle is within tree bounds
		if(bounds_min[0] > part_pos[0] || part_pos[0] > bounds_max[0]
		|| bounds_min[1] > part_pos[1] || part_pos[1] > bounds_max[1]
		|| bounds_min[2] > part_pos[2] || part_pos[2] > bounds_max[2]) {
			continue;
		}

		// Current node
		uint cur_node = (uint)0u;

		// It takes at most TREE_DEPTH_MAX iterations to find a place to insert the particle
		for(int depth = 0; depth <= TREE_DEPTH_MAX; depth++) {
			// Number of particles in the current node
			const uint cur_bodies = node[cur_node].bodies;

			node[cur_node].bodies += (uint)1u;

			// Midpoint of the currrent node
			const real bounds_mid[3] = {
				(real)0.5 * (bounds_min[0] + bounds_max[0]),
				(real)0.5 * (bounds_min[1] + bounds_max[1]),
				(real)0.5 * (bounds_min[2] + bounds_max[2])
			};

			// Add the particle's mass to the current node
			node[cur_node].params[TREE_MAS_X] += part_mas * part_pos[0];
			node[cur_node].params[TREE_MAS_Y] += part_mas * part_pos[1];
			node[cur_node].params[TREE_MAS_Z] += part_mas * part_pos[2];
			node[cur_node].params[TREE_MAS_W] += part_mas;

			// If the current node has enough room, insert the particle
			if( (uint)8u > cur_bodies ) {
				node[cur_node].branch[cur_bodies] = (uint)idx;
				break;
			}

			// If the current node is full, split it
			if( (uint)8u == cur_bodies ) {
				// Check if splitting the node would exceed the maximum number of tree nodes
				if( tree.max_nodes <= (tree.num_nodes + (uint)8u) ) {
					fprintf(stderr, "%s error: Exceeded maximum number of tree nodes!\n", "tree_build()");
					return STARFLOOD_FAILURE;
				}

				// Check if splitting the node would exceed the maximum tree depth
				if( TREE_DEPTH_MAX <= depth ) {
					fprintf(stderr, "%s error: Exceeded maximum tree depth!\n", "tree_build()");
					return STARFLOOD_FAILURE;
				}

				node_t new_node[8];

				// Initialize new_node
				{
					node_t nul_node;

					nul_node.parent = cur_node;
					nul_node.bodies = (uint)0u;

					for(int i = 0; i < 8; i++) {
						nul_node.branch[i] = (uint)0u;
					}

					for(int i = 0; i < 6; i++) {
						nul_node.bounds[i] = (real)0;
					}

					nul_node.params[TREE_MAS_X] = (real)0.0;
					nul_node.params[TREE_MAS_Y] = (real)0.0;
					nul_node.params[TREE_MAS_Z] = (real)0.0;
					nul_node.params[TREE_MAS_W] = (real)0.0;

					for(int i = 0; i < 8; i++) {
						new_node[i] = nul_node;
					}
				}

				// Loop over each octant
				for(int i = 0; i < 8; i++) {
					// Loop over each dimension
					for(int j = 0; j < 3; j++) {
						// Set the new node's bounds
						new_node[i].bounds[3*0+j] = (i >> j) & 1 ? bounds_mid[j] : bounds_min[j];
						new_node[i].bounds[3*1+j] = (i >> j) & 1 ? bounds_max[j] : bounds_mid[j];
					}
				}

				// Split the current node's children into the new nodes
				for(int i = 0; i < 8; i++) {
					// Child particle's index
					const uint child_idx = node[cur_node].branch[i];

					// Child particle's position
					const real child_pos[3] = {
						pos[3u*child_idx+0u],
						pos[3u*child_idx+1u],
						pos[3u*child_idx+2u]
					};

					// Child particle's mass
					const real child_mas = mas[child_idx];

					// Octant the child particle should be inserted into
					const int child_octant =
						(bounds_mid[0] < child_pos[0] ? 1 : 0) +
						(bounds_mid[1] < child_pos[1] ? 2 : 0) +
						(bounds_mid[2] < child_pos[2] ? 4 : 0);

					// Number of particles in the child node
					const uint child_bodies = new_node[child_octant].bodies;

					new_node[child_octant].bodies += (uint)1u;

					new_node[child_octant].branch[child_bodies] = child_idx;

					new_node[child_octant].params[TREE_MAS_X] += child_mas * child_pos[0];
					new_node[child_octant].params[TREE_MAS_Y] += child_mas * child_pos[1];
					new_node[child_octant].params[TREE_MAS_Z] += child_mas * child_pos[2];
					new_node[child_octant].params[TREE_MAS_W] += child_mas;
				}

				// Add the new nodes to the tree
				{
					for(int i = 0; i < 8; i++) {
						node[tree.num_nodes + (uint)i] = new_node[i];
					}

					for(int i = 0; i < 8; i++) {
						node[cur_node].branch[i] = tree.num_nodes + (uint)i;
					}

					tree.num_nodes += (uint)8u;
				}
			}

			// Octant the particle should be inserted into
			//
			// Less than equal to is used here, so theoretically
			// the tree can have 9 nodes in the same position
			// without exceeding the max depth (unlikely if the
			// initial conditions were generated correctly)
			int octant =
				(bounds_mid[0] <= part_pos[0] ? 1 : 0) +
				(bounds_mid[1] <= part_pos[1] ? 2 : 0) +
				(bounds_mid[2] <= part_pos[2] ? 4 : 0);

			// Set the bounds
			for(int i = 0; i < 3; i++) {
				bounds_min[i] = (octant >> i) & 1 ? bounds_mid[i] : bounds_min[i];
				bounds_max[i] = (octant >> i) & 1 ? bounds_max[i] : bounds_mid[i];
			}

			// Set the current node to the child node the particle is to be inserted into
			cur_node = node[cur_node].branch[octant];
		}
	}

	TIMING_STOP();
	TIMING_PRINT("tree_build()", "insertion");
	TIMING_START();

	for(unsigned int i = 0u; i < (unsigned int)tree.num_nodes; i++) {
		real mass = node[i].params[TREE_MAS_W];

		real inv_mass = (real)0.0 < mass ? (real)1.0 / mass : (real)1.0;

		node[i].params[TREE_MAS_X] *= inv_mass;
		node[i].params[TREE_MAS_Y] *= inv_mass;
		node[i].params[TREE_MAS_Z] *= inv_mass;
	}

	TIMING_STOP();
	TIMING_PRINT("tree_build()", "center_of_mass");

	printf("  built tree with %u nodes\n", tree.num_nodes);

	return STARFLOOD_SUCCESS;
}
