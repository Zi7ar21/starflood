// Needed for timing.h
#define _POSIX_C_SOURCE 199309L

#include "gravity.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "config.h"
#include "rng.h"
#include "timing.h"
#include "types.h"

#ifdef SOLVER_VIS
extern real matMVP[16];

int rasterize_line_3d(f32* restrict render_buffer, const unsigned int sizex, const unsigned int sizey, const f32* restrict rgba, const real* restrict xyz0, const real* restrict xyz1) {
	real p0[4] = {xyz0[0], xyz0[1], xyz0[2], (real)1};
	real p1[4] = {xyz1[0], xyz1[1], xyz1[2], (real)1};

	// p0 = matMVP * p0;
	{
		real y[4];

		for(int i = 0; i < 4; i++) {
			y[i] =
				(matMVP[4*i+0] * p0[0]) +
				(matMVP[4*i+1] * p0[1]) +
				(matMVP[4*i+2] * p0[2]) +
				(matMVP[4*i+3] * p0[3]);
		}

		for(int i = 0; i < 4; i++) {
			p0[i] = y[i];
		}
	}

	// p1 = matMVP * p1;
	{
		real y[4];

		for(int i = 0; i < 4; i++) {
			y[i] =
				(matMVP[4*i+0] * p1[0]) +
				(matMVP[4*i+1] * p1[1]) +
				(matMVP[4*i+2] * p1[2]) +
				(matMVP[4*i+3] * p1[3]);
		}

		for(int i = 0; i < 4; i++) {
			p1[i] = y[i];
		}
	}

	// Adapted from section 1.7 "Program to plot a line" of
	// "A Rasterizing Algorithm for Drawing Curves" by Alois Zingl
	// https://zingl.github.io/Bresenham.pdf

	const int x0 = (int)( (real)(-0.5) + (real)0.5 * (real)sizey * p0[0] + (real)0.5 * (real)sizex );
	const int y0 = (int)( (real)(-0.5) + (real)0.5 * (real)sizey * p0[1] + (real)0.5 * (real)sizey );

	const int x1 = (int)( (real)(-0.5) + (real)0.5 * (real)sizey * p1[0] + (real)0.5 * (real)sizex );
	const int y1 = (int)( (real)(-0.5) + (real)0.5 * (real)sizey * p1[1] + (real)0.5 * (real)sizey );

	#if 0
	printf(
		"  p0: (% f, % f, % f)\n"
		"  p1: (% f, % f, % f)\n"
		"  c0: (% d, % d)\n"
		"  c1: (% d, % d)\n",
		p0[0], p0[1], p0[2],
		p1[0], p1[1], p1[2],
		x0, y0,
		x1, y1
	);
	#endif

	int x = x0, y = y0;

	const int dx =  abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
	const int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;

	int err = dx + dy;

	for(;;) {
		if(0 < x && x < (int)sizex && 0 < y && y < (int)sizey) {
			for(unsigned int chan = 0u; chan < 4u; chan++) {
				render_buffer[4u*(unsigned int)((int)sizex*y+x)+chan] += rgba[chan];
			}
		}

		const int e2 = 2 * err;

		if(e2 >= dy) {
			if(x == x1) { break; }
			err += dy;
			x += sx;
		}

		if(e2 <= dx) {
			if(y == y1) { break; }
			err += dx;
			y += sy;
		}
	}

	return STARFLOOD_SUCCESS;
}
#endif

#if 0
//#ifdef ENABLE_GRID
// Solve gravity using the grid (particle-particle/particle-mesh)
int solve_gravity_part_mesh(sim_t* restrict sim_ptr) {
	TIMING_INIT();

	sim_t sim = *sim_ptr;

	grid_t grid = sim.grid;

	unsigned int grid_dim[3] = {grid.dim[0], grid.dim[1], grid.dim[2]};

	real bounds_min[3] = {grid.bounds_min[0], grid.bounds_min[1], grid.bounds_min[2]};
	real bounds_max[3] = {grid.bounds_max[0], grid.bounds_max[1], grid.bounds_max[2]};

	real bounds_wid[3] = {
		bounds_max[0] - bounds_min[0],
		bounds_max[1] - bounds_min[1],
		bounds_max[2] - bounds_min[2]
	};

	f32* grid_poten = (f32*)grid.mem;
	i32* grid_accum = (i32*)&grid_poten[grid_dim[2] * grid_dim[1] * grid_dim[0]];

	TIMING_START();

	for(unsigned int i = 0u; i < grid.dim[2] * grid.dim[1] * grid.dim[0]; i++) {
		#ifdef _OPENMP
		#pragma omp atomic write
		#endif
		grid_accum[i] = (i32)0;
	}

	TIMING_STOP();
	TIMING_PRINT("solve_gravity()", "grid_clear");
	TIMING_START();

	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic, 256)
	#endif
	for(unsigned int i = 0u; i < N; i++) {
		real pos_world[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		real pos_local[3] = {
			(pos_world[0] - bounds_min[0]) / bounds_wid[0],
			(pos_world[1] - bounds_min[1]) / bounds_wid[1],
			(pos_world[2] - bounds_min[2]) / bounds_wid[2]
		};

		for(unsigned int j = 0u; j < (unsigned int)GRID_SAMPLES; j++) {
			u32 s[4] = {
				(u32)0xB79ABC95u + (u32)j,
				(u32)0xE0DA3F84u + (u32)i,
				(u32)0xAB75F07Bu + (u32)step_number,
				(u32)0xCF52CA01u
			};

			pcg4d(s);
			pcg4d(s); // second round for better statistical quality

			real r[4] = {
				(real)INV_PCG32_MAX * (real)s[0],
				(real)INV_PCG32_MAX * (real)s[1],
				(real)INV_PCG32_MAX * (real)s[2],
				(real)INV_PCG32_MAX * (real)s[3]
			};

			// Box-Muller Transform
			// https://en.wikipedia.org/wiki/Box–Muller_transform
			real n[4] = {
				real_sqrt( (real)(-2.0) * real_log(r[0]) ) * real_cos( (real)TAU * r[1] ),
				real_sqrt( (real)(-2.0) * real_log(r[0]) ) * real_sin( (real)TAU * r[1] ),
				real_sqrt( (real)(-2.0) * real_log(r[2]) ) * real_cos( (real)TAU * r[3] ),
				real_sqrt( (real)(-2.0) * real_log(r[2]) ) * real_sin( (real)TAU * r[3] )
			};

			int coord[3] = {
				(int)( (real)0.500 * n[0] + (real)grid_dim[0] * pos_local[0] + (real)0.5),
				(int)( (real)0.500 * n[1] + (real)grid_dim[1] * pos_local[1] + (real)0.5),
				(int)( (real)0.500 * n[2] + (real)grid_dim[2] * pos_local[2] + (real)0.5)
			};

			if(0 <= coord[0] && coord[0] < (int)grid_dim[0]
			&& 0 <= coord[1] && coord[1] < (int)grid_dim[1]
			&& 0 <= coord[2] && coord[2] < (int)grid_dim[2]) {
				unsigned int index = grid_dim[1]*grid_dim[0]*(unsigned int)coord[2]+grid_dim[0]*(unsigned int)coord[1]+(unsigned int)coord[0];

				#ifdef _OPENMP
				#pragma omp atomic update
				#endif
				grid_accum[index] += (i32)1;
			}
		}
	}

	TIMING_STOP();
	TIMING_PRINT("solve_gravity()", "grid_dist");
	TIMING_START();

	double inv_grid_samples = 1.0 / (double)GRID_SAMPLES;

	for(unsigned int i = 0u; i < grid.dim[2u] * grid.dim[1u] * grid.dim[0u]; i++) {
		i32 val = (i32)0;

		#ifdef _OPENMP
		#pragma omp atomic read
		#endif
		val = grid_accum[i];

		grid_poten[i] = (f32)(inv_grid_samples * (double)val);
	}

	TIMING_STOP();
	TIMING_PRINT("solve_gravity()", "grid_read");
	TIMING_START();

	#ifdef ENABLE_FFT
	grid_solve_fft(&grid);
	#endif

	TIMING_STOP();
	TIMING_PRINT("solve_gravity()", "solve_fft");
}
#endif

#ifdef ENABLE_TREE
// Solve gravity using the tree (particle-particle/particle-tree method)
#ifdef SOLVER_VIS
int solve_gravity_part_tree(const sim_t* restrict sim_ptr, f32* render_buffer, unsigned int sizex, unsigned int sizey)
#else
int solve_gravity_part_tree(const sim_t* restrict sim_ptr)
#endif
{
	const sim_t sim = *sim_ptr;

	const unsigned int N = sim.N;
	const real* pos = sim_find(&sim, SIM_POS);
	const real* vel = sim_find(&sim, SIM_VEL);
	real* acc = sim_find(&sim, SIM_ACC);
	const real* mas = sim_find(&sim, SIM_MAS);
	real* pot = sim_find(&sim, SIM_POT);

	const tree_t tree = sim.tree;

	const uint num_nodes = tree.num_nodes;

	const node_t* node = tree.node;

	#ifdef _OPENMP
		#ifdef ENABLE_OFFLOAD_SIM
		// TODO: fix
		//#pragma omp target teams distribute parallel for map(tofrom: acc[:3u*N], pot[:N]) map(to: pos[:3u*N], mas[:N], node[:num_nodes])
		#else
		#pragma omp parallel for schedule(dynamic, 256)
		#endif
	#endif
	for(unsigned int idx = 0u; idx < N; idx++) {
		real acc_sum[3] = {(real)0.0, (real)0.0, (real)0.0};

		#ifdef ENABLE_KAHAN_SUMMATION_SOLVER
		real acc_com[3] = {(real)0.0, (real)0.0, (real)0.0};
		#endif

		real pot_sum = (real)0.0;

		#ifdef ENABLE_KAHAN_SUMMATION_ENERGY
		real pot_com = (real)0.0;
		#endif

		// body i's position
		const real r_i[3] = {
			pos[3u*idx+0u],
			pos[3u*idx+1u],
			pos[3u*idx+2u]
		};

		const real m_i = mas[idx]; // body i's mass

		uint next_node = (uint)0u; // next node to be searched

		int cur_depth = 0; // current depth of the search
		int tree_path[TREE_DEPTH_MAX]; // octant of the node at cur_depth we are searching

		for(int i = 0; i < TREE_DEPTH_MAX; i++) {
			tree_path[i] = -1;
		}

		//for(unsigned int j = 0u; j < 8u * (unsigned int)TREE_NODES_MAX; j++) {
		//while(1) {
		for(unsigned int search_iter = 0u; search_iter < 2u * 8u * (unsigned int)TREE_NODES_MAX; search_iter++) {
			//if(i == j) {
			//	continue;
			//}

			// check if we are done searching
			if(0 > cur_depth) {
				//printf( (N - 1) < idx ? "%u\n" : "%u ", search_iter );
				break;
			}

			const node_t this_node = node[next_node];

			// this node's bounds
			const real bounds_min[3] = {this_node.bounds[3*0+0], this_node.bounds[3*0+1], this_node.bounds[3*0+2]};
			const real bounds_max[3] = {this_node.bounds[3*1+0], this_node.bounds[3*1+1], this_node.bounds[3*1+2]};

			// check if this node is empty
			if( (uint)0u >= this_node.bodies ) {
				next_node = this_node.parent;
				cur_depth = cur_depth - 1;
				continue;
			}

			// check if this node has children (nodes)
			if( (uint)8u < this_node.bodies ) {
				// current octant being searched
				const int octant = tree_path[cur_depth];

				// should not be possible
				//if(8 <= octant) {
				//	fprintf(stderr, "%s error: glitch\n", "solve_gravity_part_tree()");
				//	break;
				//	//return STARFLOOD_FAILURE;
				//}

				// check if this node is done being searched
				if(7 <= octant) {
					next_node = this_node.parent;
					tree_path[cur_depth] = -1;
					cur_depth = cur_depth - 1;
					continue;
				}

				// check if we haven't started searching this node
				if(0 > octant) {
					// this node's midpoint
					const real bounds_mid[3] = {
						(real)0.5 * (bounds_min[0] + bounds_max[0]),
						(real)0.5 * (bounds_min[1] + bounds_max[1]),
						(real)0.5 * (bounds_min[2] + bounds_max[2])
					};

					// vector to this node's midpoint from the current particle
					const real bounds_off[3] = {
						bounds_mid[0] - r_i[0],
						bounds_mid[1] - r_i[1],
						bounds_mid[2] - r_i[2]
					};

					// distance squared to this node's midpoint from the current particle
					const real node_dist2 = (bounds_off[0]*bounds_off[0])+(bounds_off[1]*bounds_off[1])+(bounds_off[2]*bounds_off[2]);

					// this node's width
					const real node_width = (real)(1.0 / 3.0) * (
						(bounds_max[0] - bounds_min[0]) +
						(bounds_max[1] - bounds_min[1]) +
						(bounds_max[2] - bounds_min[2])
					);

					// width / distance
					const real dist_ratio = (real)0.0 < node_dist2 ? (node_width * node_width) / node_dist2 : (real)0.0;

					// draw bounding box
					#ifdef SOLVER_VIS
					if( (unsigned int)SOLVER_VIS == idx ) {
						const f32 rgba[4] = {
							( (f32)1.0 / (f32)(cur_depth + 1) ) * (0 == (cur_depth % 3) ? (f32)0.800 : (f32)0.100),
							( (f32)1.0 / (f32)(cur_depth + 1) ) * (1 == (cur_depth % 3) ? (f32)0.800 : (f32)0.100),
							( (f32)1.0 / (f32)(cur_depth + 1) ) * (2 == (cur_depth % 3) ? (f32)0.800 : (f32)0.100),
							(f32)1.000
						};

						{
							// axis 0
							{
								const real p0[3] = {bounds_min[0], bounds_min[1], bounds_min[2]}, p1[3] = {bounds_max[0], bounds_min[1], bounds_min[2]};
								rasterize_line_3d(render_buffer, sizex, sizey, rgba, p0, p1);
							}
							{
								const real p0[3] = {bounds_min[0], bounds_max[1], bounds_min[2]}, p1[3] = {bounds_max[0], bounds_max[1], bounds_min[2]};
								rasterize_line_3d(render_buffer, sizex, sizey, rgba, p0, p1);
							}
							{
								const real p0[3] = {bounds_min[0], bounds_min[1], bounds_max[2]}, p1[3] = {bounds_max[0], bounds_min[1], bounds_max[2]};
								rasterize_line_3d(render_buffer, sizex, sizey, rgba, p0, p1);
							}
							{
								const real p0[3] = {bounds_min[0], bounds_max[1], bounds_max[2]}, p1[3] = {bounds_max[0], bounds_max[1], bounds_max[2]};
								rasterize_line_3d(render_buffer, sizex, sizey, rgba, p0, p1);
							}

							// axis 1
							{
								const real p0[3] = {bounds_min[0], bounds_min[1], bounds_min[2]}, p1[3] = {bounds_min[0], bounds_max[1], bounds_min[2]};
								rasterize_line_3d(render_buffer, sizex, sizey, rgba, p0, p1);
							}
							{
								const real p0[3] = {bounds_max[0], bounds_min[1], bounds_min[2]}, p1[3] = {bounds_max[0], bounds_max[1], bounds_min[2]};
								rasterize_line_3d(render_buffer, sizex, sizey, rgba, p0, p1);
							}
							{
								const real p0[3] = {bounds_min[0], bounds_min[1], bounds_max[2]}, p1[3] = {bounds_min[0], bounds_max[1], bounds_max[2]};
								rasterize_line_3d(render_buffer, sizex, sizey, rgba, p0, p1);
							}
							{
								const real p0[3] = {bounds_max[0], bounds_min[1], bounds_max[2]}, p1[3] = {bounds_max[0], bounds_max[1], bounds_max[2]};
								rasterize_line_3d(render_buffer, sizex, sizey, rgba, p0, p1);
							}

							// axis 2
							{
								const real p0[3] = {bounds_min[0], bounds_min[1], bounds_min[2]}, p1[3] = {bounds_min[0], bounds_min[1], bounds_max[2]};
								rasterize_line_3d(render_buffer, sizex, sizey, rgba, p0, p1);
							}
							{
								const real p0[3] = {bounds_max[0], bounds_min[1], bounds_min[2]}, p1[3] = {bounds_max[0], bounds_min[1], bounds_max[2]};
								rasterize_line_3d(render_buffer, sizex, sizey, rgba, p0, p1);
							}
							{
								const real p0[3] = {bounds_min[0], bounds_max[1], bounds_min[2]}, p1[3] = {bounds_min[0], bounds_max[1], bounds_max[2]};
								rasterize_line_3d(render_buffer, sizex, sizey, rgba, p0, p1);
							}
							{
								const real p0[3] = {bounds_max[0], bounds_max[1], bounds_min[2]}, p1[3] = {bounds_max[0], bounds_max[1], bounds_max[2]};
								rasterize_line_3d(render_buffer, sizex, sizey, rgba, p0, p1);
							}
						}
					}
					#endif

					// check if we are to visit this node's children
					if( (real)TREE_THETA <= dist_ratio ) {
						next_node = this_node.branch[0];
						tree_path[cur_depth] = 0;
						cur_depth = cur_depth + 1;
						continue;
					}
				} else {
					// we are searching this node, move to the next octant
					next_node = this_node.branch[octant + 1];
					tree_path[cur_depth] = octant + 1;
					cur_depth = cur_depth + 1;
					continue;
				}
			}

			// body j's position
			real r_j[3] = {
				this_node.params[TREE_MAS_X],
				this_node.params[TREE_MAS_Y],
				this_node.params[TREE_MAS_Z]
			};

			real m_j = this_node.params[TREE_MAS_W]; // body j's mass

			// factor out the current particle
			if(bounds_min[0] <= r_i[0] && r_i[0] <= bounds_max[0]
			&& bounds_min[1] <= r_i[1] && r_i[1] <= bounds_max[1]
			&& bounds_min[2] <= r_i[2] && r_i[2] <= bounds_max[2]
			&& (uint)1u < this_node.bodies) {
				for(int k = 0; k < 3; k++) {
					r_j[k] *= m_j;
				}

				for(int k = 0; k < 3; k++) {
					r_j[k] -= m_i * r_i[k];
				}

				m_j -= m_i;

				for(int k = 0; k < 3; k++) {
					r_j[k] *= (real)0.0 < m_j ? (real)1.0 / m_j : (real)1.0;
				}
			}

			// vector to body i's position from body j's position
			const real r_ij[3] = {
				r_i[0] - r_j[0],
				r_i[1] - r_j[1],
				r_i[2] - r_j[2]
			};

			// squared distance between bodies i and j
			const real r2 = (r_ij[0]*r_ij[0])+(r_ij[1]*r_ij[1])+(r_ij[2]*r_ij[2]);

			#ifdef EPSILON
			const real inv_r2 = (real)1.0 /          ( r2 + (real)(EPSILON*EPSILON) );
			const real inv_r1 = (real)1.0 / real_sqrt( r2 + (real)(EPSILON*EPSILON) );
			#else
			const real inv_r2 = (real)0.0 < r2 ? (real)1.0 /          (r2) : (real)0.0;
			const real inv_r1 = (real)0.0 < r2 ? (real)1.0 / real_sqrt(r2) : (real)0.0;
			#endif

			// gravitational potential of body j
			// (G and m_i are taken into account later)
			const real pot_j = m_j * inv_r1;

			const real F[3] = {
				pot_j * inv_r2 * r_ij[0],
				pot_j * inv_r2 * r_ij[1],
				pot_j * inv_r2 * r_ij[2]
			};

			for(unsigned int k = 0u; k < 3u; k++) {
				#ifdef ENABLE_KAHAN_SUMMATION_SOLVER
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				const real y = F[k] - acc_com[k];
				volatile real t = acc_sum[k] + y;
				volatile real z = t - acc_sum[k];
				acc_com[k] = z - y;
				acc_sum[k] = t;
				#else
				// Naïve summation
				acc_sum[k] += F[k];
				#endif
			}

			{
				#ifdef ENABLE_KAHAN_SUMMATION_ENERGY
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				const real y = pot_j - pot_com;
				volatile real t = pot_sum + y;
				volatile real z = t - pot_sum;
				pot_com = z - y;
				pot_sum = t;
				#else
				// Naïve summation
				pot_sum += pot_j;
				#endif
			}

			#ifdef SOLVER_VIS
			if( (unsigned int)SOLVER_VIS == idx ) {
				#if 0
				printf(
					"debug\n"
					"  idx: %u\n"
					"  search_iter: %u\n"
					"  next_node: %u\n"
					"  cur_depth: %d\n"
					"  r_j: (% f, % f, % f)\n"
					"  m_j: % f\n",
					idx,
					search_iter,
					next_node,
					cur_depth,
					r_j[0], r_j[1], r_j[2],
					m_j
				);
				#endif

				const f32 rgba[4] = {0.8f, 0.1f, 0.8f, 1.0f};

				rasterize_line_3d(render_buffer, sizex, sizey, rgba, r_i, r_j);
			}
			#endif

			// we are moving back up after this
			next_node = this_node.parent;
			cur_depth = cur_depth - 1;
		}

		/*
		for(unsigned int j = 0u; j < 0u; j++) {
			//if(i == j) {
			//	continue;
			//}

			real m_j = mas[j]; // body j's mass

			// body j's position
			real r_j[3] = {
				pos[3u*j+0u],
				pos[3u*j+1u],
				pos[3u*j+2u]
			};

			// vector to body i's position from body j's position
			real r_ij[3] = {
				r_i[0] - r_j[0],
				r_i[1] - r_j[1],
				r_i[2] - r_j[2]
			};

			// squared distance between bodies i and j
			real r2 = (r_ij[0]*r_ij[0])+(r_ij[1]*r_ij[1])+(r_ij[2]*r_ij[2]);

			#ifdef EPSILON
			real inv_r2 = (real)1.0 /          ( r2 + (real)(EPSILON*EPSILON) );
			real inv_r1 = (real)1.0 / real_sqrt( r2 + (real)(EPSILON*EPSILON) );
			#else
			real inv_r2 = (real)0.0 < r2 ? (real)1.0 /          (r2) : (real)0.0;
			real inv_r1 = (real)0.0 < r2 ? (real)1.0 / real_sqrt(r2) : (real)0.0;
			#endif

			// gravitational potential of body j
			// (G and m_i are taken into account later)
			real pot_j = m_j * inv_r1;

			real F[3] = {
				pot_j * inv_r2 * r_ij[0],
				pot_j * inv_r2 * r_ij[1],
				pot_j * inv_r2 * r_ij[2]
			};

			for(unsigned int k = 0u; k < 3u; k++) {
				#ifdef ENABLE_KAHAN_SUMMATION_SOLVER
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				const real y = F[k] - acc_com[k];
				volatile real t = acc_sum[k] + y;
				volatile real z = t - acc_sum[k];
				acc_com[k] = z - y;
				acc_sum[k] = t;
				#else
				// Naïve summation
				acc_sum[k] += F[k];
				#endif
			}

			{
				#ifdef ENABLE_KAHAN_SUMMATION_ENERGY
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				const real y = pot_j - pot_com;
				volatile real t = pot_sum + y;
				volatile real z = t - pot_sum;
				pot_com = z - y;
				pot_sum = t;
				#else
				// Naïve summation
				pot_sum += pot_j;
				#endif
			}
		}
		*/

		acc[3u*idx+0u] = (real)(-G) * m_i * acc_sum[0u];
		acc[3u*idx+1u] = (real)(-G) * m_i * acc_sum[1u];
		acc[3u*idx+2u] = (real)(-G) * m_i * acc_sum[2u];
		pot[   idx   ] = (real)(-G) * pot_sum;
	}

	return STARFLOOD_SUCCESS;
}
#endif

// Solve gravity using brute-force method (particle-particle method)
int solve_gravity_part_part(const sim_t* restrict sim_ptr) {
	TIMING_INIT();

	sim_t sim = *sim_ptr;

	unsigned int N = sim.N;
	real* pos = sim_find(&sim, SIM_POS);
	real* vel = sim_find(&sim, SIM_VEL);
	real* acc = sim_find(&sim, SIM_ACC);
	real* mas = sim_find(&sim, SIM_MAS);
	real* pot = sim_find(&sim, SIM_POT);

	#ifdef PAIRWISE_SOLVER_DECIMATION
	const unsigned int j_length = (          N / (unsigned int)PAIRWISE_SOLVER_DECIMATION);
	const unsigned int j_offset = (step_number % (unsigned int)PAIRWISE_SOLVER_DECIMATION) * j_length;
	#endif

	TIMING_START();

	#ifdef _OPENMP
		#ifdef ENABLE_OFFLOAD_SIM
		#pragma omp target teams distribute parallel for map(tofrom: acc[:3u*N], pot[:N]) map(to: pos[:3u*N], mas[:N])
		#else
		#pragma omp parallel for schedule(dynamic, 256)
		#endif
	#endif
	for(unsigned int i = 0u; i < N; i++) {
		real acc_sum[3] = {(real)0.0, (real)0.0, (real)0.0};

		#ifdef ENABLE_KAHAN_SUMMATION_SOLVER
		real acc_com[3] = {(real)0.0, (real)0.0, (real)0.0};
		#endif

		real pot_sum = (real)0.0;

		#ifdef ENABLE_KAHAN_SUMMATION_ENERGY
		real pot_com = (real)0.0;
		#endif

		#ifdef ENABLE_SPH
		real rho_sum = (real)0.0;
		#ifdef ENABLE_KAHAN_SUMMATION_SOLVER
		real rho_com = (real)0.0;
		#endif
		#endif

		const real m_i = mas[i]; // body i's mass

		// body i's position
		const real r_i[3] = {
			pos[3u*i+0u],
			pos[3u*i+1u],
			pos[3u*i+2u]
		};

		#ifdef PAIRWISE_SOLVER_DECIMATION
		for(unsigned int j = j_offset; j < (j_offset+j_length); j++)
		#else
		for(unsigned int j = 0u; j < N; j++)
		#endif
		{
			/*
			if(i == j) {
				continue;
			}
			*/

			const real m_j = mas[j]; // body j's mass

			// body j's position
			const real r_j[3] = {
				pos[3u*j+0u],
				pos[3u*j+1u],
				pos[3u*j+2u]
			};

			// vector to body i's position from body j's position
			const real r_ij[3] = {
				r_i[0] - r_j[0],
				r_i[1] - r_j[1],
				r_i[2] - r_j[2]
			};

			// squared distance between bodies i and j
			const real r2 = (r_ij[0]*r_ij[0])+(r_ij[1]*r_ij[1])+(r_ij[2]*r_ij[2]);

			#ifdef EPSILON
			const real inv_r2 = (real)1.0 /          ( r2 + (real)(EPSILON*EPSILON) );
			const real inv_r1 = (real)1.0 / real_sqrt( r2 + (real)(EPSILON*EPSILON) );
			#else
			const real inv_r2 = (real)0.0 < r2 ? (real)1.0 /          (r2) : (real)0.0;
			const real inv_r1 = (real)0.0 < r2 ? (real)1.0 / real_sqrt(r2) : (real)0.0;
			#endif

			// gravitational potential of body j
			// (G and m_i are taken into account later)
			const real pot_j = m_j * inv_r1;

			const real F[3] = {
				pot_j * inv_r2 * r_ij[0],
				pot_j * inv_r2 * r_ij[1],
				pot_j * inv_r2 * r_ij[2]
			};

			for(unsigned int k = 0u; k < 3u; k++) {
				#ifdef ENABLE_KAHAN_SUMMATION_SOLVER
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				const real y = F[k] - acc_com[k];
				volatile real t = acc_sum[k] + y;
				volatile real z = t - acc_sum[k];
				acc_com[k] = z - y;
				acc_sum[k] = t;
				#else
				// Naïve summation
				acc_sum[k] += F[k];
				#endif
			}

			{
				#ifdef ENABLE_KAHAN_SUMMATION_ENERGY
				// Kahan summation
				// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				const real y = pot_j - pot_com;
				volatile real t = pot_sum + y;
				volatile real z = t - pot_sum;
				pot_com = z - y;
				pot_sum = t;
				#else
				// Naïve summation
				pot_sum += pot_j;
				#endif
			}
		}

		#ifdef PAIRWISE_SOLVER_DECIMATION
		acc[3u*i+0u] = (real)PAIRWISE_SOLVER_DECIMATION * (real)(-G) * m_i * acc_sum[0u];
		acc[3u*i+1u] = (real)PAIRWISE_SOLVER_DECIMATION * (real)(-G) * m_i * acc_sum[1u];
		acc[3u*i+2u] = (real)PAIRWISE_SOLVER_DECIMATION * (real)(-G) * m_i * acc_sum[2u];
		pot[   i   ] = (real)PAIRWISE_SOLVER_DECIMATION * (real)(-G) * pot_sum;
		#else
		acc[3u*i+0u] = (real)(-G) * m_i * acc_sum[0u];
		acc[3u*i+1u] = (real)(-G) * m_i * acc_sum[1u];
		acc[3u*i+2u] = (real)(-G) * m_i * acc_sum[2u];
		pot[   i   ] = (real)(-G) * pot_sum;
		#endif
	}

	TIMING_STOP();
	TIMING_PRINT("solve_gravity()", "brute_force");

	return STARFLOOD_SUCCESS;
}
