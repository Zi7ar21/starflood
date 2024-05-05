#include <barnes-hut.hpp>

#include <cmath>
#include <cstdio>
#include <rng.hpp>
//#include <omp.h>

#ifdef STARFLOOD_RENDER_INTERACTIONS
#include <graphics.hpp>
#endif

void BarnesHut(std::vector<Node> &tree, float* image, int w, int h, int* ids, real* mas, real* pos, real* acc, real* pen, int N, int step_num) {
	for(int i = 0; i < N; i++) ids[i] = -1;

	tree.clear();

	#ifndef TREE_FIT
	tree.push_back(Node(-1,(real)-1.0,(real)-1.0,(real)1.0,(real)1.0));
	#else
	tree.push_back(Node(-1,pos[0],pos[1],pos[2],pos[0],pos[1],pos[2]));
	#endif

	{
		const real extrspace = 0.1;

		#ifdef TREE_FIT
		for(int i = 1; i < N; i++) {
			tree[0].x_min = fmin(tree[0].x_min,pos[3*i+0]);
			tree[0].y_min = fmin(tree[0].y_min,pos[3*i+1]);
			tree[0].z_min = fmin(tree[0].z_min,pos[3*i+2]);
			tree[0].x_max = fmax(tree[0].x_max,pos[3*i+0]);
			tree[0].y_max = fmax(tree[0].y_max,pos[3*i+1]);
			tree[0].z_max = fmax(tree[0].z_max,pos[3*i+2]);
		}


		tree[0].x_min -= extrspace;
		tree[0].y_min -= extrspace;
		tree[0].z_min -= extrspace;
		tree[0].x_max += extrspace;
		tree[0].y_max += extrspace;
		tree[0].z_max += extrspace;
		#endif

		#ifdef JITTER_TREE
		uint32_t ns = (uint32_t)step_num+(uint32_t)37; // set the random number generator seed

		float z0 = urand(&ns);
		float z1 = urand(&ns);
		float z2 = urand(&ns);

		tree[0].x_min += extrspace*(real)(z0-0.5f);
		tree[0].y_min += extrspace*(real)(z1-0.5f);
		tree[0].z_min += extrspace*(real)(z2-0.5f);
		tree[0].x_max += extrspace*(real)(z0-0.5f);
		tree[0].y_max += extrspace*(real)(z1-0.5f);
		tree[0].z_max += extrspace*(real)(z2-0.5f);
		#endif
	}

	/*
	Quadrants:
	+---+---+
	| 2 | 3 |
	+---+---+
	| 0 | 1 |
	+---+---+
	*/

	// Build the Barnes-Hut Tree
	for(int i = 0; i < N; i++) {
		int current_node = 0;

		int depth = 0;

		while(depth <= 12) {
			real x_min = tree[current_node].x_min;
			real y_min = tree[current_node].y_min;
			real z_min = tree[current_node].z_min;
			real x_max = tree[current_node].x_max;
			real y_max = tree[current_node].y_max;
			real z_max = tree[current_node].z_max;

			// don't continue if particle is outside of bounds
			if((x_max < pos[3*i+0]) || (pos[3*i+0] < x_min)
			|| (y_max < pos[3*i+1]) || (pos[3*i+1] < y_min)
			|| (z_max < pos[3*i+2]) || (pos[3*i+2] < z_min)) break;

			real hx = ((real)0.5*tree[current_node].x_min)+((real)0.5*tree[current_node].x_max);
			real hy = ((real)0.5*tree[current_node].y_min)+((real)0.5*tree[current_node].y_max);
			real hz = ((real)0.5*tree[current_node].z_min)+((real)0.5*tree[current_node].z_max);

			if(tree[current_node].child0 != -1
			|| tree[current_node].child1 != -1
			|| tree[current_node].child2 != -1
			|| tree[current_node].child3 != -1
			|| tree[current_node].child4 != -1
			|| tree[current_node].child5 != -1
			|| tree[current_node].child6 != -1
			|| tree[current_node].child7 != -1) {
				// node is split
				int q = (pos[3*i+0] > hx ? 1 : 0)+(pos[3*i+1] > hy ? 2 : 0)+(pos[3*i+2] > hz ? 4 : 0);

				switch(q) {
					case 0:
					current_node = tree[current_node].child0;
					break;
					case 1:
					current_node = tree[current_node].child1;
					break;
					case 2:
					current_node = tree[current_node].child2;
					break;
					case 3:
					current_node = tree[current_node].child3;
					break;
					case 4:
					current_node = tree[current_node].child4;
					break;
					case 5:
					current_node = tree[current_node].child5;
					break;
					case 6:
					current_node = tree[current_node].child6;
					break;
					case 7:
					current_node = tree[current_node].child7;
					break;
					default:
					printf("Severe error!\n");
					break;
				}

				depth++;

				continue;
			}

			// node is not split
			if(tree[current_node].id != -1) {
				/*
				Quadrants:
				+---+---+
				| 2 | 3 |
				+---+---+
				| 0 | 1 |
				+---+---+
				*/

				// node has particle
				tree.push_back(Node(current_node,x_min,y_min,z_min,   hx,   hy,   hz));
				tree[current_node].child0 = tree.size() - 1;

				tree.push_back(Node(current_node,   hx,y_min,z_min,x_max,   hy,   hz));
				tree[current_node].child1 = tree.size() - 1;

				tree.push_back(Node(current_node,x_min,   hy,z_min,   hx,y_max,   hz));
				tree[current_node].child2 = tree.size() - 1;

				tree.push_back(Node(current_node,   hx,   hy,z_min,x_max,y_max,   hz));
				tree[current_node].child3 = tree.size() - 1;

				tree.push_back(Node(current_node,x_min,y_min,   hz,   hx,   hy,z_max));
				tree[current_node].child4 = tree.size() - 1;

				tree.push_back(Node(current_node,   hx,y_min,   hz,x_max,   hy,z_max));
				tree[current_node].child5 = tree.size() - 1;

				tree.push_back(Node(current_node,x_min,   hy,   hz,   hx,y_max,z_max));
				tree[current_node].child6 = tree.size() - 1;

				tree.push_back(Node(current_node,   hx,   hy,   hz,x_max,y_max,z_max));
				tree[current_node].child7 = tree.size() - 1;

				int id = tree[current_node].id;

				tree[current_node].id = -1;

				int q = (pos[3*id+0] > hx ? 1 : 0)+(pos[3*id+1] > hy ? 2 : 0)+(pos[3*id+2] > hz ? 4 : 0);

				switch(q) {
					case 0:
					tree[tree[current_node].child0].id = id;
					tree[tree[current_node].child0].np = 1;
					break;
					case 1:
					tree[tree[current_node].child1].id = id;
					tree[tree[current_node].child1].np = 1;
					break;
					case 2:
					tree[tree[current_node].child2].id = id;
					tree[tree[current_node].child2].np = 1;
					break;
					case 3:
					tree[tree[current_node].child3].id = id;
					tree[tree[current_node].child3].np = 1;
					break;
					case 4:
					tree[tree[current_node].child4].id = id;
					tree[tree[current_node].child4].np = 1;
					break;
					case 5:
					tree[tree[current_node].child5].id = id;
					tree[tree[current_node].child5].np = 1;
					break;
					case 6:
					tree[tree[current_node].child6].id = id;
					tree[tree[current_node].child6].np = 1;
					break;
					case 7:
					tree[tree[current_node].child7].id = id;
					tree[tree[current_node].child7].np = 1;
					break;
					default:
					printf("Severe error!\n");
					break;
				}

				continue;
			}

			// node doesn't have particle
			tree[current_node].np = 1;
			tree[current_node].id = i;

			int tmpid = current_node;
			while(tmpid != -1) {
				tree[tmpid].m += mas[i];
				tree[tmpid].x += mas[i]*pos[3*i+0];
				tree[tmpid].y += mas[i]*pos[3*i+1];
				tree[tmpid].z += mas[i]*pos[3*i+2];
				tmpid = tree[tmpid].parent;
			}

			break;
		}

		if(depth == 8) {
			//failures++;
			//printf("A %d\n", failures);
			//return EXIT_FAILURE;
		}

		//break;
	}

	// Calculate Center of Mass
	for(size_t i = 0; i < tree.size(); i++) {
		tree[i].x = (tree[i].m != (real)0) ? (tree[i].x / tree[i].m) : tree[i].x;
		tree[i].y = (tree[i].m != (real)0) ? (tree[i].y / tree[i].m) : tree[i].y;
		tree[i].z = (tree[i].m != (real)0) ? (tree[i].z / tree[i].m) : tree[i].z;
	}

	{
		std::vector<Node> ntree;

		std::stack<int> to_visit;

		to_visit.push(0);

		while(!to_visit.empty()) {
				int cur_node = to_visit.top();

				to_visit.pop();

				ntree.push_back(tree[cur_node]);
		}

		tree = ntree;
	}

	for(int i = 0; i < N; i++) pen[i] = 0;

	// Compute Forces
	#pragma omp parallel for schedule(dynamic,1024)
	for(int i = 0; i < N; i++) {
		// don't continue if particle is outside of bounds
		if((tree[0].x_max < pos[3*i+0]) || (pos[3*i+0] < tree[0].x_min)
		|| (tree[0].y_max < pos[3*i+1]) || (pos[3*i+1] < tree[0].y_min)
		|| (tree[0].z_max < pos[3*i+2]) || (pos[3*i+2] < tree[0].z_min)) continue;

		acc[3*i+0] = (real)0;
		acc[3*i+1] = (real)0;
		acc[3*i+2] = (real)0;

		real u_i = pos[3*i+0];
		real v_i = pos[3*i+1];
		real w_i = pos[3*i+2];

		std::stack<int> to_visit;

		to_visit.push(0);

		int depth = 0;

		while(!to_visit.empty()) {
			//real hx = ((real)0.5*x_min)+((real)0.5*x_max);
			//real hy = ((real)0.5*x_min)+((real)0.5*x_max);

			int cur_node = to_visit.top();

			to_visit.pop();

			// if this is empty or contains itself, we need not go any further
			if((tree[cur_node].np < 1) || (tree[cur_node].id == i)) continue;

			real x_min = tree[cur_node].x_min, y_min = tree[cur_node].y_min, z_min = tree[cur_node].z_min, x_max = tree[cur_node].x_max, y_max = tree[cur_node].y_max, z_max = tree[cur_node].z_max;

			#ifdef STARFLOOD_RENDER_INTERACTIONS
			if(i == 0) {
				drawLineUV(image, w, h, 1., 1., 0., 1., x_min, y_min, x_max, y_min);
				drawLineUV(image, w, h, 1., 1., 0., 1., x_min, y_min, x_min, y_max);
				drawLineUV(image, w, h, 1., 1., 0., 1., x_min, y_max, x_max, y_max);
				drawLineUV(image, w, h, 1., 1., 0., 1., x_max, y_min, x_max, y_max);
			}
			#endif

			real width = ((real)0.5*(x_max - x_min))+((real)0.5*(y_max - y_min))+((real)0.5*(z_max - z_min));

			real dx = tree[cur_node].x;
			real dy = tree[cur_node].y;
			real dz = tree[cur_node].z;
			real dm = tree[cur_node].m;

			// probably empty, idk just to be safe
			if(dm < (real)0.001) continue;

			if(x_min < u_i && u_i < x_max
			&& y_min < v_i && v_i < y_max
			&& z_min < w_i && w_i < z_max) {
				// this particle is contained in the node we are evaluating so let's
				// calculate what it would be like with it removed
				dx *= dm;
				dy *= dm;
				dz *= dm;
				dm -= mas[i];
				dx -= mas[i]*pos[3*i+0];
				dy -= mas[i]*pos[3*i+1];
				dz -= mas[i]*pos[3*i+2];
				dx /= dm;
				dy /= dm;
				dz /= dm;
			}

			dx = dx-u_i;
			dy = dy-v_i;
			dz = dz-w_i;

			real dist = sqrtf((dx*dx)+(dy*dy)+(dz*dz));

			// again, probably empty, just to be safe
			if(dm < (real)0.001) continue;

			const real theta = THETA;

			if((width/dist) > theta) {
				if(tree[cur_node].child0 != -1) to_visit.push(tree[cur_node].child0);
				if(tree[cur_node].child1 != -1) to_visit.push(tree[cur_node].child1);
				if(tree[cur_node].child2 != -1) to_visit.push(tree[cur_node].child2);
				if(tree[cur_node].child3 != -1) to_visit.push(tree[cur_node].child3);
				if(tree[cur_node].child4 != -1) to_visit.push(tree[cur_node].child4);
				if(tree[cur_node].child5 != -1) to_visit.push(tree[cur_node].child5);
				if(tree[cur_node].child6 != -1) to_visit.push(tree[cur_node].child6);
				if(tree[cur_node].child7 != -1) to_visit.push(tree[cur_node].child7);
				continue;
			}

			#ifdef STARFLOOD_RENDER_INTERACTIONS
			if(i == 0) {
				//drawLineUV(image, w, h, 1., 1., 0., 1., x_min, y_min, x_max, y_min);
				//drawLineUV(image, w, h, 1., 1., 0., 1., x_min, y_min, x_min, y_max);
				//drawLineUV(image, w, h, 1., 1., 0., 1., x_min, y_max, x_max, y_max);
				//drawLineUV(image, w, h, 1., 1., 0., 1., x_max, y_min, x_max, y_max);
				drawLineUV(image, w, h, 1., 0., 0., 1., u_i, v_i, tree[cur_node].x, tree[cur_node].y);
			}
			#endif

			// F = G*m_i*m_j/r^2
			acc[3*i+0] += ((dx/dist)*mas[i]*dm)/((dx*dx)+(dy*dy)+(dz*dz)+(real)0.0001);
			acc[3*i+1] += ((dy/dist)*mas[i]*dm)/((dx*dx)+(dy*dy)+(dz*dz)+(real)0.0001);
			acc[3*i+2] += ((dz/dist)*mas[i]*dm)/((dx*dx)+(dy*dy)+(dz*dz)+(real)0.0001);
			pen[  i  ] += (mas[i]*dm)/sqrt((dx*dx)+(dy*dy)+(dz*dz)+0.0001); // U = G*m_i*m_j/r_ij

			//int quad = (p.x > hx ? 1 : 0) + (p.y > hy ? 2 : 0);

			/*
			switch(quad) {
				case 0:
				bound_max = h;
				break;

				case 1:
				bound_min.x = h.x;
				bound_max.y = h.y;
				break;

				case 2:
				bound_max.x = h.x;
				bound_min.y = h.y;
				break;

				case 3:
				bound_min = h;
				break;

				default:
				break;
			}
			*/
		}
	}
}