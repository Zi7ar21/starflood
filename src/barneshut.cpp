#include "barneshut.h"
#include "constants.h"

#include <stack>
#include <vector>

class Node {
	public:
		int parent, child[4], id, bodies;
		real m, mx, my;

		Node() {
			parent = -1;
			child[0] = -1;
			child[1] = -1;
			child[2] = -1;
			child[3] = -1;
			id = -1;
			bodies = 0;
			m = (real)0;
			mx = (real)0;
			my = (real)0;
		}

		Node(int _parent, int _id, real _m, real _mx, real _my) {
			parent = _parent;
			child[0] = -1;
			child[1] = -1;
			child[2] = -1;
			child[3] = -1;
			id = _id;
			bodies = 1;
			m = _m;
			mx = _mx;
			my = _my;
		}
};

void barnes_hut(real* a, real* x, real* m, int N) {
	real bound_min_x = x[0], bound_max_x = x[0], bound_min_y = x[1], bound_max_y = x[1];

	for(int i = 1; i < N; i++) {
		real xi[2] = {x[2 * i + 0], x[2 * i + 1]};

		bound_min_x = xi[0] < bound_min_x ? xi[0] : bound_min_x;
		bound_max_x = xi[0] > bound_max_x ? xi[0] : bound_max_x;
		bound_min_y = xi[1] < bound_min_y ? xi[1] : bound_min_y;
		bound_max_y = xi[1] > bound_max_y ? xi[1] : bound_max_y;
	}

	bound_min_x -= (real)0.3;
	bound_max_x += (real)0.3;
	bound_min_y -= (real)0.3;
	bound_max_y += (real)0.3;

	std::vector<Node> tree;

	// build tree
	{
	tree.push_back(Node(-1, 0, m[0], x[2 * 0 + 0] * m[0], x[2 * 0 + 1] * m[0]));

	// insert particles
	for(int i = 1; i < N; i++) {
		float xi[2] = {x[2 * i + 0], x[2 * i + 1]}; // position of the body we are inserting

		real cur_bound_min_x = bound_min_x;
		real cur_bound_max_x = bound_max_x;
		real cur_bound_min_y = bound_min_y;
		real cur_bound_max_y = bound_max_y;

		int cur_node = 0;

		while(true) {
			real cur_bound_center_x = (real)0.5 * (cur_bound_min_x + cur_bound_max_x);
			real cur_bound_center_y = (real)0.5 * (cur_bound_min_y + cur_bound_max_y);

			if(tree[cur_node].child[0] >= 0
			|| tree[cur_node].child[1] >= 0
			|| tree[cur_node].child[2] >= 0
			|| tree[cur_node].child[3] >= 0) {
				// node is subdivided
				int quadrant = 0;

				if(xi[1] > cur_bound_center_y) {
					if(xi[0] > cur_bound_center_x) {
						quadrant = 1;
					} else {
						quadrant = 0;
					}
				} else {
					if(xi[0] > cur_bound_center_x) {
						quadrant = 3;
					} else {
						quadrant = 2;
					}
				}

				if(tree[cur_node].child[quadrant] < 0) {
					// node isn't allocated
					tree[cur_node].child[quadrant] = tree.size();

					tree.push_back(Node(cur_node, i, m[i], m[i] * xi[0], m[i] * xi[1]));

					while(cur_node >= 0) {
						tree[cur_node].bodies += 1;
						tree[cur_node].m += m[i];
						tree[cur_node].mx += m[i] * xi[0];
						tree[cur_node].my += m[i] * xi[1];

						cur_node = tree[cur_node].parent;
					}

					break;
				} else {
					// node is allocated
					cur_node = tree[cur_node].child[quadrant];

					if(quadrant == 0) {
						cur_bound_max_x = cur_bound_center_x;
						cur_bound_min_y = cur_bound_center_y;
					}
					if(quadrant == 1) {
						cur_bound_min_x = cur_bound_center_x;
						cur_bound_min_y = cur_bound_center_y;
					}
					if(quadrant == 2) {
						cur_bound_max_x = cur_bound_center_x;
						cur_bound_max_y = cur_bound_center_y;
					}
					if(quadrant == 3) {
						cur_bound_min_x = cur_bound_center_x;
						cur_bound_max_y = cur_bound_center_y;
					}
				}
			} else {
				// node isn't subdivided
				if(tree[cur_node].id < 0) {
					// no particle contained
					// actually wait wtf this shouldn't happen
					break;
				} else {
					// node contains a particle
					int tmp_id = tree[cur_node].id;

					tree[cur_node].id = -1;

					int quadrant = 0;

					if(x[2 * tmp_id + 1] > cur_bound_center_y) {
						if(x[2 * tmp_id + 0] > cur_bound_center_x) {
							quadrant = 1;
						} else {
							quadrant = 0;
						}
					} else {
						if(x[2 * tmp_id + 0] > cur_bound_center_x) {
							quadrant = 3;
						} else {
							quadrant = 2;
						}
					}

					tree.push_back(Node(cur_node, tmp_id, m[tmp_id], m[tmp_id] * x[2 * tmp_id + 0], m[tmp_id] * x[2 * tmp_id + 1]));
				}
			}
		}
	}

	// finish center of mass calculation
	for(int i = 0; i < tree.size(); i++) {
			tree[i].mx /= tree[i].m != (real)0 ? tree[i].m : (real)1;
			tree[i].my /= tree[i].m != (real)0 ? tree[i].m : (real)1;
	}
	}

	// calculate forces
	{
		const real theta = (real)0.5;

		#pragma omp simd
		for(int i = 0; i < 2 * N; i++) a[i] = (real)0;

		for(int i = 0; i < N; i++) {
			std::stack<int> to_visit;

			to_visit.push(0);

			float rx = x[2 * i + 0];
			float ry = x[2 * i + 1];

			//if(bound_min_x > rx || rx > bound_max_x || bound_min_y > ry || ry > bound_max_y) continue;

			//width / distance falls below a chosen threshold (a parameter named theta),

			float domain_width = (bound_max_x - bound_min_x);

			while(!to_visit.empty()) {
				int depth = 0;

				int cur_node = to_visit.top();

				to_visit.pop();

				real dx = x[2 * tree[cur_node].id + 0] - x[2 * i + 0], dy = x[2 * tree[cur_node].id + 1] - x[2 * i + 1];

				real r2 = dx*dx+dy*dy;

				real width = exp2f(-(float)depth) * domain_width;

				if(tree[cur_node].id >= 0 || tree[cur_node].bodies == 1 || (width/sqrtf(r2)) < theta || (tree[cur_node].child[0] >= 0 && tree[cur_node].child[1] >= 0 && tree[cur_node].child[2] >= 0 && tree[cur_node].child[3] >= 0)) {
					if(tree[cur_node].id == i) continue; // don't calculate forces with itself	// calculate forces
					//std::cout << "guh" << std::endl;
					a[2 * i + 0] += (dx / sqrtf(r2)) * tree[cur_node].m * m[i] * (1.0 / r2);
					a[2 * i + 1] += (dy / sqrtf(r2)) * tree[cur_node].m * m[i] * (1.0 / r2);
				} else {
					// visit child nodes
					// if node exists, add it
					if(tree[cur_node].child[0] >= 0) to_visit.push(tree[cur_node].child[0]);
					if(tree[cur_node].child[1] >= 0) to_visit.push(tree[cur_node].child[1]);
					if(tree[cur_node].child[2] >= 0) to_visit.push(tree[cur_node].child[2]);
					if(tree[cur_node].child[3] >= 0) to_visit.push(tree[cur_node].child[3]);
				}
			}
		}

		#pragma omp simd
		for(int i = 0; i < 2 * N; i++) a[i] *= G * m[i / 2];
	}
}
