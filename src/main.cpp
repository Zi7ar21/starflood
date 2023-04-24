#include <iostream>
#include <cstdint>
//#include <cstdlib>
#include <cmath>
#include <stack>
#include <vector>

#include "common.hpp"
#include "config.h"
#include "constants.h"
#include "image.hpp"

// Real number type (float or double)
typedef float Real;

// Barnes hut tree node
class Node {
	public:
		int id; // node identifier
		int parent; // <0: this is the root node
		int child1; // <0: not allocated
		int child2; // <0: not allocated
		int child3; // <0: not allocated
		int child4; // <0: not allocated
		int depth; // 0: root

		int particle_id; // <0: no particles, >=0: particle id

		int particles_contained; // number of particles contained in this and child nodes
		float mass; // mass of particles in this and child nodes
		float mass_x; // center of mass x
		float mass_y; // center of mass y

		Node(int _id): id(_id), parent(-1), child1(-1), child2(-1), child3(-1), child4(-1), depth(-1), particle_id(-1), particles_contained(0), mass(0.0), mass_x(0.0), mass_y(0.0) {}
};

class BarnesHut {
	public:
		unsigned int N; // number of particles in the simulation
		Real* m; // particle mass buffer
		Real* r; // particle position buffer
		Real* v; // particle velocity buffer
		Real* a; // particle acceleration buffer

		// Tree bounds
		Real bound_min_x;
		Real bound_min_y;
		Real bound_max_x;
		Real bound_max_y;

		std::vector<Node> node;

		BarnesHut(double _bound_min_x, double _bound_min_y, double _bound_max_x, double _bound_max_y): bound_min_x(_bound_min_x), bound_min_y(_bound_min_y), bound_max_x(_bound_max_x), bound_max_y(_bound_max_y) {
			//node.push_back(Node(0));
		}

		void Build() {
			node.clear(); // reset tree

			node.push_back(Node(0)); // create root node

			// initialize the root node
			node[0].depth = 0;
			node[0].particle_id = 0;
			node[0].particles_contained = 1;

			node[0].mass = m[0];
			node[0].mass_x = m[0] * r[2 * 0 + 0];
			node[0].mass_y = m[0] * r[2 * 0 + 1];

			//std::cout << "i = 0" << std::endl;
			//std::cout << "rx = " << r[2*0+0] << std::endl;
			//std::cout << "ry = " << r[2*0+1] << std::endl;

			// insert the rest of the particles
			for(int i = 1; i < N; i++) {
				float rx = r[2 * i + 0];
				float ry = r[2 * i + 1];
				//std::cout << "i = " << i << std::endl;
				//std::cout << "rx = " << rx << std::endl;
				//std::cout << "ry = " << ry << std::endl;

				if((bound_min_x > rx) || (rx > bound_max_x) || (bound_min_y > ry) || (ry > bound_max_y)) continue; // skip if out of bounds

				int cur_node = 0; // current node identifier

				// current node bounds
				Real current_bound_min_x = bound_min_x;
				Real current_bound_min_y = bound_min_y;
				Real current_bound_max_x = bound_max_x;
				Real current_bound_max_y = bound_max_y;

				while(true) {
					//std::cout << "INSIDE NODE " << cur_node << std::endl;

					// Find the center of the current node
					double current_bound_center_x = 0.5 * (current_bound_max_x + current_bound_min_x);
					double current_bound_center_y = 0.5 * (current_bound_max_y + current_bound_min_y);

					//std::cout << "current_bound_min_x: " << current_bound_min_x << std::endl;
					//std::cout << "current_bound_min_y: " << current_bound_min_y << std::endl;
					//std::cout << "current_bound_max_x: " << current_bound_max_x << std::endl;
					//std::cout << "current_bound_max_y: " << current_bound_max_y << std::endl;
					//std::cout << "current_bound_center_x: " << current_bound_center_x << std::endl;
					//std::cout << "current_bound_center_y: " << current_bound_center_y << std::endl;

					//std::cout << "children: [" << node[cur_node].child1 << ", "<< node[cur_node].child2 << ", " << node[cur_node].child3 << ", " << node[cur_node].child4 << "]" << std::endl;
					if((node[cur_node].child1 < 0) && (node[cur_node].child2 < 0) && (node[cur_node].child3 < 0) && (node[cur_node].child4 < 0)) {
						// node is not subdivided
						//std::cout << "NODE " << cur_node << " IS NOT SUBDIVIDED" << std::endl;
						if(node[cur_node].particle_id < 0) {
							// node doesn't contain a particle
							//std::cout << "NODE " << cur_node << " DOESN'T CONTAIN A PARTICLE" << std::endl;

							node[cur_node].particle_id = i;

							// update mass in this and parent nodes
							while(true) {
								if(node[cur_node].parent < 0) {
									// this node is the root node
									break;
								} else {
									// this node is not the root node
									node[cur_node].particles_contained += 1;
									node[cur_node].mass += m[i];
									node[cur_node].mass_x += m[i] * r[2 * i + 0];
									node[cur_node].mass_y += m[i] * r[2 * i + 1];
									cur_node = node[cur_node].parent;
								}
							}

							break;
						} else {
							// node contains a particle
							//std::cout << "NODE " << cur_node << " CONTAINS A PARTICLE" << std::endl;

							int new_node_id = node.size();
							//std::cout << "Creating node " << new_node_id << std::endl;
							node.push_back(Node(new_node_id));
							node[new_node_id].depth = node[cur_node].depth + 1;

							int new_node_id_1 = node.size();
							//std::cout << "Creating node " << new_node_id_1 << std::endl;
							node.push_back(Node(new_node_id_1));
							node[new_node_id_1].depth = node[cur_node].depth + 1;

							int new_node_id_2 = node.size();
							//std::cout << "Creating node " << new_node_id_2 << std::endl;
							node.push_back(Node(new_node_id_2));
							node[new_node_id_2].depth = node[cur_node].depth + 1;

							int new_node_id_3 = node.size();
							//std::cout << "Creating node " << new_node_id_3 << std::endl;
							node.push_back(Node(new_node_id_3));
							node[new_node_id_3].depth = node[cur_node].depth + 1;

							int cur_node_particle = node[cur_node].particle_id;

							double cx = r[2 * cur_node_particle + 0];
							double cy = r[2 * cur_node_particle + 1];

							node[cur_node].particle_id = -1;
							node[new_node_id].particle_id = cur_node_particle;
							node[new_node_id].particles_contained = 1;
							node[new_node_id].mass = m[cur_node_particle];
							node[new_node_id].mass_x = m[cur_node_particle] * r[2 * cur_node_particle + 0];
							node[new_node_id].mass_y = m[cur_node_particle] * r[2 * cur_node_particle + 1];
							node[new_node_id].parent = cur_node;

							int quadrant = cx >= current_bound_center_x ? (cy >= current_bound_center_y ? 1 : 4) : (cy >= current_bound_center_y ? 2 : 3);

							if(quadrant == 1) {
								node[cur_node].child1 = new_node_id;
								node[cur_node].child2 = new_node_id_1;
								node[cur_node].child3 = new_node_id_2;
								node[cur_node].child4 = new_node_id_3;
							}
							if(quadrant == 2) {
								node[cur_node].child1 = new_node_id_1;
								node[cur_node].child2 = new_node_id;
								node[cur_node].child3 = new_node_id_2;
								node[cur_node].child4 = new_node_id_3;
							}
							if(quadrant == 3) {
								node[cur_node].child1 = new_node_id_1;
								node[cur_node].child2 = new_node_id_2;
								node[cur_node].child3 = new_node_id;
								node[cur_node].child4 = new_node_id_3;
							}
							if(quadrant == 4) {
								node[cur_node].child1 = new_node_id_1;
								node[cur_node].child2 = new_node_id_2;
								node[cur_node].child3 = new_node_id_3;
								node[cur_node].child4 = new_node_id;
							}
						}
					} else {
						// node is subdivided
						//std::cout << "NODE " << cur_node << " IS SUBDIVIDED" << std::endl;

						// move into appropriate child
						int quadrant = rx >= current_bound_center_x ? (ry >= current_bound_center_y ? 1 : 4) : (ry >= current_bound_center_y ? 2 : 3);
						if(quadrant == 1) {
							cur_node = node[cur_node].child1;
							current_bound_min_x = current_bound_center_x;
							current_bound_min_y = current_bound_center_y;
						}
						if(quadrant == 2) {
							cur_node = node[cur_node].child2;
							current_bound_max_x = current_bound_center_x;
							current_bound_min_y = current_bound_center_y;
						}
						if(quadrant == 3) {
							cur_node = node[cur_node].child3;
							current_bound_max_x = current_bound_center_x;
							current_bound_max_y = current_bound_center_y;
						}
						if(quadrant == 4) {
							cur_node = node[cur_node].child4;
							current_bound_min_x = current_bound_center_x;
							current_bound_max_y = current_bound_center_y;
						}
					}
				}
			}

			for(int i = 0; i < node.size(); i++) {
				// divide the sum of position times mass by the number of masses (final step in calculating center of mass)
				node[i].mass_x /= (node[i].particles_contained > 1) ? (double)node[i].particles_contained : 1.0;
				node[i].mass_y /= (node[i].particles_contained > 1) ? (double)node[i].particles_contained : 1.0;
			}

			//std::cout << "Done builidng tree" << std::endl;
		}

		void UpdateAcceleration() {
			const float theta = 1.0;

			for(int i = 0; i < (2 * N); i++) a[i] = (float)0;

			for(int i = 0; i < N; i++) {
				//std::cout << "a" << i << std::endl;
				std::stack<int> to_visit;

				to_visit.push(0);

				float rx = r[2 * i + 0];
				float ry = r[2 * i + 1];

				//if(bound_min_x > rx || rx > bound_max_x || bound_min_y > ry || ry > bound_max_y) continue;

				//width / distance falls below a chosen threshold (a parameter named theta),

				float domain_width = (bound_max_x - bound_min_x);

				while(!to_visit.empty()) {
					int cur_node = to_visit.top();

					//std::cout << "visiting " << cur_node << std::endl;
					to_visit.pop();

					Real width = exp2f((float)-node[cur_node].depth) * domain_width;

					Real dx = (node[cur_node].mass_x - rx);
					Real dy = (node[cur_node].mass_y - ry);
					Real d2 = (dx*dx)+(dy*dy);
					if(d2 < (Real)(0.01 * 0.01)) continue; // too close
					Real d1 = sqrtf(d2);

					if(((width / d1) < theta) || (node[cur_node].child1 == -1 && node[cur_node].child2 == -1 && node[cur_node].child3 == -1 && node[cur_node].child4 == -1)) {
						if(node[cur_node].particle_id > -1) {
							// calculate forces
							//std::cout << "guh" << std::endl;
							a[2 * i + 0] += (dx / d1) * G_SI * node[cur_node].mass * m[i] * (1.0 / d2);
							a[2 * i + 1] += (dy / d1) * G_SI * node[cur_node].mass * m[i] * (1.0 / d2);
						}
					} else {
						// visit child nodes
						// if node exists, add it

						//std::cout << "gah" << std::endl;
						if(node[cur_node].child1 != -1) to_visit.push(node[cur_node].child1);
						if(node[cur_node].child2 != -1) to_visit.push(node[cur_node].child2);
						if(node[cur_node].child3 != -1) to_visit.push(node[cur_node].child3);
						if(node[cur_node].child4 != -1) to_visit.push(node[cur_node].child4);
					}
				}

				//std::cout << "b" << std::endl;
			}
			//std::cout << "c" << std::endl;
		}
};

/*
// exact bias: 0.020888578919738908
uint32_t triple32(uint32_t x) {
	x ^= x >> (uint32_t)17U;
	x *= (uint32_t)0xED5AD4BBU;
	x ^= x >> (uint32_t)11U;
	x *= (uint32_t)0xAC4C1B51U;
	x ^= x >> (uint32_t)15U;
	x *= (uint32_t)0x31848BABU;
	x ^= x >> (uint32_t)14U;
	return x;
}
*/


/*
def W( x, y, z, h ):
	"""
	Gausssian Smoothing kernel (3D)
	x     is a vector/matrix of x positions
	y     is a vector/matrix of y positions
	z     is a vector/matrix of z positions
	h     is the smoothing length
	w     is the evaluated smoothing function
	"""
	
	r = np.sqrt(x**2 + y**2 + z**2)
	
	w = (1.0 / (h*np.sqrt(np.pi)))**3 * np.exp( -r**2 / h**2)
	
	return w
*/

//float W(float x, float y, float z) {
//}

void UpdateAcceleration(unsigned int N, Real* m, Real* r, Real* v, Real* a) {
	for(unsigned int i = 0; i < (2 * N); i++) a[i] = (Real)0;

	// Compute Forces (O(N^2), very slow...)
	{
		for(unsigned int i = 0; i < N; i++) {
			float m1 = m[i];

			//#pragma omp parallel for
			for(unsigned int j = 0; j < N; j++) {
				if(i == j) continue;

				float m2 = m[j];

				float d[2] = {
					{r[2 * j + 0] - r[2 * i + 0]},
					{r[2 * j + 1] - r[2 * i + 1]},
				};

				float r2 = (d[0] * d[0]) + (d[1] * d[1]);
				float r1  = sqrtf(r2);

				a[2 * i + 0] += (d[0] / r1) * G_SI * m1 * m2 * (1.0 / r2);
				a[2 * i + 1] += (d[1] / r1) * G_SI * m1 * m2 * (1.0 / r2);
			}
		}
	}
}

int main(int argc, char** argv) {
	if(argc > 1) {
		std::cerr << "Invalid argument: \"" << argv[1] << "\"" << std::endl;

		return EXIT_FAILURE;
	}

	unsigned int N = 10000; // number of particles

	double dt = 0.01;
	//double dt = 0.66666666666666666666666666666667;
	//double dt = 1.0;

	unsigned int ImageSizeX = 512;
	unsigned int ImageSizeY = 512;

	std::cout
	<< "+---------------+\n"
	<< "|   Starflood   |\n"
	<< "+---------------+\n"
	<< "\n"
	<< "Parameters:\n"
	<< "{\n"
	<< "  N: " << N << "\n"
	<< "  dt: " << dt << "\n"
	<< "  ImageSizeX: " << ImageSizeX << "\n"
	<< "  ImageSizeY: " << ImageSizeY << "\n"
	<< "}\n"
	<< std::endl;

	//return EXIT_SUCCESS;

	Buffer2D<RGBA32F> image(ImageSizeX, ImageSizeY);

	//simulation.Initialize();

	// row-major
	Real* m = (Real*)malloc(sizeof(Real) * (size_t)N * (size_t)1); // Particle masses
	Real* r = (Real*)malloc(sizeof(Real) * (size_t)N * (size_t)2); // Particle positions
	Real* v = (Real*)malloc(sizeof(Real) * (size_t)N * (size_t)2); // Particle velocities
	Real* a = (Real*)malloc(sizeof(Real) * (size_t)N * (size_t)2); // Particle accelerations

	BarnesHut barneshut(-1.0, -1.0, 1.0, 1.0);

	barneshut.N = N;
	barneshut.m = m;
	barneshut.r = r;
	barneshut.v = v;
	barneshut.a = a;

	// initialize particle mass
	//for(unsigned int i = 0; i < N; i++) m[i] = 10.0 * ((2.0 * ((double)std::rand() / (double)RAND_MAX)) - 1.0);
	for(unsigned int i = 0; i < N; i++) m[i] = (Real)((long double)500.0 / (long double)N);

	std::srand(1);

	std::rand();

	// initialize particle position
	for(unsigned int i = 0; i < (2 * N); i++) r[i] = 1.0 * ((2.0 * ((long double)std::rand() / (long double)RAND_MAX)) - 1.0);

	std::srand(2);

	std::rand();

	// initialize particle velocity
	for(unsigned int i = 0; i < (2 * N); i++) v[i] = 0.01 * ((2.0 * ((long double)std::rand() / (long double)RAND_MAX)) - 1.0);
	//for(unsigned int i = 0; i < (2 * N); i++) v[i] = (Real)0;
	//for(unsigned int i = 0; i < (2 * N); i++) v[i] = 0.0;
	//for(unsigned int i = 0; i < (2 * N); i++) {
	//	v[2 * i + 0] =  0.1 * r[2 * i + 1] / sqrtf(r[2 * i + 0]*r[2 * i + 0]+r[2 * i + 1]*r[2 * i + 1]);
	//	v[2 * i + 1] = -001 * r[2 * i + 0] / sqrtf(r[2 * i + 0]*r[2 * i + 0]+r[2 * i + 1]*r[2 * i + 1]);
	//}
	for(unsigned int i = 0; i < (2 * N); i++) {
		v[2 * i + 0] =  0.1 * r[2 * i + 1] / (r[2 * i + 0]*r[2 * i + 0]+r[2 * i + 1]*r[2 * i + 1]);
		v[2 * i + 1] = -0.1 * r[2 * i + 0] / (r[2 * i + 0]*r[2 * i + 0]+r[2 * i + 1]*r[2 * i + 1]);
	}

	//for(unsigned int i = 0; i < (2 * N); i++) a[i] = (Real)0;

	//UpdateAcceleration(N, m, r, v, a); // update acceleration

	std::srand(1); std::rand();
	//for(unsigned int i = 0; i < (2 * N); i++) if(r[i] != r[i]) r[i] = 1.0 * ((2.0 * ((long double)std::rand() / (long double)RAND_MAX)) - 1.0);

	std::srand(2); std::rand();
	//for(unsigned int i = 0; i < (2 * N); i++) if(v[i] != v[i]) v[i] = 0.1 * ((2.0 * ((long double)std::rand() / (long double)RAND_MAX)) - 1.0);

	//for(unsigned int i = 0; i < (2 * N); i++) if(a[i] != a[i]) a[i] = (Real)0;

	Real min_x = r[2 * 0 + 0], max_x = r[2 * 0 + 0], min_y = r[2 * 0 + 1], max_y = r[2 * 0 + 1];

	for(unsigned int i = 1; i < N; i++) {
		min_x = std::min(min_x, r[2 * i + 0]);
		max_x = std::max(max_x, r[2 * i + 0]);
		min_y = std::min(min_y, r[2 * i + 1]);
		max_y = std::max(max_y, r[2 * i + 1]);
	}

	barneshut.bound_min_x = std::max(min_x, (Real)-10.0);
	barneshut.bound_max_x = std::min(max_x, (Real) 10.0);
	barneshut.bound_min_y = std::max(min_y, (Real)-10.0);
	barneshut.bound_max_y = std::min(max_y, (Real) 10.0);

	//std::cout << "d" << std::endl;
	barneshut.Build();

	//std::cout << "e" << std::endl;

	barneshut.UpdateAcceleration();
	//std::cout << "f" << std::endl;

	Real t = 0.0;

	#define FRAME_INTERVAL 4

	for(int n = 1; n < 1200; n++) {
		if(n % FRAME_INTERVAL == 0) {
		//#pragma omp simd
		for(unsigned int i = 0; i < (ImageSizeX * ImageSizeY); i++) image.data[i] = RGBA32F(0.0, 0.0, 0.0, 1.0); // Clear image buffer

		// Render the image
		for(unsigned int i = 0; i < N; i++) {
			// vec2 uv = (fragCoord - 0.5 * iResolution.xy) / min(iResolution.x, iResolution.y);
			// vec2 fragCoord = (min(iResolution.x, iResolution.y) * uv) + (0.5 * iResolution.xy);

			//Eigen::Vector3f p = r.row(i);

			//double _t = 0.0125 * t;

			//Eigen::Vector3f _p(
			//(Real)(p(0) * cos(2.0 * acos(-1.0) * _t)) - (p(2) * sin(2.0 * acos(-1.0) * _t)),
			//(Real)(p(0) * sin(2.0 * acos(-1.0) * _t)) + (p(2) * sin(2.0 * acos(-1.0) * _t)),
			//p(1)
			//);

			//p = _p;

			float px = r[2 * i + 0];
			float py = r[2 * i + 1];

			float uv_x = (fminf(ImageSizeX, ImageSizeY) * (3e-1 * px)) + (0.5 * ImageSizeX);
			float uv_y = (fminf(ImageSizeX, ImageSizeY) * (3e-1 * py)) + (0.5 * ImageSizeY);

			int coord_x = (int)uv_x;
			int coord_y = (int)uv_y;

			if(0 <= coord_x && coord_x < ImageSizeX && 0 <= coord_y && coord_y < ImageSizeY) {
				image.data[coord_x + ImageSizeX * coord_y].r += 0.5;
				image.data[coord_x + ImageSizeX * coord_y].g += 0.5;
				image.data[coord_x + ImageSizeX * coord_y].b += 0.5;
			}
		}

		// Write the image buffer
		{
			char filename[80];

			sprintf(filename, "out/IMG_%04d.tga", n / FRAME_INTERVAL);
			//sprintf(filename, "out/IMG_%04d.hdr", n);

			WriteImageBuffer(image, filename);
			//WriteImageBufferHDR(image, filename);
		}
		}

		//UpdateSimulation(simulation);

		//#pragma omp simd
		for(unsigned int i = 0; i < (2 * N); i++) v[i] += (Real)0.5 * dt * a[i]; // 1/2 kick
		//v = v + (0.5 * dt * a); // 1/2 kick

		//for(unsigned int i = 0; i < (2 * N); i++) v[i] *= 0.99;

		//#pragma omp simd
		for(unsigned int i = 0; i < (2 * N); i++) r[i] += dt * v[i]; // drift
		//r = r + (dt * v); // drift

		//UpdateAcceleration(N, m, r, v, a); // update acceleration

		min_x = r[2 * 0 + 0];
		max_x = r[2 * 0 + 0];
		min_y = r[2 * 0 + 1];
		max_y = r[2 * 0 + 1];

		for(unsigned int i = 1; i < N; i++) {
			min_x = std::min(min_x, r[2 * i + 0]);
			max_x = std::max(max_x, r[2 * i + 0]);
			min_y = std::min(min_y, r[2 * i + 1]);
			max_y = std::max(max_y, r[2 * i + 1]);
		}

		barneshut.bound_min_x = std::max(min_x, (Real)-10.0);
		barneshut.bound_max_x = std::min(max_x, (Real) 10.0);
		barneshut.bound_min_y = std::max(min_y, (Real)-10.0);
		barneshut.bound_max_y = std::min(max_y, (Real) 10.0);

		//std::cout << "d" << std::endl;
		barneshut.Build();

		barneshut.UpdateAcceleration();
		//std::cout << "e" << std::endl;

		// Re-initialize NaN's
		std::srand(1); std::rand();
		//for(unsigned int i = 0; i < (2 * N); i++) if(r[i] != r[i]) r[i] = 1.0 * ((2.0 * ((long double)std::rand() / (long double)RAND_MAX)) - 1.0);

		std::srand(2); std::rand();
		//for(unsigned int i = 0; i < (2 * N); i++) if(v[i] != v[i]) v[i] = 0.1 * ((2.0 * ((long double)std::rand() / (long double)RAND_MAX)) - 1.0);

		//for(unsigned int i = 0; i < (2 * N); i++) if(a[i] != a[i]) a[i] = (Real)0;

		//#pragma omp simd
		for(unsigned int i = 0; i < (2 * N); i++) v[i] += (Real)0.5 * dt * a[i]; // 1/2 kick
		//v = v + (0.5 * dt * a); // 1/2 kick

		//for(unsigned int i = 0; i < (2 * N); i++) v[i] *= 0.99;

		t += dt;
	}

	std::cout << "Freeing allocated memory..." << std::endl;

	free(m);
	free(r);
	free(v);
	free(a);

	std::cout << "Finished!" << std::endl;

	return EXIT_SUCCESS;
}
