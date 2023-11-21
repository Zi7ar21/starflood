#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <stack>
#include <cstdint>
#include <omp.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STARFLOOD_ENABLE_PROFILING

typedef float real;

const real     PI = 3.1415926535897932384626433832795028841971693993751058209749445923;
const real INV_PI = 0.3183098861837906715377675267450287240689192914809128974953346881;

const real     TAU = 6.2831853071795864769252867665590057683943387987502116419498891846;
const real INV_TAU = 0.1591549430918953357688837633725143620344596457404564487476673441;

// https://nullprogram.com/blog/2018/07/31/
uint32_t triple32(uint32_t x) {
	// exact bias: 0.020888578919738908
	x ^= x >> (uint32_t)17u;
	x *= (uint32_t)0xED5AD4BBu;
	x ^= x >> (uint32_t)11u;
	x *= (uint32_t)0xAC4C1B51u;
	x ^= x >> (uint32_t)15u;
	x *= (uint32_t)0x31848BABu;
	x ^= x >> (uint32_t)14u;
	return x;
}

float urand(uint32_t* ns) {
	uint32_t _ns = *ns;

	_ns = triple32(_ns);

	*ns = _ns;

	return (2.32830629776081821092e-10f)*(float)_ns;
}

void drawLine(float* image, int w, int h, float r, float g, float b, float a, int x0, int y0, int x1, int y1) {
	// Bresenham's Line Algorithm

	int dx = abs(x1 - x0);
	int sx = x0 < x1 ? 1 : -1;
	int dy = -abs(y1 - y0);
	int sy = y0 < y1 ? 1 : -1;
	int error = dx + dy;

	while(true) {
		if((0 <= x0) && (x0 < w) && (0 <= y0) && (y0 < h)) {
			int index = (w*y0)+x0;
			image[4*index+0] += r;
			image[4*index+1] += g;
			image[4*index+2] += b;
			image[4*index+3] += a;
		}

		if(x0 == x1 && y0 == y1) break;

		int e2 = 2 * error;

		if(e2 >= dy) {
			if(x0 == x1) break;
			error = error + dy;
			x0 = x0 + sx;
		}

		if(e2 <= dx) {
			if(y0 == y1) break;
			error = error + dx;
			y0 = y0 + sy;
		}
	}
}

class Node {
	public:
		int parent;
		int child0;
		int child1;
		int child2;
		int child3;

		int id;
		int np;
		float m;
		float x;
		float y;

		Node(): parent(-1), child0(-1), child1(-1), child2(-1), child3(-1), id(-1), np(0), m(0.0f), x(0.0f), y(0.0f) {}

		Node(int _parent): parent(_parent), child0(-1), child1(-1), child2(-1), child3(-1), id(-1), np(0), m(0.0f), x(0.0f), y(0.0f) {}
};

/*
void Remove() {
}
*/

void drawLineUV(float* image, int w, int h, float r, float g, float b, float a, real x0, real y0, real x1, real y1) {
	drawLine(image, w, h, r, g, b, a,
	(int)((real)h * x0 + (real)0.5 * (real)w),
	(int)((real)h * y0 + (real)0.5 * (real)h),
	(int)((real)h * x1 + (real)0.5 * (real)w),
	(int)((real)h * y1 + (real)0.5 * (real)h));
}

void BarnesHut(std::vector<Node> tree, float* image, int w, int h, int* ids, real* pos, real* acc, real bx_min, real by_min, real bx_max, real by_max, int N) {
	for(int i = 0; i < N; i++) ids[i] = -1;

	tree.clear();

	tree.push_back(-1);

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

		real min_x = bx_min;
		real min_y = by_min;
		real max_x = bx_max;
		real max_y = by_max;

		int depth = 1;

		while(depth <= 12) {
			real hx = ((real)0.5*min_x)+((real)0.5*max_x);
			real hy = ((real)0.5*min_y)+((real)0.5*max_y);

			if(tree[current_node].child0 != -1
			|| tree[current_node].child1 != -1
			|| tree[current_node].child2 != -1
			|| tree[current_node].child3 != -1) {
				// node is split
				int q = (pos[2*i+0] > hx ? 1 : 0)+(pos[2*i+1] > hy ? 2 : 0);

				switch(q) {
					case 0:
					max_x = hx;
					max_y = hy;
					current_node = tree[current_node].child0;
					break;
					case 1:
					min_x = hx;
					max_y = hy;
					current_node = tree[current_node].child1;
					break;
					case 2:
					max_x = hx;
					min_y = hy;
					current_node = tree[current_node].child2;
					break;
					case 3:
					min_x = hx;
					min_y = hy;
					current_node = tree[current_node].child3;
					break;
					default:
					printf("Severe error!\n");
					break;
				}

				depth++;

				continue;
			} 

			//if(depth < 6) {
			if(true) {
				drawLineUV(image, w, h, 0., .5, 0., 1., min_x, min_y, max_x, min_y);
				drawLineUV(image, w, h, 0., .5, 0., 1., min_x, min_y, min_x, max_y);
				drawLineUV(image, w, h, 0., .5, 0., 1., min_x, max_y, max_x, max_y);
				drawLineUV(image, w, h, 0., .5, 0., 1., max_x, min_y, max_x, max_y);
				drawLineUV(image, w, h, 0., .5, 0., 1., min_x,    hy, max_x,    hy);
				drawLineUV(image, w, h, 0., .5, 0., 1.,    hx, min_y,    hx, max_y);
			}

			// node is not split
			if(tree[current_node].id != -1) {
				// node has particle
				tree.push_back(current_node);
				tree[current_node].child0 = tree.size() - 1;

				tree.push_back(current_node);
				tree[current_node].child1 = tree.size() - 1;

				tree.push_back(current_node);
				tree[current_node].child2 = tree.size() - 1;

				tree.push_back(current_node);
				tree[current_node].child3 = tree.size() - 1;

				int id = tree[current_node].id;

				tree[current_node].id = -1;

				int q = (pos[2*id+0] > hx ? 1 : 0)+(pos[2*id+1] > hy ? 2 : 0);

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
					default:
					printf("Severe error!\n");
					break;
				}

				continue;
			}

			// node doesn't have particle
			tree[current_node].np = 1;
			tree[current_node].id = i;
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
	for(int i = 0; i < tree.size(); i++) {
		tree[i].x = (tree[i].m != (real)0) ? (tree[i].x / tree[i].m) : tree[i].x;
		tree[i].y = (tree[i].m != (real)0) ? (tree[i].y / tree[i].m) : tree[i].y;
	}

	// Compute Forces
	for(int i = 0; i < N; i++) {
		std::stack<int> to_visit;

		to_visit.push(0);

		int depth = 0;

		real x_min = bx_min;
		real y_min = by_min;//width/dist;
		real x_max = bx_max;
		real y_max = by_max;

		while(!to_visit.empty()) {
			real hx = ((real)0.5*x_min)+((real)0.5*x_max);
			real hy = ((real)0.5*x_min)+((real)0.5*x_max);

			int cur_node = to_visit.top();

			to_visit.pop();

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

	// is this thing working?
	drawLine(image, w, h, 1., 0., 0., 1., 15+0, 15-1, 15+4, 15-1); drawLine(image, w, h, 1., 0., 0., 1., 15+0, 15+1, 15+4, 15+1);
	drawLine(image, w, h, 1., 1., 0., 1., 15-1, 15+0, 15-1, 15+4); drawLine(image, w, h, 1., 1., 0., 1., 15+1, 15+0, 15+1, 15+4);
	drawLine(image, w, h, 0., 1., 0., 1., 15+0, 15-1, 15-4, 15-1); drawLine(image, w, h, 0., 1., 0., 1., 15+0, 15+1, 15-4, 15+1);
	drawLine(image, w, h, 0., 0., 1., 1., 15-1, 15+0, 15-1, 15-4); drawLine(image, w, h, 0., 0., 1., 1., 15+1, 15+0, 15+1, 15-4);
}

int main(int argc, char** argv) {
	if(argc > 1) {
		printf("%s: Invalid Parameter '%s'\n", argv[0], argv[1]);

		return EXIT_FAILURE;
	}

	printf("Hello, world!\n");

	fflush(stdout);

	// timestamps
	double t0, t1;

	// initialize timestamps so it is cached or something idk
	t0 = omp_get_wtime();
	t1 = omp_get_wtime();
	int w = 960; // Image Width
	int h = 540; // Image Height
	const long long int N = 384; // number of bodies in the simulation

	const real dt = 0.0314159265358979323846264338327950288419716939937510582097494459; // Simulation Timestep

	float* image; // Image Buffer, RGBA32F
	real* sim; // Simulation Buffer
	real* mas; // Mass Buffer
	real* pos; // Position Buffer
	real* vel; // Velocity Buffer
	real* acc; // Acceleration Buffer
	int* ids; // id buffer, contains ids of node paticles are in

	stbi_flip_vertically_on_write(1);
	size_t image_size = sizeof(float)*(4*w*h); // Image Buffer Size, RGBA32F
	size_t mas_size = (N * 1) * sizeof(real); // mass buffer size
	size_t pos_size = (N * 2) * sizeof(real); // position buffer size
	size_t vel_size = (N * 2) * sizeof(real); // velocity buffer size
	size_t acc_size = (N * 2) * sizeof(real); // acceleration buffer size
	size_t sim_size = mas_size+pos_size+vel_size+acc_size; // simulation total buffer size
	size_t ids_size = (N * 1) * sizeof(int);

	printf("%zu\n",sim_size);

	std::vector<Node> tree;

	// Allocate Memory
	{
		image = (float*)malloc(image_size); // Allocate image buffer

		if(image == NULL) {
			printf("Error: malloc() returned a null pointer!\n");

			fflush(stdout);

			return EXIT_FAILURE;
		}

		sim = (real*)malloc(sim_size); // Allocate simulation memory

		if(sim == NULL) {
			printf("Error: malloc() returned a null pointer!\n");

			fflush(stdout);

			free(image);

			return EXIT_FAILURE;
		}

		ids = (int*)malloc(ids_size);

		mas = sim; // mass buffer
		pos = mas+(mas_size/sizeof(real)); // position buffer
		vel = pos+(pos_size/sizeof(real)); // velocity buffer
		acc = vel+(vel_size/sizeof(real)); // acceleration buffer

		printf("mas: %p\n",mas);
		printf("pos: %p\n",pos);
		printf("vel: %p\n",vel);
		printf("acc: %p\n",acc);
	}

	real tmin_x = (real)-1.;
	real tmin_y = (real)-1.;
	real tmax_x = (real)1.;
	real tmax_y = (real)1.;

	{

	}

	// Initialize Simulation
	{
		for(int i = 0; i < (N * 1); i++) mas[i] = (real)1;
		for(int i = 0; i < (N * 2); i++) acc[i] = (real)0;
		for(int i = 0; i < (N * 2); i++) vel[i] = (real)0;

		// Generate Initial Conditions
		uint32_t ns = (uint32_t)1u;

		// randomly selected positions and velocities
		for(int i = 0; i < N; i++) {
			ns = (uint32_t)i+(uint32_t)42u; // set the random number generator seed

			float z0 = urand(&ns);
			float z1 = urand(&ns);
			float z2 = urand(&ns);
			float z3 = urand(&ns);

			// normal distribution
			pos[2*i+0] = 0.1*sqrt(-2.0*log(z0))*cos(TAU*z2);
			pos[2*i+1] = 0.1*sqrt(-2.0*log(z0))*sin(TAU*z2);
			//pos[2*i+2] = 0.1*sqrt(-2.0*log(z1))*cos(TAU*z3);
		}

		for(int i = 0; i < N; i++) {
			vel[2*i+0] = -.1*pos[2*i+1];
			vel[2*i+1] =  .1*pos[2*i+0];
		}

		// calculate initial gravitational accelerations
		BarnesHut(tree, image, w, h, ids, pos, acc, tmin_x, tmin_y, tmax_x, tmax_y, N);
	}

	#ifdef STARFLOOD_ENABLE_PROFILING
	FILE* diagfile = fopen("./log.csv", "w");

	fprintf(diagfile, "\"n\",\"Clear Image (ms)\",\"Render Image (ms)\",\"Write Image (ms)\",\"First Kick (ms)\",\"Drift (ms)\",\"Update Acceleration (ms)\",\"Second Kick (ms)\"\n");

	fflush(diagfile);
	#endif

	const int num_steps = 150;

	for(int step_num = 0; step_num < num_steps; step_num++) {
		printf("\rRunning simulation, %d/%d (%.2f%) completed...", step_num, num_steps, 100.0*((double)step_num/(double)num_steps));

		fflush(stdout);

		#ifdef STARFLOOD_ENABLE_PROFILING
		fprintf(diagfile, "%d,", step_num);
		#endif

		// render and write the image
		{
			stbi_flip_vertically_on_write(1);

			// clear the image buffer
			{
				#ifdef STARFLOOD_ENABLE_PROFILING
				t0 = omp_get_wtime();
				#endif

				// VERY FAST
				//#pragma omp simd
				for(int i = 0; i < (4*w*h); i++) image[i] = (float)0;
			
				#ifdef STARFLOOD_ENABLE_PROFILING
				t1 = omp_get_wtime();
				fprintf(diagfile, "%.6f,", 1000.0*(t1-t0));
				#endif
			}

			// render the image
			{
				#ifdef STARFLOOD_ENABLE_PROFILING
				t0 = omp_get_wtime();
				#endif
			
				for(int i = 0; i < N; i++) {
					// vec2 uv = (fragCoord - 0.5 * resolution) / resolution.y
					real t = 0.04*(real)step_num;
					//t = 0.0;
					//real uv[2] = {r[2 * i + 0]*cos(t)+r[2 * i + 2]*sin(t), r[2 * i + 1]};
					real uv[2] = {pos[2*i+0], pos[2*i+1]};

					uv[0] *= (real)1.0; // 0.05
					uv[1] *= (real)1.0;

					int coord[2] = {
					(int)((real)h * uv[0] + (real)0.5 * (real)w),
					(int)((real)h * uv[1] + (real)0.5 * (real)h)};

					if((0 <= coord[0]) && (coord[0] < w) && (0 <= coord[1]) && (coord[1] < h)) {
						image[4*(w*coord[1]+coord[0])+0] += 0.5f;
						image[4*(w*coord[1]+coord[0])+1] += 0.5f;
						image[4*(w*coord[1]+coord[0])+2] += 0.5f;
						image[4*(w*coord[1]+coord[0])+3] += 1.0f;
					}
				}

				#ifdef STARFLOOD_ENABLE_PROFILING
				t1 = omp_get_wtime();
				fprintf(diagfile, "%.6f,", 1000.0*(t1-t0));
				#endif
			}

			BarnesHut(tree, image, w, h, ids, pos, acc, tmin_x, tmin_y, tmax_x, tmax_y, N); // if you want tree rendered

			for(int i = 0; i < (w * h); i++) {
				image[4*i+0] /= (image[4*i+3] != (float)0) ? image[4*i+3] : (float)1;
				image[4*i+1] /= (image[4*i+3] != (float)0) ? image[4*i+3] : (float)1;
				image[4*i+2] /= (image[4*i+3] != (float)0) ? image[4*i+3] : (float)1;
				image[4*i+3]  = (float)1;
			}

			// write the image
			{
				char file_name[FILENAME_MAX];

				snprintf(file_name, FILENAME_MAX, "out/step_%04d.hdr", step_num);
				#ifdef STARFLOOD_ENABLE_PROFILING
				t0 = omp_get_wtime();
				#endif

				if(stbi_write_hdr(file_name, w, h, 4, image) == 0) {
					printf("Failure on step %d\n", step_num);
				}

				#ifdef STARFLOOD_ENABLE_PROFILING
				t1 = omp_get_wtime();
				fprintf(diagfile, "%.6f,", 1000.0*(t1-t0));
				#endif
			}
		}

		// kick
		{
			#ifdef STARFLOOD_ENABLE_PROFILING
			t0 = omp_get_wtime();
			#endif

			// VERY FAST
			//#pragma omp simd
			for(int i = 0; i < (N * 2); i++) vel[i] += (real)0.5 * dt * acc[i];

			#ifdef STARFLOOD_ENABLE_PROFILING
			t1 = omp_get_wtime();
			fprintf(diagfile, "%.6f,", 1000.0*(t1-t0));
			#endif
		}

		// drift
		{
			#ifdef STARFLOOD_ENABLE_PROFILING
			t0 = omp_get_wtime();
			#endif

			// VERY FAST
			//#pragma omp simd
			for(int i = 0; i < (N * 2); i++) pos[i] += dt * vel[i];

			#ifdef STARFLOOD_ENABLE_PROFILING
			t1 = omp_get_wtime();
			fprintf(diagfile, "%.6f,", 1000.0*(t1-t0));
			#endif
		}

		if(step_num == 10) printf("%f\n", acc[30]);

		// update acceleration
		{
			for(int i = 0; i < (N * 2); i++) acc[i] = (real)0;

			#ifdef STARFLOOD_ENABLE_PROFILING
			t0 = omp_get_wtime();
			#endif

			BarnesHut(tree, image, w, h, ids, pos, acc, tmin_x, tmin_y, tmax_x, tmax_y, N);

			#ifdef STARFLOOD_ENABLE_PROFILING
			t1 = omp_get_wtime();
			fprintf(diagfile, "%.6f,", 1000.0*(t1-t0));
			#endif
		}

		//for(int i = 0; i < N; i++) {
		//	a[2*i+0] += -1.0*r[2*i+2];
		//	a[2*i+2] +=  1.0*r[2*i+0];
		//}

		// kick
		{
			#ifdef STARFLOOD_ENABLE_PROFILING
			t0 = omp_get_wtime();
			#endif

			// VERY FAST
			//#pragma omp simd
			for(int i = 0; i < (N * 2); i++) vel[i] += (real)0.5 * dt * acc[i];

			#ifdef STARFLOOD_ENABLE_PROFILING
			t1 = omp_get_wtime();
			fprintf(diagfile, "%.6f\n", 1000.0*(t1-t0));
			#endif
		}
	}

	printf("\rRunning simulation, %d/%d (100.00%) completed... Done!\n\n", num_steps, num_steps);

	fflush(stdout);

	#ifdef STARFLOOD_ENABLE_PROFILING
	fclose(diagfile);
	#endif

	free(image);
	free(sim);

	return EXIT_SUCCESS;
}