#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <stack>
#include <cstdint>
#include <omp.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#define FRAME_INTERVAL 1
#define FRAMES 600
#define TIMESTEP (2.0*0.0314159265358979323846264338327950288419716939937510582097494459)
#define TREE_FIT
#define JITTER_TREE
#define THETA 0.5

// starflood version numbering (major.minor.patch)
#ifndef STARFLOOD_VERSION_MAJOR
#define STARFLOOD_VERSION_MAJOR 0
#endif
#ifndef STARFLOOD_VERSION_MINOR
#define STARFLOOD_VERSION_MINOR 0
#endif
#ifndef STARFLOOD_VERSION_PATCH
#define STARFLOOD_VERSION_PATCH 0
#endif

//#define STARFLOOD_RENDER_INTERACTS

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

void drawLineUV(float* image, int w, int h, float r, float g, float b, float a, real x0, real y0, real x1, real y1) {
	drawLine(image, w, h, r, g, b, a,
	(int)((real)h * x0 + (real)0.5 * (real)w),
	(int)((real)h * y0 + (real)0.5 * (real)h),
	(int)((real)h * x1 + (real)0.5 * (real)w),
	(int)((real)h * y1 + (real)0.5 * (real)h));
}

class Node {
	public:
		/*
		Quadrants:
		+---+---+
		| 2 | 3 |
		+---+---+
		| 0 | 1 |
		+---+---+
		*/

		int parent; // index of this node's parent, -1: this is the root node
		int child0; // index of this node's child in quadrant 0
		int child1; // index of this node's child in quadrant 1
		int child2; // index of this node's child in quadrant 2
		int child3; // index of this node's child in quadrant 3

		int id; // id of the particle contained
		int np; // number of particles contained in this node and children

		real x_min; // this node's axis-aligned bounding-box x-minimum
		real y_min; // this node's axis-aligned bounding-box y-minimum
		real x_max; // this node's axis-aligned bounding-box x-maximum
		real y_max; // this node's axis-aligned bounding-box y-maximum

		real m; // mass contained in this node and children
		real x; // this node's center of mass in the x-axis
		real y; // this node's center of mass in the y-axis

		Node() {}

		Node(int _parent, real _x_min, real _y_min, real _x_max, real _y_max):
		parent(_parent), child0(-1), child1(-1), child2(-1), child3(-1), id(-1), np(0),
		x_min(_x_min), y_min(_y_min), x_max(_x_max), y_max(_y_max), m(0), x(0), y(0) {}
};

void BarnesHut(std::vector<Node> &tree, float* image, int w, int h, int* ids, real* mas, real* pos, real* acc, int N, int step_num) {
	for(int i = 0; i < N; i++) ids[i] = -1;

	tree.clear();

	#ifndef TREE_FIT
	tree.push_back(Node(-1,(real)-1.0,(real)-1.0,(real)1.0,(real)1.0));
	#else
	tree.push_back(Node(-1,pos[0],pos[1],pos[0],pos[1]));
	#endif

	{
		const real extrspace = 0.1;

		#ifdef TREE_FIT
		for(int i = 1; i < N; i++) {
			tree[0].x_min = fmin(tree[0].x_min,pos[2*i+0]);
			tree[0].y_min = fmin(tree[0].y_min,pos[2*i+1]);
			tree[0].x_max = fmax(tree[0].x_max,pos[2*i+0]);
			tree[0].y_max = fmax(tree[0].y_max,pos[2*i+1]);
		}


		tree[0].x_min -= extrspace;
		tree[0].y_min -= extrspace;
		tree[0].x_max += extrspace;
		tree[0].y_max += extrspace;
		#endif

		#ifdef JITTER_TREE
		uint32_t ns = (uint32_t)step_num+(uint32_t)37; // set the random number generator seed

		float z0 = urand(&ns);
		float z1 = urand(&ns);

		tree[0].x_min += extrspace*(real)(z0-0.5f);
		tree[0].y_min += extrspace*(real)(z1-0.5f);
		tree[0].x_max += extrspace*(real)(z0-0.5f);
		tree[0].y_max += extrspace*(real)(z1-0.5f);
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
			real x_max = tree[current_node].x_max;
			real y_max = tree[current_node].y_max;

			// don't continue if particle is outside of bounds
			if((x_max < pos[2*i+0]) || (pos[2*i+0] < x_min) || (y_max < pos[2*i+1]) || (pos[2*i+1] < y_min)) break;

			real hx = ((real)0.5*tree[current_node].x_min)+((real)0.5*tree[current_node].x_max);
			real hy = ((real)0.5*tree[current_node].y_min)+((real)0.5*tree[current_node].y_max);

			if(tree[current_node].child0 != -1
			|| tree[current_node].child1 != -1
			|| tree[current_node].child2 != -1
			|| tree[current_node].child3 != -1) {
				// node is split
				int q = (pos[2*i+0] > hx ? 1 : 0)+(pos[2*i+1] > hy ? 2 : 0);

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
				tree.push_back(Node(current_node,x_min,y_min,   hx,   hy));
				tree[current_node].child0 = tree.size() - 1;

				tree.push_back(Node(current_node,   hx,y_min,x_max,   hy));
				tree[current_node].child1 = tree.size() - 1;

				tree.push_back(Node(current_node,x_min,   hy,   hx,y_max));
				tree[current_node].child2 = tree.size() - 1;

				tree.push_back(Node(current_node,   hx,   hy,x_max,y_max));
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

			int tmpid = current_node;
			while(tmpid != -1) {
				tree[tmpid].m += mas[i];
				tree[tmpid].x += mas[i]*pos[2*i+0];
				tree[tmpid].y += mas[i]*pos[2*i+1];
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
	for(int i = 0; i < tree.size(); i++) {
		tree[i].x = (tree[i].m != (real)0) ? (tree[i].x / tree[i].m) : tree[i].x;
		tree[i].y = (tree[i].m != (real)0) ? (tree[i].y / tree[i].m) : tree[i].y;
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

	// Compute Forces
	#pragma omp parallel for schedule(dynamic,256)
	for(int i = 0; i < N; i++) {
		// don't continue if particle is outside of bounds
		if((tree[0].x_max < pos[2*i+0]) || (pos[2*i+0] < tree[0].x_min) || (tree[0].y_max < pos[2*i+1]) || (pos[2*i+1] < tree[0].y_min)) continue;

		acc[2*i+0] = (real)0;
		acc[2*i+1] = (real)0;

		real u_i = pos[2*i+0];
		real v_i = pos[2*i+1];

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

			real x_min = tree[cur_node].x_min, y_min = tree[cur_node].y_min, x_max = tree[cur_node].x_max, y_max = tree[cur_node].y_max;

			#ifdef STARFLOOD_RENDER_INTERACTS
			if(i == 0) {
				drawLineUV(image, w, h, 1., 1., 0., 1., x_min, y_min, x_max, y_min);
				drawLineUV(image, w, h, 1., 1., 0., 1., x_min, y_min, x_min, y_max);
				drawLineUV(image, w, h, 1., 1., 0., 1., x_min, y_max, x_max, y_max);
				drawLineUV(image, w, h, 1., 1., 0., 1., x_max, y_min, x_max, y_max);
			}
			#endif

			real width = ((real)0.5*(x_max - x_min))+((real)0.5*(y_max - y_min));

			real dx = tree[cur_node].x;
			real dy = tree[cur_node].y;
			real dm = tree[cur_node].m;

			// probably empty, idk just to be safe
			if(dm < (real)0.001) continue;

			if(x_min < u_i && u_i < x_max && y_min < v_i && v_i < y_max) {
				// this particle is contained in the node we are evaluating so let's
				// calculate what it would be like with it removed
				dx *= dm;
				dy *= dm;
				dm -= mas[i];
				dx -= mas[i]*pos[2*i+0];
				dy -= mas[i]*pos[2*i+1];
				dx /= dm;
				dy /= dm;
			}

			dx = dx-u_i;
			dy = dy-v_i;

			real dist = sqrtf((dx*dx)+(dy*dy));

			// again, probably empty, just to be safe
			if(dm < (real)0.001) continue;

			const real theta = THETA;

			if((width/dist) > theta) {
				if(tree[cur_node].child0 != -1) to_visit.push(tree[cur_node].child0);
				if(tree[cur_node].child1 != -1) to_visit.push(tree[cur_node].child1);
				if(tree[cur_node].child2 != -1) to_visit.push(tree[cur_node].child2);
				if(tree[cur_node].child3 != -1) to_visit.push(tree[cur_node].child3);
				continue;
			}

			#ifdef STARFLOOD_RENDER_INTERACTS
			if(i == 0) {
				//drawLineUV(image, w, h, 1., 1., 0., 1., x_min, y_min, x_max, y_min);
				//drawLineUV(image, w, h, 1., 1., 0., 1., x_min, y_min, x_min, y_max);
				//drawLineUV(image, w, h, 1., 1., 0., 1., x_min, y_max, x_max, y_max);
				//drawLineUV(image, w, h, 1., 1., 0., 1., x_max, y_min, x_max, y_max);
				drawLineUV(image, w, h, 1., 0., 0., 1., u_i, v_i, tree[cur_node].x, tree[cur_node].y);
			}
			#endif

			acc[2*i+0] += ((dx/dist)*mas[i]*dm)/((dx*dx)+(dy*dy)+(real)0.0001);
			acc[2*i+1] += ((dy/dist)*mas[i]*dm)/((dx*dx)+(dy*dy)+(real)0.0001);

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

int main(int argc, char** argv) {
	if(argc > 1) {
		if(strcmp(argv[1],"-h") || strcmp(argv[1],"--help")) {
			printf(
			"Usage: %s [OPTION]...\n"
			"  -h, --help       display this help and exit\n"
			//"      --version    output version information and exit\n"
			);

			return EXIT_SUCCESS;
		} else {
			printf("%s: unrecognized option '%s'\nTry '%s -h' or '%s --help' for more information.\n", argv[0], argv[1]);

			return EXIT_FAILURE;
		}
	}

	printf("Starflood %d.%d.%d\n", STARFLOOD_VERSION_MAJOR, STARFLOOD_VERSION_MINOR, STARFLOOD_VERSION_PATCH); // print version information

	#ifndef STARFLOOD_NO_SPLASH
	printf(
	"\n"
	"╭───┬─╴╭─╮┌─╮┌──┐  ╭─╮╭─╮┌─╮\n"
	"│   │  │ ││ ││  │  │ ││ ││ │\n"
	"╰─╮ │  ├─┤├─╯├─╴│  │ ││ ││ │\n"
	"│ │ │  │ ││ ││  │  │ ││ ││ │\n"
	"╰─╯ ╵  ╵ ╵╵ ╵╵  └──╰─╯╰─╯└─╯\n"
	"\n");
	#endif

	fflush(stdout);

	// timestamps
	double t0, t1;

	// call omp_get_wtime() so it is cached or something idk
	t0 = omp_get_wtime();
	t1 = omp_get_wtime();

	size_t w = 1920; // image width
	size_t h = 1080; // image height
	size_t N = 2*131072; // number of bodies in the simulation

	if(N < 1) {
		printf("Error: N was set to %zu, it must be greater than 0!\n", N);

		return EXIT_FAILURE;
	}

	const real dt = TIMESTEP; // Simulation Timestep

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

	std::vector<Node> tree;

	// Allocate Memory
	{
		image = (float*)malloc(image_size); // Allocate image buffer

		if(image == NULL) {
			printf("Error: malloc() returned a null pointer!\n");

			fflush(stdout);

			perror("Error");

			return EXIT_FAILURE;
		}

		sim = (real*)malloc(sim_size); // Allocate simulation memory

		if(sim == NULL) {
			printf("Error: malloc() returned a null pointer!\n");

			fflush(stdout);

			perror("Error");

			free(image);

			return EXIT_FAILURE;
		}

		ids = (int*)malloc(ids_size);

		mas = sim; // mass buffer
		pos = mas+(mas_size/sizeof(real)); // position buffer
		vel = pos+(pos_size/sizeof(real)); // velocity buffer
		acc = vel+(vel_size/sizeof(real)); // acceleration buffer
	}

	// Initialize Simulation
	{
		for(int i = 0; i < (N * 1); i++) mas[i] = (real)(((long double)10)/((long double)N));
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
			//pos[2*i+0] = 0.1*sqrt(-2.0*log(z0))*cos(TAU*z2);
			//pos[2*i+1] = 0.1*sqrt(-2.0*log(z0))*sin(TAU*z2);
			//pos[2*i+2] = 0.1*sqrt(-2.0*log(z1))*cos(TAU*z3);
			pos[2*i+0] = 0.2*sqrt(z1)*cos(TAU*z0);
			pos[2*i+1] = 0.2*sqrt(z1)*sin(TAU*z0);
		}

		for(int i = 0; i < N; i++) {
			vel[2*i+0] = -.3*pos[2*i+1];
			vel[2*i+1] =  .3*pos[2*i+0];
		}

		// calculate initial gravitational accelerations
		BarnesHut(tree, image, w, h, ids, mas, pos, acc, N, 0);
	}

	#ifdef STARFLOOD_ENABLE_PROFILING
	FILE* diagfile = fopen("./log.csv", "w");

	fprintf(diagfile, "\"n\",\"Clear Image (ms)\",\"Render Bodies (ms)\",\"Render Quadtree (ms)\",\"Write Image (ms)\",\"First Kick (ms)\",\"Drift (ms)\",\"Update Acceleration (ms)\",\"Second Kick (ms)\"\n");

	fflush(diagfile);
	#endif

	const int num_steps = FRAMES;

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

				if((step_num % FRAME_INTERVAL) == 0) {
					// VERY FAST
					//#pragma omp simd
					for(int i = 0; i < (4*w*h); i++) image[i] = (float)0;
				}

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

				if((step_num % FRAME_INTERVAL) == 0){
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
							image[4*(w*coord[1]+coord[0])+0] += 0.05f;
							image[4*(w*coord[1]+coord[0])+1] += 0.05f;
							image[4*(w*coord[1]+coord[0])+2] += 0.05f;
							image[4*(w*coord[1]+coord[0])+3] += 0.0f;
						}
					}

					for(int i = 0; i < (w * h); i++) {
						image[4*i+3] = (real)1;
					}
				}

				#ifdef STARFLOOD_ENABLE_PROFILING
				t1 = omp_get_wtime();
				fprintf(diagfile, "%.6f,", 1000.0*(t1-t0));
				#endif
			}

			{
				#ifdef STARFLOOD_ENABLE_PROFILING
				t0 = omp_get_wtime();
				#endif

				if((step_num % FRAME_INTERVAL) == 0) {
					/*
					for(int i = 0; i < tree.size(); i++) {
						real x_min = tree[i].x_min, y_min = tree[i].y_min, x_max = tree[i].x_max, y_max = tree[i].y_max;

						if((x_max-x_min) < (real)0.1) continue;

						drawLineUV(image, w, h, 0., .01, 0., .1, x_min, y_min, x_max, y_min);
						drawLineUV(image, w, h, 0., .01, 0., .1, x_min, y_min, x_min, y_max);
						drawLineUV(image, w, h, 0., .01, 0., .1, x_min, y_max, x_max, y_max);
						drawLineUV(image, w, h, 0., .01, 0., .1, x_max, y_min, x_max, y_max);
					}
					*/
				}

				#ifdef STARFLOOD_ENABLE_PROFILING
				t1 = omp_get_wtime();
				fprintf(diagfile, "%.6f,", 1000.0*(t1-t0));
				#endif
			}

			if((step_num % FRAME_INTERVAL) == 0){
				#ifdef STARFLOOD_RENDER_INTERACTS
				BarnesHut(tree, image, w, h, ids, mas, pos, acc, N, step_num);
				#endif

				// test, test, is this thing working?
				drawLine(image, w, h, 1., 0., 0., 1., 15+0, 15-1, 15+4, 15-1); drawLine(image, w, h, 1., 0., 0., 1., 15+0, 15+1, 15+4, 15+1);
				drawLine(image, w, h, 1., 1., 0., 1., 15-1, 15+0, 15-1, 15+4); drawLine(image, w, h, 1., 1., 0., 1., 15+1, 15+0, 15+1, 15+4);
				drawLine(image, w, h, 0., 1., 0., 1., 15+0, 15-1, 15-4, 15-1); drawLine(image, w, h, 0., 1., 0., 1., 15+0, 15+1, 15-4, 15+1);
				drawLine(image, w, h, 0., 0., 1., 1., 15-1, 15+0, 15-1, 15-4); drawLine(image, w, h, 0., 0., 1., 1., 15+1, 15+0, 15+1, 15-4);

				for(int i = 0; i < (w * h); i++) {
					image[4*i+0] /= (image[4*i+3] != (float)0) ? image[4*i+3] : (float)1;
					image[4*i+1] /= (image[4*i+3] != (float)0) ? image[4*i+3] : (float)1;
					image[4*i+2] /= (image[4*i+3] != (float)0) ? image[4*i+3] : (float)1;
					image[4*i+3]  = (float)1;
				}
			}

			// write the image
			{
				char file_name[FILENAME_MAX];

				snprintf(file_name, FILENAME_MAX, "out/step_%04d.hdr", step_num/FRAME_INTERVAL);

				#ifdef STARFLOOD_ENABLE_PROFILING
				t0 = omp_get_wtime();
				#endif

				if((step_num % FRAME_INTERVAL) == 0){
					if(stbi_write_hdr(file_name, w, h, 4, image) == 0) {
						printf("Failure on step %d\n", step_num);
					}
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

		// update acceleration
		{
			for(int i = 0; i < (N * 2); i++) acc[i] = (real)0;

			#ifdef STARFLOOD_ENABLE_PROFILING
			t0 = omp_get_wtime();
			#endif

			BarnesHut(tree, image, w, h, ids, mas, pos, acc, N, step_num);

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

		/*
		for(int i = 0; i < N; i++) {
			if(abs(pos[2*i+0] > 1.) || abs(pos[2*i+0]) > 1.) {
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
			}
		}
		*/
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