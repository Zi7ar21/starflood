#include <utilities/graphics.hpp>



#include <cmath>
//#include <omp.h>

// Draw a line given pixel coordinates
void drawLine(float* image, int w, int h, float r, float g, float b, float a, int x0, int y0, int x1, int y1) {
	// Bresenham's Line Algorithm

	int dx = abs(x1 - x0);
	int sx = x0 < x1 ? 1 : -1;
	int dy = -abs(y1 - y0);
	int sy = y0 < y1 ? 1 : -1;
	int error = dx + dy;

	while(true) {
		if((0 <= x0) && (x0 < w) && (0 <= y0) && (y0 < h)) {
			int index = w*y0+x0;
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

// Draw a line given screen-space UV coordinates
void drawLineUV(float* image, int w, int h, float r, float g, float b, float a, real x0, real y0, real x1, real y1) {
	drawLine(image, w, h, r, g, b, a,
	(int)((real)h * x0 + (real)0.5 * (real)w),
	(int)((real)h * y0 + (real)0.5 * (real)h),
	(int)((real)h * x1 + (real)0.5 * (real)w),
	(int)((real)h * y1 + (real)0.5 * (real)h));
}

// vout[3]: projected point, vert[3] point to project
void projection(real* vout, real* vert, real* cam_pos, real* cam_dir, real focal_length, real time) {
	const double rotation_theta = -0.078125*PI;

	//real rcam[3] = {
	//	cam_dir[0],
	//	cam_dir[1]*(real)cos(rotation_theta)-cam_dir[2]*(real)sin(rotation_theta),
	//	cam_dir[1]*(real)sin(rotation_theta)+cam_dir[2]*(real)cos(rotation_theta)
	//};


	real diff[3] = {
		vert[0] - cam_pos[0],
		vert[1] - cam_pos[1],
		vert[2] - cam_pos[2]
	};

	{
		real rdif[3] = {
			diff[0],
			diff[1]*(real)cos(rotation_theta)-diff[2]*(real)sin(rotation_theta),
			diff[1]*(real)sin(rotation_theta)+diff[2]*(real)cos(rotation_theta)
		};

		diff[0] = rdif[0];
		diff[1] = rdif[1];
		diff[2] = rdif[2];
	}

	vout[0] = (real)( 0.0);
	vout[1] = (real)( 0.0);
	vout[2] = (real)(-1.0);

	// clip particles behind the camera
	if((cam_dir[0]*diff[0])+(cam_dir[1]*diff[1])+(cam_dir[2]*diff[2]) < (real)0.0) return;

	real dist = (real)sqrt((diff[0]*diff[0])+(diff[1]*diff[1])+(diff[2]*diff[2]));

	/*
	real obj_dif[3] = {
	pos[3*i+0]-cam_pos[0],
	pos[3*i+1]-cam_pos[1],
	pos[3*i+2]-cam_pos[2]};
	*/

	//vout[0] = (diff[0]*(real)1.)/(dist*focal_length);
	//vout[1] = (diff[1]*(real)1.)/(dist*focal_length);
	vout[0] = (diff[0]*(real)1.)/(diff[2]*focal_length);
	vout[1] = (diff[1]*(real)1.)/(diff[2]*focal_length);
	vout[2] = dist;
}

void drawLine3D(float* image, int w, int h, float r, float g, float b, float a, real* cam_pos, real* cam_dir, real focal_length, real time, real* p0, real* p1) {
	real proj0[3];
	real proj1[3];

	p0[2] = fmaxf(p0[2],cam_pos[2]);
	p1[2] = fmaxf(p1[2],cam_pos[2]);

	projection(proj0, p0, cam_pos, cam_dir, focal_length, time);
	projection(proj1, p1, cam_pos, cam_dir, focal_length, time);

	if((proj0[2] <= (real)0.0) && (proj1[2] <= (real)0.0)) return; // skip drawing if either point is behind camera

	drawLineUV(image, w, h, r, g, b, a, proj0[0], proj0[1], proj1[0], proj1[1]);
}

