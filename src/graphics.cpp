#include <graphics.hpp>

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

// Draw a line given screen-space UV coordinates
void drawLineUV(float* image, int w, int h, float r, float g, float b, float a, real x0, real y0, real x1, real y1) {
	drawLine(image, w, h, r, g, b, a,
	(int)((real)h * x0 + (real)0.5 * (real)w),
	(int)((real)h * y0 + (real)0.5 * (real)h),
	(int)((real)h * x1 + (real)0.5 * (real)w),
	(int)((real)h * y1 + (real)0.5 * (real)h));
}
