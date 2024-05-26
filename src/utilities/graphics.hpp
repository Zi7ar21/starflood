#pragma once

#include <utilities/common.h>

// Draw a line given pixel coordinates
void drawLine(float* image, int w, int h, float r, float g, float b, float a, int x0, int y0, int x1, int y1);

// Draw a line given screen-space UV coordinates
void drawLineUV(float* image, int w, int h, float r, float g, float b, float a, real x0, real y0, real x1, real y1);

void projection(real* vout, real* vert, real* cam_pos, real* cam_dir, real focal_length, real time);

void drawLine3D(float* image, int w, int h, float r, float g, float b, float a, real* cam_pos, real* cam_dir, real focal_length, real time, real* p0, real* p1);
