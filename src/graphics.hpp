#pragma once

#include <common.h>

// Draw a line given pixel coordinates
void drawLine(float* image, int w, int h, float r, float g, float b, float a, int x0, int y0, int x1, int y1);

// Draw a line given screen-space UV coordinates
void drawLineUV(float* image, int w, int h, float r, float g, float b, float a, real x0, real y0, real x1, real y1);
