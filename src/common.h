#pragma once

#include <cstdint>

// Real number type (currently does not toggle precision of special math functions, e.g. sqrtf!)
typedef float real;
//typedef double real;

// Default render width and height (pixels)
//
// Image data is RGB32F, so 1 frame is 3*4*RENDER_W*RENDER_H bytes!
// Starflood only keeps 1 frame at a time in memory as it writes them
// to the disk
//
// Common Screen Sizes:
//    640× 360:   2,700 KiB/frame
//   1920×1080:  24,300 KiB/frame
//   2560×1440:  43,200 KiB/frame
//   3840×2160:  97,200 KiB/frame
//   7680×4320: 388,800 KiB/frame
// 2^n Sizes:
//   1024×1024:  12 MiB/frame
//   2048×2048:  48 MiB/frame
//   4096×4096: 192 MiB/frame
//   8192×8192: 768 MiB/frame
#define RENDER_W 1920
#define RENDER_H 1080

// Default number of bodies in the simulation
//
#define NUM_BODIES 65536

// Default Barnes-Hut parameter "theta", higher: fewer interactions, less accurate, lower: more interactions, more accurate
//#define THETA 0.5
#define THETA 0.5

// Starflood will render 1 out of every FRAME_INTERVAL frames by default
#define FRAME_INTERVAL 1

// Default number of timesteps
#define FRAMES 600

// Timestep parameter, sometimes called delta t
//#define TIMESTEP 0.01
#define TIMESTEP (1.0*0.0314159265358979323846264338327950288419716939937510582097494459)

// Fits the Barnes-Hut tree to the bodies
#define TREE_FIT

// Randomly jitter the tree, can help stochastically eliminate artifacts
#define JITTER_TREE

// Render Barnes-Hut interactions by default
//#define STARFLOOD_RENDER_INTERACTS

// Enable self-profiling by default
#define STARFLOOD_ENABLE_PROFILING

// Disable writing images (for quick debugging)
//#define STARFLOOD_DISABLE_IMAGE_WRITE

// Starflood version numbering (major.minor.patch)
#ifndef STARFLOOD_VERSION_MAJOR
#define STARFLOOD_VERSION_MAJOR 0
#endif
#ifndef STARFLOOD_VERSION_MINOR
#define STARFLOOD_VERSION_MINOR 0
#endif
#ifndef STARFLOOD_VERSION_PATCH
#define STARFLOOD_VERSION_PATCH 0
#endif

// Toss a few extra digits to ensure optimal floats, see http://www.mimirgames.com/articles/programming/digits-of-pi-needed-for-floating-point-numbers/
const real     PI = 3.1415926535897932384626433832795028841971693993751058209749445923;
const real INV_PI = 0.3183098861837906715377675267450287240689192914809128974953346881;

const real     TAU = 6.2831853071795864769252867665590057683943387987502116419498891846;
const real INV_TAU = 0.1591549430918953357688837633725143620344596457404564487476673441;
