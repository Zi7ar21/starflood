#pragma once

// Mathematical Constants

// 64 decimal places (exlcuding first digit)

#define      PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define  INV_PI 0.3183098861837906715377675267450287240689192914809128974953346881
#define SQRT_PI 1.7724538509055160272981674833411451827975494561223871282138077899

#define      TAU 6.2831853071795864769252867665590057683943387987502116419498891846
#define  INV_TAU 0.1591549430918953357688837633725143620344596457404564487476673441
#define SQRT_TAU 2.5066282746310005024157652848110452530069867406099383166299235763

// Physical Constants

// 2018 CODATA-recommended value: 6.67430(15)*10^-11 m^3 kg^-1 s^-2
//#define G_SI 6.6743015e-11
#define G_SI 1.0

// ##### Constants #####

/*
https://float.exposed/
http://www.mimirgames.com/articles/programming/digits-of-pi-needed-for-floating-point-numbers/
*/

/*
It's convenient to make constant floating-point variables with the best possible precision so you don't need to type them everywhere,
and precomputing their inverses allows you to change division (which is slow) into multiplication e.g. x / PI = INV_PI * x
*/

//const float     PI = glm::uintBitsToFloat(0x40490FDBU); // Pi
//const float INV_PI = glm::uintBitsToFloat(0x3EA2F983U); // 1 / Pi
//const float     TAU = glm::uintBitsToFloat(0x40C90FDBU); // 2 * Pi = "Tau"
//const float INV_TAU = glm::uintBitsToFloat(0x3E22F983U); // 1 / (2 * Pi)

//const float     PI = 3.1415926535897932384626433832795028841971693993751058209749445923; // Pi
//const float INV_PI = 0.3183098861837906715377675267450287240689192914809128974953346881; // 1 / Pi

//const float     TAU = 6.2831853071795864769252867665590057683943387987502116419498891846; // 2 * Pi = "Tau"
//const float INV_TAU = 0.1591549430918953357688837633725143620344596457404564487476673441; // 1 / (2 * Pi)

// 2^32 (but slightly different, so 32-bit RNG doesn't ever reach 0.0 or 1.0)
// I know I said "C++11 lacks a real way to do this" earlier but there's no real way to get the closest
// but still *slightly* different float besides using hexadecimal
//const float     RNG32_MAX = uintBitsToFloat(0x4F7FFFFFU); // Actual: 0x4F800000U
//const float INV_RNG32_MAX = uintBitsToFloat(0x2F7FFFFFU); // Actual: 0x2F800000U
const float INV_RNG32_MAX = 2.32830629776081821092e-10f;

// Units
// https://www.seas.upenn.edu/~amyers/NaturalUnits.pdf

//const float constant_c    = 299792458.0f; // Speed of Light (m/s, Meters per Second)
//const float constant_c_km = 299792.4580f; // Speed of Light (km/s, Kilometers per Second)
//const float G = 6.6743015e-11; // Gravitational Constant
