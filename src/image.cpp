#include "image.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

//#define VERBOSE_LOGGING

/*
#include <glm/common.hpp>
using namespace glm;
*/

/*
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
*/

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

float clamp(float a, float b, float c) {
	return fminf(fmaxf(a, b), c);
}

void WriteImageBuffer(Buffer2D<RGBA32F>& image, std::string filename) {
	#ifdef VERBOSE_LOGGING
	std::cout
	<< "Writing image to disk:\n"
	<< "  File: " << filename << "\n"
	<< "  Size: " << image.size_x << "x" << image.size_y << "\n"
	<< "  Form: " << "RGB8\n"
	<< "  Type: " << "Portable Network Graphics (.png)\n"
	<< std::endl;
	#endif

	uint8_t* raw_image_buffer = new uint8_t[3u * image.size_x * image.size_y];

	stbi_flip_vertically_on_write(1);

	for(uint32_t i = (uint32_t)0u; i < (image.size_x * image.size_y); i++) {
		raw_image_buffer[3u * i + 0u] = (uint8_t)(255.0f * clamp(image.data[i].r, 0.0f, 1.0f));
		raw_image_buffer[3u * i + 1u] = (uint8_t)(255.0f * clamp(image.data[i].g, 0.0f, 1.0f));
		raw_image_buffer[3u * i + 2u] = (uint8_t)(255.0f * clamp(image.data[i].b, 0.0f, 1.0f));
	}

	int write_success = stbi_write_png(filename.c_str(), image.size_x, image.size_y, 3u, raw_image_buffer, 3u * image.size_x);

	delete[] raw_image_buffer;

	if(write_success != 0) {
		#ifdef VERBOSE_LOGGING
		std::cout << " " << filename << " written successfully.\n" << std::endl;
		#endif
	} else {
		std::cerr << " " << filename << " failed to write!\n" << std::endl;
	}
}

void WriteImageBufferHDR(Buffer2D<RGBA32F>& image, std::string filename) {
	#ifdef VERBOSE_LOGGING
	std::cout
	<< "Writing image to disk:\n"
	<< "  File: " << filename << "\n"
	<< "  Size: " << image.size_x << "x" << image.size_y << "\n"
	<< "  Form: " << "RGB32F\n"
	<< "  Type: " << "Radiance HDR (.hdr)\n"
	<< std::endl;
	#endif

	float* raw_image_buffer = new float[3u * image.size_x * image.size_y];

	stbi_flip_vertically_on_write(1);

	for(uint32_t i = (uint32_t)0u; i < (image.size_x * image.size_y); i++) {
		raw_image_buffer[3u * i + 0u] = image.data[i].r;
		raw_image_buffer[3u * i + 1u] = image.data[i].g;
		raw_image_buffer[3u * i + 2u] = image.data[i].b;
	}

	int write_success = stbi_write_hdr(filename.c_str(), image.size_x, image.size_y, 3u, raw_image_buffer);

	delete[] raw_image_buffer;

	if(write_success != 0) {
		#ifdef VERBOSE_LOGGING
		std::cout << " " << filename << " written successfully.\n" << std::endl;
		#endif
	} else {
		std::cout << " " << filename << " failed to write!\n" << std::endl;
	}
}

/*
void read_image(image_RGB32F &image, std::string filename) {
	std::cout
	<< "Loading image from disk:\n"
	<< "  File: " << filename << "\n"
	<< std::endl;

	int width, height, channels;

	//stbi_set_flip_vertically_on_load(1);

	float* data = stbi_loadf(filename.c_str(), &width, &height, &channels, 3);

	if(data == nullptr) {
		std::cerr << "Error: Failed to load " << filename << std::endl;

		return;
	}

	image.size = uvec2(width, height);

	#ifdef VERBOSE_LOGGING
	std::cout << "Creating a "
	<< image.size.x << "x" << image.size.y
	<< " (" << image.size.x * image.size.y * sizeof(vec3) << " byte)"
	<< " image buffer...\n" << std::endl;
	#endif

	image.data = new vec3[image.size.x * image.size.y];

	if(image.data == nullptr) {
		std::cerr << "Error: Failed to allocate an image buffer!" << std::endl;

		return;
	}

	uvec2 pixel_coord;

	//for(pixel_coord.y = 0u; pixel_coord.y < image.size.y; pixel_coord.y++) {
	//for(pixel_coord.x = 0u; pixel_coord.x < image.size.x; pixel_coord.x++) {
	//	image.data[pixel_coord.x + image.size.x * pixel_coord.y] = glm::vec3(
	//	data[3u * (pixel_coord.x + pixel_coord.y * image.size.x) + 0u],
	//	data[3u * (pixel_coord.x + pixel_coord.y * image.size.x) + 1u],
	//	data[3u * (pixel_coord.x + pixel_coord.y * image.size.x) + 2u]
	//	);
	//}
	//}

	for(unsigned int i = 0u; i < image.size.x * image.size.y; i++) image.data[i] = glm::vec3(data[3u * i + 0u], data[3u * i + 1u], data[3u * i + 2u]);

	#ifdef VERBOSE_LOGGING
	std::cout << "Closing " << filename << "..." << std::endl;
	#endif

	stbi_image_free(data);

	std::cout << "Loaded " << filename << " successfully.\n" << std::endl;
}

vec3 texel_fetch(image_RGB32F& image, uvec2 coord) {
	if(image.size.x < coord.x || image.size.y < coord.y) return vec3(0);

	return image.data[coord.x + image.size.x * coord.y];
}

vec3 texture2D_RGB32F(image_RGB32F &image, vec2 uv) {
	uv *= vec2(image.size.x, image.size.y);

	 vec2 f = glm::fract(uv);
	uvec2 i = glm::uvec2(uv);

	//vec2 k = smoothstep(0.0f, 1.0f, f);
	vec2 k = f;

	//uvec2 coord_0 = mod(i + uvec2(0u, 0u), image.size);
	//uvec2 coord_1 = mod(i + uvec2(1u, 0u), image.size);
	//uvec2 coord_2 = mod(i + uvec2(0u, 1u), image.size);
	//uvec2 coord_3 = mod(i + uvec2(1u, 1u), image.size);

	// For some reason, the above code does not work even if glm/gtc/integer.hpp is included, so we do it the hard way:
	uvec2 coord_0 = uvec2((i.x + 0u) % image.size.x, (i.y + 0u) % image.size.y);
	uvec2 coord_1 = uvec2((i.x + 1u) % image.size.x, (i.y + 0u) % image.size.y);
	uvec2 coord_2 = uvec2((i.x + 0u) % image.size.x, (i.y + 1u) % image.size.y);
	uvec2 coord_3 = uvec2((i.x + 1u) % image.size.x, (i.y + 1u) % image.size.y);

	vec3 t_0 = texel_fetch(image, coord_0);
	vec3 t_1 = texel_fetch(image, coord_1);
	vec3 t_2 = texel_fetch(image, coord_2);
	vec3 t_3 = texel_fetch(image, coord_3);

	return glm::mix(glm::mix(t_0, t_1, k.x), glm::mix(t_2, t_3, k.x), k.y);
}
*/
