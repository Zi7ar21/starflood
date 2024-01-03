#include "image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

void write_image_hdr(char* filename, int size_x, int size_y, float* data) {
	stbi_flip_vertically_on_write(1);

	stbi_write_hdr(filename, size_x, size_y, 3, data);
}
