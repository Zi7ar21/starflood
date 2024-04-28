#include <cstdint>
#include <cstdio>
#include <cstdlib>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

void sphereIntersect(float ro[3], float rd[3], float ce[3], float ra) {
	float oc[3] = {ro[0] - ce[0], ro[1] - };
}

// sphere of size ra centered at point ce
vec2 sphIntersect(vec3 ro, vec3 rd, vec3 ce, float ra) {
	vec3 oc = ro - ce;
	float b = dot( oc, rd );
	vec3 qc = oc - b*rd;
	float h = ra*ra - dot( qc, qc );
	if( h<0.0 ) return vec2(-1.0); // no intersection
	h = sqrt( h );
	return vec2( -b-h, -b+h );
}

int main(int argc, char **argv) {
	if(argc > 1) {
		printf("%s: unrecognized option '%s'\n", argv[0], argv[1]);

		return EXIT_FAILURE;
	}

	size_t render_w = 384, render_h = 256;

	size_t render_buffer_size = sizeof(float) * 3 * render_w * render_h;

	float* render_buffer = (float*)malloc(render_buffer_size);

	if(render_buffer == NULL) {
		perror("Error");

		return EXIT_FAILURE;
	}

	for(size_t i = 0; i < (3 * render_w * render_h); i++) render_buffer[i] = 0;

	for(size_t y = 0; y < render_h; y++) {
		for(size_t x = 0; x < render_w; x++) {
			float uv[2] = {
				(((float)x+0.5f)-(0.5f*(float)render_w))/(0.5f*(float)render_h),
				(((float)y+0.5f)-(0.5f*(float)render_h))/(0.5f*(float)render_h)
			};

			float ro[3] = { 0.0f,  0.0f,  5.0f};
			float rd[3] = {uv[0], uv[1], -1.0f};

			float rd_r = 1.0f / sqrtf(rd[0]*rd[0]+rd[1]*rd[1]+1.0f);

			rd[0] *= rd_r;
			rd[1] *= rd_r;
			rd[2] *= rd_r;

			render_buffer[3*(render_w*y+x)+0] = rd[0];
			render_buffer[3*(render_w*y+x)+1] = rd[1];
			render_buffer[3*(render_w*y+x)+2] = rd[2];
		}
	}

	for(size_t i = 0; i < (3 * render_w * render_h); i++) render_buffer[i] = fmaxf(render_buffer[i], 0.0f);

	stbi_flip_vertically_on_write(1);

	const char* output_filename = "out/output.hdr";

	if(stbi_write_hdr(output_filename, (int)render_w, (int)render_h, 3, render_buffer) == 0) {
		printf("Error: stb_write_hdr(\"%s\", %d, %d, %d, render_buffer) returned 0!\n", output_filename, (int)render_w, (int)render_h, 3);

		free(render_buffer);

		return EXIT_FAILURE;
	}

	free(render_buffer);

	return EXIT_SUCCESS;
}
