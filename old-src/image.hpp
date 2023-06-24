#pragma once

#include <cstdint>
#include <cstdlib>
#include <string>

class RGB32F {
	public:
		float r, g, b;

		RGB32F(                            ): r( 0), g( 0), b( 0) {}
		RGB32F(float _r                    ): r(_r), g(_r), b(_r) {}
		RGB32F(float _r, float _g, float _b): r(_r), g(_g), b(_b) {}
};

class RGBA32F {
	public:
		float r, g, b, a;

		RGBA32F(                                      ): r( 0), g( 0), b( 0), a( 0) {}
		RGBA32F(float _r                              ): r(_r), g(_r), b(_r), a(_r) {}
		RGBA32F(float _r, float _g, float _b, float _a): r(_r), g(_g), b(_b), a(_a) {}
};

template<typename T> class Buffer2D {
	public:
		T* data;
		uint32_t size_x, size_y;

		Buffer2D(): size_x(0), size_y(0) {
			data = nullptr;
		}

		Buffer2D(uint32_t _size_x, uint32_t _size_y): size_x(_size_x), size_y(_size_y) {
			data = (T*)calloc(size_x * size_y, sizeof(T));
		}

		// destructor
		~Buffer2D() {
			free(data);
		}

		// https://en.cppreference.com/w/cpp/language/rule_of_three

		// copy constructor
		Buffer2D(const Buffer2D& other) {
			size_x = other.size_x;
			size_y = other.size_y;

			data = (T*)calloc(size_x * size_y, sizeof(T));

			for(uint32_t i = (uint32_t)0u; i < (size_x * size_y); i++) data[i] = other.data[i];
		}

		// copy assignment
		Buffer2D &operator=(const Buffer2D& other) {
			if(this == &other) return *this;

			size_x = other.size_x;
			size_y = other.size_y;

			free(data);

			data = (T*)calloc(size_x * size_y, sizeof(T));

			for(uint32_t i = (uint32_t)0u; i < (size_x * size_y); i++) data[i] = other.data[i];

			return *this;
		}
};

void WriteImageBuffer   (Buffer2D<RGBA32F>& image, std::string filename); // Saves a RGBA32F image as an 8-BPC PNG
void WriteImageBufferHDR(Buffer2D<RGBA32F>& image, std::string filename); // Saves a RGBA32F image as a 32-BPC Floating-Point HDR

//void ReadImage(Image2D<RGBA32F>& image, std::string filename);

//RGBA32F TexelFetch(Image2D<RGBA32F>& image, uint32_t coord_x, uint32_t coord_y);

//RGBA32F Texture2D(Image2D<RGBA32F>& image, float uv_x, float uv_y);
