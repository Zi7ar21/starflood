#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

int main(int argc, char** argv) {
	if(argc > 1) {
		if(!strcmp(argv[1],"-h") || !strcmp(argv[1],"--help")) {
			printf(
			"Usage: %s [OPTION]...\n"
			"  -h, --help       display this help and exit\n"
			//"      --version    output version information and exit\n"
			, argv[0]);

			return EXIT_SUCCESS;
		} else {
			printf("%s: unrecognized option '%s'\nTry '%s -h' or '%s --help' for more information.\n", argv[0], argv[1], argv[0], argv[0]);

			return EXIT_FAILURE;
		}
	}

	return EXIT_SUCCESS;
}