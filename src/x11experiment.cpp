/*
	Simple Xlib application for creating a window and drawing a box in it.
	gcc input.c -o output -lX11
*/

#include <X11/Xlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
	XEvent event;
	char str[256];
	snprintf(str, 256, "Hello, world!");
	int screen;
 
	// open connection to the server
	Display* display = XOpenDisplay(NULL);

	if(display == NULL) {
		printf("Error: Cannot XOpenDisplay() returned NULL, cannot open display!\n");

		return EXIT_FAILURE;
	}

	screen = DefaultScreen(display);

	// create window
	Window window = XCreateSimpleWindow(display, RootWindow(display, screen), 10, 10, 200, 200, 1, BlackPixel(display, screen), WhitePixel(display, screen));

	// select kind of events we are interested in
	XSelectInput(display, window, ExposureMask | ButtonPressMask | KeyPressMask);

	// map (show) the window
	XMapWindow(display, window);

	// event loop
	while(true) {
		XNextEvent(display, &event);
 
		// draw or redraw the window
		if(event.type == Expose) {
			XFillRectangle(display, window, DefaultGC(display, screen), 20, 20, 10, 10);
			XDrawString(display, window, DefaultGC(display, screen), 50, 50, str, strlen(msg));
		}

		if(event.type == ButtonPress) {
			XDrawLine(display, window, DefaultGC(display, screen), 0, 0, 100, 100);
		}

		// exit on key press
		if(event.type == KeyPress) break;
	}

	// close connection to the server
	XCloseDisplay(display);

	return EXIT_SUCCESS;
}