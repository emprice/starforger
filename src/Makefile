CFLAGS = -Wall -O3 -I gsl -ferror-limit=5 -std=gnu99
LDFLAGS = -s WASM=1 -s EXPORTED_FUNCTIONS="['_malloc', '_free']" \
		  -s ALLOW_MEMORY_GROWTH=1 -s MODULARIZE=1 -s EXPORT_NAME=compute \
		  -L gsl/.libs -L gsl/cblas/.libs -lgsl -lgslcblas -lm

default: ../js/compute.js

../js/compute.js: blackbody.o mandelagol.o orbit.o extras.o
	emcc $(CFLAGS) $(LDFLAGS) $^ -o $@

%.o: %.c Makefile
	emcc $(CFLAGS) -c $< -o $@

# vim: set ft=make:
