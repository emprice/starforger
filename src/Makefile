default: ../js/mandelagol.js

../js/mandelagol.js: mandelagol.c Makefile
	emcc -O3 -s WASM=1 -s EXPORTED_FUNCTIONS="['_malloc', '_free']" -s ALLOW_MEMORY_GROWTH=1 -s MODULARIZE=1 -s EXPORT_NAME=mandelagol -I gsl -L gsl/.libs -lgsl -lm $< -o $@

# vim: set ft=make:
