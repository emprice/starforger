starforger
==========

A transit light curve simulator and animator with no server-side dependencies.

# Doing math fast in the browser environment

Javascript was not natively designed to do the "heavy" math operations that
are necessary for computing the Mandel & Agol analytic quadratic transit
model. Javascript provides elementary math functions like `exp` and `pow`
through the `Math` namespace, but evaluating elliptic integrals is far
beyond the scope of typical Javascript.

One way to address this limitation is to use WebAssembly, an optimized,
compiled web binary format, to do the expensive calculations needed to
compute the transit flux. The GNU Scientific Library (GSL) provides elliptic
integral functions written in C, so we use
[Emscripten](https://emscripten.org/index.html) to port that code into
a WebAssembly (WASM) library. After installing and activating the Emscripten SDK
following [these instructions](https://emscripten.org/docs/getting_started/downloads.html#installation-instructions-using-the-emsdk-recommended),
the following shell commands, executed under the `src/gsl` subdirectory,
 will compile the full GSL into WASM.

```sh
./autogen.sh
emconfigure ./configure
emmake make LDFLAGS=-all-static
```

<!-- vim: set ft=markdown: -->
