# Neurodynamics
Code to generate data for the manuscript "Breeding novel solutions in the brain: a model of Darwinian neurodynamics".
Other results in the paper were done by András Szilágy (found at: [TBA]).

Files in this repository:
  - `meta7.c`: program code in C to simulate selection and evolution within a periodically changing environment. Can also be used to test a ton of different things with attractor neural networks.
     - `randomGenerator.c`: function library to be included with `meta7.c` containing custom pseudorandom number generators.
     - `colors.c`: function library to be included with `meta7.c` containing some higlighting functionality.
  - `meta.nb`: program code in Wolfram Language (to be run with WRI's **Mathematica**) to simulate different attractor networks (abstract metamodel).
