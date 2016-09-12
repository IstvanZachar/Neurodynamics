# Neurodynamics
Code to generate data for the manuscript "Breeding novel solutions in the brain: A model of Darwinian neurodynamics"

Files in this repository:
  - meta7.c   Program code in C to simulate selection and evolution with periodically changing environment.
     - randomGenerator.c   Function library to be included with `meta7.c` containing custom pseudorandom number generators.
     - colors.c            Function library to be included with `meta7.c` containing some higlighting functionality.
  - meta.nb   Program code in Wolfram Language (to be run with WRI's Mathematica) to simulate different attractor networks (abstract metamodel).
