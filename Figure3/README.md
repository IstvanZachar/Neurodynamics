Generating data and image for Figure 3.

1. Compile `storkey_retrain.c` with accompanied `makefile`.
2. Run the resulting  binary`aann` to generate data to the standard output.
   Successive runs should be generated by modifying the macro for the number of retrains `#define RETRAINNUM` in `storkey_retrain.c`.
3. Generated data should be saved in `.dat` files, e.g. `retrain.dat`. Pre-generated data files are provided as `retrain*.dat`.
5. Run the *Mathematica* code `abstract.nb` (optimized for *Mathematica* 11.0.0.0),
   it generates simulated abstract data and creates a figure with both the C-code generated and *Mathematica* generated data. Pre-generated image file `abstract.tiff` is provided.

By István Zachar
2016
