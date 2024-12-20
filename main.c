#include "include/mmserv.h"

#include <stdio.h>  /* for printf */
#include <stdlib.h> /* for exit */

void load_data(
  IN char* filepath,
  IN size_t size,
  OUT void* buff)
{
  FILE* f = fopen(filepath, "r");
  if (!f) {
    fprintf(stderr, "Error: could not open file %s\n", filepath);
    exit(8);
  }
  if ((fread(buff, sizeof(data_t), size, f)) != size) {
    fprintf(stderr, "Error: could not read file %s\n", filepath);
    exit(8);
  }
  fclose(f);
}

int main() {
  /* Load the data */
  complex x[NUM_TX_ANT][NUM_SC];              /* Transmitted signal */
  complex H[NUM_RX_ANT][NUM_TX_ANT][NUM_SC];  /* Channel */
  complex R[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];  /* Noise covariance matrix */
  complex y[NUM_RX_ANT][NUM_SC];              /* Received signal */

  load_data("data/x.bin", NUM_TX_ANT * NUM_SC * 2, x);
  load_data("data/H.bin", NUM_RX_ANT * NUM_TX_ANT * NUM_SC * 2, H);
  load_data("data/R.bin", NUM_TX_ANT * NUM_TX_ANT * NUM_SC * 2, R);
  load_data("data/y.bin", NUM_RX_ANT * NUM_SC * 2, y);

  /* Calculate the MMSE approximation */
  complex x_MMSE[NUM_TX_ANT][NUM_SC];
  mmse(H, y, R, x_MMSE);
  printf("%f\n", mse(x, x_MMSE));

  /* Save the result */
  FILE* f = fopen("out/x_mmse.bin", "w");
  if (!f) {
    fprintf(stderr, "Error: could not open file out/x_MMSE.bin\n");
    exit(8);
  }
  if ((fwrite(x_MMSE, sizeof(data_t), NUM_TX_ANT * NUM_SC * 2, f)) != NUM_TX_ANT * NUM_SC * 2) {
    fprintf(stderr, "Error: could not write file out/x_MMSE.bin\n");
    exit(8);
  }
  fclose(f);
}
