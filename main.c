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
  size_t i, j, k;

  /* Load the data */
  complex x[NUM_TX_ANT][NUM_SC];              /* Transmitted signal */
  complex H[NUM_RX_ANT][NUM_TX_ANT][NUM_SC];  /* Channel */
  complex R[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];  /* Noise covariance matrix */
  complex y[NUM_RX_ANT][NUM_SC];              /* Received signal */

  data_t x_t[NUM_TX_ANT][NUM_SC];
  load_data("data/x_re.bin", NUM_TX_ANT * NUM_SC, x_t);
  for (i = 0; i < NUM_TX_ANT; i++)
    for (j = 0; j < NUM_SC; j++)
      x[i][j].re = x_t[i][j];
  load_data("data/x_im.bin", NUM_TX_ANT * NUM_SC, x_t);
  for (i = 0; i < NUM_TX_ANT; i++)
    for (j = 0; j < NUM_SC; j++)
      x[i][j].im = x_t[i][j];
  data_t H_t[NUM_RX_ANT][NUM_TX_ANT][NUM_SC];
  load_data("data/H_re.bin", NUM_RX_ANT * NUM_TX_ANT * NUM_SC, H_t);
  for (i = 0; i < NUM_RX_ANT; i++)
    for (j = 0; j < NUM_TX_ANT; j++)
      for (k = 0; k < NUM_SC; k++)
        H[i][j][k].re = H_t[i][j][k];
  load_data("data/H_im.bin", NUM_RX_ANT * NUM_TX_ANT * NUM_SC, H_t);
  for (i = 0; i < NUM_RX_ANT; i++)
    for (j = 0; j < NUM_TX_ANT; j++)
      for (k = 0; k < NUM_SC; k++)
        H[i][j][k].im = H_t[i][j][k];
  data_t R_t[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  load_data("data/R_re.bin", NUM_TX_ANT * NUM_TX_ANT * NUM_SC, R_t);
  for (i = 0; i < NUM_TX_ANT; i++)
    for (j = 0; j < NUM_TX_ANT; j++)
      for (k = 0; k < NUM_SC; k++)
        R[i][j][k].re = R_t[i][j][k];
  load_data("data/R_im.bin", NUM_TX_ANT * NUM_TX_ANT * NUM_SC, R_t);
  for (i = 0; i < NUM_TX_ANT; i++)
    for (j = 0; j < NUM_TX_ANT; j++)
      for (k = 0; k < NUM_SC; k++)
        R[i][j][k].im = R_t[i][j][k];
  data_t y_t[NUM_RX_ANT][NUM_SC];
  load_data("data/y_re.bin", NUM_RX_ANT * NUM_SC, y_t);
  for (i = 0; i < NUM_RX_ANT; i++)
    for (j = 0; j < NUM_SC; j++)
      y[i][j].re = y_t[i][j];
  load_data("data/y_im.bin", NUM_RX_ANT * NUM_SC, y_t);
  for (i = 0; i < NUM_RX_ANT; i++)
    for (j = 0; j < NUM_SC; j++)
      y[i][j].im = y_t[i][j];

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
