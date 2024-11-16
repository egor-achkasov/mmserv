#include "include/mmserv.h"

#include "printf.h"

extern data_t x_raw[NUM_TX_ANT][NUM_SC][2];              /* Transmitted signal */
extern data_t H_raw[NUM_RX_ANT][NUM_TX_ANT][NUM_SC][2];  /* Channel */
extern data_t R_raw[NUM_TX_ANT][NUM_TX_ANT][NUM_SC][2];  /* Noise covariance matrix */
extern data_t y_raw[NUM_RX_ANT][NUM_SC][2];              /* Received signal */

int main() {
  uint32_t i, j, k;

  /* Cast the data into complex data structures */
  complex x[NUM_TX_ANT][NUM_SC];
  complex H[NUM_RX_ANT][NUM_TX_ANT][NUM_SC];
  complex y[NUM_RX_ANT][NUM_SC];
  complex R[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < NUM_SC; ++j){
      x[i][j].re = x_raw[i][j][0];
      x[i][j].im = x_raw[i][j][1];
    }
  for (i = 0; i < NUM_RX_ANT; ++i)
    for (j = 0; j < NUM_TX_ANT; ++j)
      for (k = 0; k < NUM_SC; ++k){
        H[i][j][k].re = H_raw[i][j][k][0];
        H[i][j][k].im = H_raw[i][j][k][1];
      }
  for (i = 0; i < NUM_RX_ANT; ++i)
    for (j = 0; j < NUM_SC; ++j){
      y[i][j].re = y_raw[i][j][0];
      y[i][j].im = y_raw[i][j][1];
    }
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < NUM_TX_ANT; ++j)
      for (k = 0; k < NUM_SC; ++k){
        R[i][j][k].re = R_raw[i][j][k][0];
        R[i][j][k].im = R_raw[i][j][k][1];
      }

  /* Calculate the MMSE approximation */
  complex x_MMSE[NUM_TX_ANT][NUM_SC];
  mmse(H, y, R, x_MMSE);

  /* Print MSE */
  printf("%f\n", mse(x, x_MMSE));
}
