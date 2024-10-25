#include "inc/mmserv.h"

int main() {
  /* Load the data */
  uint16_t x_raw[NUM_TX_ANT][NUM_SC][2];              /* Transmitted signal */
  uint16_t H_raw[NUM_RX_ANT][NUM_TX_ANT][NUM_SC][2];  /* Channel */
  uint16_t R_raw[NUM_TX_ANT][NUM_TX_ANT][NUM_SC][2];  /* Noise covariance matrix */
  uint16_t y_raw[NUM_RX_ANT][NUM_SC][2];              /* Received signal */
  load_data("data/x.bin", NUM_TX_ANT*NUM_SC*2, x_raw);
  load_data("data/H.bin", NUM_RX_ANT*NUM_TX_ANT*NUM_SC*2, H_raw);
  load_data("data/y.bin", NUM_RX_ANT*NUM_SC*2, y_raw);
  load_data("data/R.bin", NUM_RX_ANT*NUM_SC*2, R_raw);

  /* Cast the data into complex data structures */
  complex x[NUM_TX_ANT][NUM_SC];
  complex H[NUM_RX_ANT][NUM_TX_ANT][NUM_SC];
  complex y[NUM_RX_ANT][NUM_SC];
  complex R[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  for (uint32_t i = 0; i < NUM_TX_ANT; ++i)
    for (uint32_t j = 0; j < NUM_SC; ++j){
      x[i][j].re = x_raw[i][j][0];
      x[i][j].im = x_raw[i][j][1];
    }
  for (uint32_t i = 0; i < NUM_RX_ANT; ++i)
    for (uint32_t j = 0; j < NUM_TX_ANT; ++j)
      for (uint32_t k = 0; k < NUM_SC; ++k){
        H[i][j][k].re = H_raw[i][j][k][0];
        H[i][j][k].im = H_raw[i][j][k][1];
      }
  for (uint32_t i = 0; i < NUM_RX_ANT; ++i)
    for (uint32_t j = 0; j < NUM_SC; ++j){
      y[i][j].re = y_raw[i][j][0];
      y[i][j].im = y_raw[i][j][1];
    }
  for (uint32_t i = 0; i < NUM_TX_ANT; ++i)
    for (uint32_t j = 0; j < NUM_TX_ANT; ++j)
      for (uint32_t k = 0; k < NUM_SC; ++k){
        R[i][j][k].re = R_raw[i][j][k][0];
        R[i][j][k].im = R_raw[i][j][k][1];
      }

  /* Calculate the MMSE approximation */
  complex x_MMSE[NUM_TX_ANT][NUM_SC];
  mmse_nosqrt(H, y, R, x_MMSE);

  /* Interleave the result */
  uint16_t res[NUM_TX_ANT][NUM_SC][2];
  for (uint32_t i = 0; i < NUM_TX_ANT; ++i)
    for (uint32_t j = 0; j < NUM_SC; ++j){
      res[i][j][0] = x_MMSE[i][j].re;
      res[i][j][1] = x_MMSE[i][j].im;
    }

  /* Save the result */
  save_data("out/x_mmse.bin", res);
}
