#include "include/mmserv.h"

#include "../common/printf.h"

/* Transmitted signal */
extern data_t x_re[NUM_TX_ANT][NUM_SC];
extern data_t x_im[NUM_TX_ANT][NUM_SC];
/* Channel */
extern data_t H_re[NUM_RX_ANT][NUM_TX_ANT][NUM_SC];
extern data_t H_im[NUM_RX_ANT][NUM_TX_ANT][NUM_SC];
/* Noise covariance matrix */
extern data_t R_re[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
extern data_t R_im[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
/* Received signal */
extern data_t y_re[NUM_RX_ANT][NUM_SC];
extern data_t y_im[NUM_RX_ANT][NUM_SC];

int main() {
  /* Cast the data to vcomplex */
  vcomplex x, H, R, y;
  x.re = (data_t *)x_re;
  x.im = (data_t *)x_im;
  H.re = (data_t *)H_re;
  H.im = (data_t *)H_im;
  R.re = (data_t *)R_re;
  R.im = (data_t *)R_im;
  y.re = (data_t *)y_re;
  y.im = (data_t *)y_im;

  /* Calculate the MMSE approximation */
  data_t x_MMSE_re[NUM_TX_ANT][NUM_SC];
  data_t x_MMSE_im[NUM_TX_ANT][NUM_SC];
  vcomplex x_MMSE;
  x_MMSE.re = (data_t *)x_MMSE_re;
  x_MMSE.im = (data_t *)x_MMSE_im;
  mmse(&H, &y, &R, &x_MMSE);
  printf("%f\n", mse(&x, &x_MMSE));
}
