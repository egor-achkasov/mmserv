#include "include/define.h"

#include "../common/printf.h"

/* extern functions */
extern void mmse();
extern acc_t mse();

/* Raw data */
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

/* Same data but casted to vcomplex */
vcomplex g_x = { .re = (data_t *)x_re, .im = (data_t *)x_im };
vcomplex g_H = { .re = (data_t *)H_re, .im = (data_t *)H_im };
vcomplex g_R = { .re = (data_t *)R_re, .im = (data_t *)R_im };
vcomplex g_y = { .re = (data_t *)y_re, .im = (data_t *)y_im };

/* MMSE approximation will be stored in this global*/
data_t x_MMSE_re[NUM_TX_ANT][NUM_SC];
data_t x_MMSE_im[NUM_TX_ANT][NUM_SC];
vcomplex g_x_MMSE = { .re = (data_t *)x_MMSE_re, .im = (data_t *)x_MMSE_im };

int main() {
  /* Calculate the MMSE approximation */
  mmse();
  printf("MSE: %f\n", mse());

  printf("Shutting down\n");
}
