#include "include/mmserv.h"

#include "../common/util.h"
#include "printf.h"

/* Transmitted signal */
extern complex x[NUM_TX_ANT][NUM_SC];
/* Channel */
extern complex H[NUM_RX_ANT][NUM_TX_ANT][NUM_SC];
/* Noise covariance matrix */
extern complex R[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
/* Received signal */
extern complex y[NUM_RX_ANT][NUM_SC];

int main() {
  complex x_MMSE[NUM_TX_ANT][NUM_SC];
  mmse(H, y, R, x_MMSE);
  printf("MSE: %f\n", mse(x, x_MMSE));
}
