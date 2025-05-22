#ifndef COMMON_H
#define COMMON_H

#include <stdint.h> /* for uint64_t, uint32_t */


/*
 * Typedefs
 */

#if defined(DATA_TYPE_float)
typedef float data_t;
typedef float acc_t;
#elif defined(DATA_TYPE_fixed)
typedef int32_t data_t;
typedef int64_t acc_t;
#define FP_Q 31
#else
#error "Please define DATA_TYPE_float or DATA_TYPE_fixed"
#endif

typedef struct {
  data_t *re;
  data_t *im;
} vcomplex;


/*
 * Global variables
 */

/* Raw data */
/* Transmitted signal */
extern data_t x_re[NUM_TX * NUM_SC];
extern data_t x_im[NUM_TX * NUM_SC];
/* Channel */
extern data_t H_re[NUM_RX * NUM_TX * NUM_SC];
extern data_t H_im[NUM_RX * NUM_TX * NUM_SC];
/* Noise covariance matrix */
extern data_t R_re[NUM_TX * NUM_TX * NUM_SC];
extern data_t R_im[NUM_TX * NUM_TX * NUM_SC];
/* Received signal */
extern data_t y_re[NUM_RX * NUM_SC];
extern data_t y_im[NUM_RX * NUM_SC];
/* MMSE raw data */
extern data_t G_re[NUM_TX * NUM_TX * NUM_SC];
extern data_t G_im[NUM_TX * NUM_TX * NUM_SC];
extern data_t L_re[NUM_TX * NUM_TX * NUM_SC];
extern data_t L_im[NUM_TX * NUM_TX * NUM_SC];
extern data_t g_D[NUM_TX * NUM_SC]; /* no imaginary part in D */
extern data_t HHy_re[NUM_TX * NUM_SC];
extern data_t HHy_im[NUM_TX * NUM_SC];
extern data_t z_re[NUM_TX * NUM_SC];
extern data_t z_im[NUM_TX * NUM_SC];
/* Result of MMSE approximation */
extern data_t x_MMSE_re[NUM_TX * NUM_SC];
extern data_t x_MMSE_im[NUM_TX * NUM_SC];

/* Same data but casted to vcomplex */
extern vcomplex g_x;
extern vcomplex g_H;
extern vcomplex g_R;
extern vcomplex g_y;
extern vcomplex g_G;
extern vcomplex g_L;
extern vcomplex g_HHy;
extern vcomplex g_z;
extern vcomplex g_x_MMSE;


/*
 * Complex matrix operations
 */

extern void cmatgram();
extern void ccholesky();
extern void cmatvecmul();
extern void cforwardsub();
extern void cbackwardsub();

#endif
