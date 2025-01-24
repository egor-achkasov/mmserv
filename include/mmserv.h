#ifndef __MMSERV_H
#define __MMSERV_H

#include "define.h"

#include <stddef.h> /* for size_t */

/*
 * Typedefs
 */

typedef struct {
  data_t* re;
  data_t* im;
} vcomplex;

/*
 * MMSE
 */

/** Calculate MMSE estimation of x in y = H*x + n
 * \param H matrix of channel coefficients. Shape [NUM_RX_ANT][NUM_TX_ANT][NUM_SC]
 * \param y received signal. Shape [NUM_RX_ANT][NUM_SC]
 * \param R noise covariance matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param x_MMSE output MMSE estimation of x. Shape [NUM_TX_ANT][NUM_SC]
 */
void mmse(
  IN vcomplex *H,
  IN vcomplex *y,
  IN vcomplex *R,
  OUT vcomplex *x_MMSE);

/** Calculate MMSE estimation of x in y = H*x + n. Uses a sqrt-free Cholesky decomposition.
 * \param H matrix of channel coefficients. Shape [NUM_RX_ANT][NUM_TX_ANT][NUM_SC]
 * \param y received signal. Shape [NUM_RX_ANT][NUM_SC]
 * \param R noise covariance matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param x_MMSE output MMSE estimation of x. Shape [NUM_TX_ANT][NUM_SC]
 */
void mmse_nosqrt(
  IN vcomplex *H,
  IN vcomplex *y,
  IN vcomplex *R,
  OUT vcomplex *x_MMSE);

/*
 * MSE
 */

/** Calculate mean squared error between the original x and the estimated x
 * \param x original x. Shape [NUM_TX_ANT][NUM_SC]
 * \param x_MMSE estimated x. Shape [NUM_TX_ANT][NUM_SC]
 * \return mean squared error
 */
acc_t mse(
  IN vcomplex *x,
  IN vcomplex *x_MMSE);

#endif /* __MMSERV_H */
