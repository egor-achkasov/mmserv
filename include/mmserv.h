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
 * Complex matrix operations
 */

/** Calculate the Hermitian transpose of a matrix A
 * \param A input matrix. Shape [NUM_RX_ANT][NUM_TX_ANT][NUM_SC]
 * \param AH output matrix. Shape [NUM_TX_ANT][NUM_RX_ANT][NUM_SC]
 */
void cmat_hermitian_transpose_RxTx(
  IN vcomplex *A,
  OUT vcomplex *AH);

/** Calculate the Hermitian transpose of a matrix A
 * \param A input matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param AH output matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 */
void cmat_hermitian_transpose_TxTx(
  IN vcomplex *A,
  OUT vcomplex *AH);

/** Multiply two matrices A and B
 * \param A input matrix. Shape [NUM_TX_ANT][NUM_RX_ANT][NUM_SC]
 * \param B input matrix. Shape [NUM_RX_ANT][NUM_TX_ANT][NUM_SC]
 * \param result output matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 */
void cmatmul_TxRx_RxTx(
  IN vcomplex *A,
  IN vcomplex *B,
  OUT vcomplex *result);

/** Multiply two matrices A and B
 * \param A input matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param B input matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param result output matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 */
void cmatmul_TxTx_TxTx(
  IN vcomplex *A,
  IN vcomplex *B,
  OUT vcomplex *result);

/** Add two matrices A and B
 * \param A input matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param B input matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param result output matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 */
void cmatadd_TxTx(
  IN vcomplex *A,
  IN vcomplex *B,
  OUT vcomplex *result);

/** Perfom Cholesky decomposition of a square matrix A
 * such that A = L*L^H with Cholesky–Banachiewicz algorithm
 * \param A input square matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param L output lower triangular matrix L. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 */
void ccholesky_TxTx(
  IN vcomplex *A,
  OUT vcomplex *L);

/** Perfom Cholesky decomposition of a square matrix A
 * such that A = L*D*L^H with Cholesky–Banachiewicz algorithm,
 * where D is a diagonal matrix and
 * L is a lower triangular matrix with ones on the diagonal.
 * This function does not calculate use csqrt and sqrt functions.
 * \param A input square matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param L output lower triangular matrix L. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param D output diagonal of the diagonal matrix D. Shape [NUM_TX_ANT][NUM_SC]
 */
void ccholesky_nosqrt_TxTx(
  IN vcomplex *A,
  OUT vcomplex *L,
  OUT vcomplex *D);

/** Multiply a matrix A with a vector b
 * \param A input matrix. Shape [NUM_TX_ANT][NUM_RX_ANT][NUM_SC]
 * \param b input vector. Shape [NUM_RX_ANT][NUM_SC]
 * \param result output vector. Shape [NUM_TX_ANT][NUM_SC]
 */
void cmatvecmul_TxRx(
  IN vcomplex *A,
  IN vcomplex *b,
  OUT vcomplex *result);

/** Find z from L*z = b with a forward substitution
 * \param L lower triangular matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param b rhs vector. Shape [NUM_TX_ANTS][NUM_SC]
 * \param result output vector z. Shape [NUM_TX_ANT][NUM_SC]
 */
void cforwardsub_TxTx(
  IN vcomplex *L,
  IN vcomplex *b,
  OUT vcomplex *result);

/** Find x from U*x = b with a backward substitution
 * \param U upper triangular matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param b rhs vector. Shape [NUM_TX_ANTS][NUM_SC]
 * \param result output vector x. Shape [NUM_TX_ANT][NUM_SC]
 */
void cbackwardsub_TxTx(
  IN vcomplex *U,
  IN vcomplex *b,
  OUT vcomplex *result);

/** Calculate mean squared error between the original x and the estimated x
 * \param x original x. Shape [NUM_TX_ANT][NUM_SC]
 * \param x_MMSE estimated x. Shape [NUM_TX_ANT][NUM_SC]
 * \return mean squared error
 */
acc_t mse(
  IN vcomplex *x,
  IN vcomplex *x_MMSE);

#endif /* __MMSERV_H */
