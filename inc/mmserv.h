#ifndef __MMSERV_H
#define __MMSERV_H

#include "define.h"

#include <stdint.h> // for int16_t, int32_t
#include <stddef.h> // for size_t
#include <stdio.h>  // for FILE

/*
 * Typedefs
 */

typedef struct {
    data_t re;
    data_t im;
} complex;

/** Calculate MMSE estimation of x in y = H*x + n
 * \param H matrix of channel coefficients. Shape [NUM_RX_ANT][NUM_TX_ANT][NUM_SC]
 * \param y received signal. Shape [NUM_RX_ANT][NUM_SC]
 * \param R noise covariance matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param x_MMSE output MMSE estimation of x. Shape [NUM_TX_ANT][NUM_SC]
 */
void mmse(
  IN complex H[NUM_RX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex y[NUM_RX_ANT][NUM_SC],
  IN complex R[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex x_MMSE[NUM_TX_ANT][NUM_SC]);

/** Calculate MMSE estimation of x in y = H*x + n. Uses a sqrt-free Cholesky decomposition.
 * \param H matrix of channel coefficients. Shape [NUM_RX_ANT][NUM_TX_ANT][NUM_SC]
 * \param y received signal. Shape [NUM_RX_ANT][NUM_SC]
 * \param R noise covariance matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param x_MMSE output MMSE estimation of x. Shape [NUM_TX_ANT][NUM_SC]
 */
void mmse_nosqrt(
  IN complex H[NUM_RX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex y[NUM_RX_ANT][NUM_SC],
  IN complex R[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex x_MMSE[NUM_TX_ANT][NUM_SC]);

/*
 * Complex operations
 */

/** Create a complex number from real and imaginary parts */
complex cmake(IN data_t re, IN data_t im);

/** Complex multiplication */
complex cmul(IN complex a, IN complex b);

/** Complex absolute value squared */
data_t cabs2(IN complex a);

/** Complex addition */
complex cadd(IN complex a, IN complex b);
void cadd_acc(IN complex a, IN complex b);

/** Complex subtraction */
complex csub(IN complex a, IN complex b);

/** Complex conjugate */
complex cconj(IN complex a);

/** Complex division */
complex cdiv(IN complex a, IN complex b);

/** Square root of a natural number */
data_t sqrt(IN data_t x);

/** Complex root */
complex csqrt(IN complex a);

/*
 * Complex matrix operations
 */

/** Calculate the Hermitian transpose of a matrix A
 * \param A input matrix. Shape [NUM_RX_ANT][NUM_TX_ANT][NUM_SC]
 * \param AH output matrix. Shape [NUM_TX_ANT][NUM_RX_ANT][NUM_SC]
 */
void cmat_hermitian_transpose_RxTx(
  IN complex A[NUM_RX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex AH[NUM_TX_ANT][NUM_RX_ANT][NUM_SC]);

/** Calculate the Hermitian transpose of a matrix A
 * \param A input matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param AH output matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 */
void cmat_hermitian_transpose_TxTx(
  IN complex A[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex AH[NUM_TX_ANT][NUM_TX_ANT][NUM_SC]);

/** Multiply two matrices A and B
 * \param A input matrix. Shape [NUM_TX_ANT][NUM_RX_ANT][NUM_SC]
 * \param B input matrix. Shape [NUM_RX_ANT][NUM_TX_ANT][NUM_SC]
 * \param result output matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 */
void cmatmul_TxRx_RxTx(
  IN complex A[NUM_TX_ANT][NUM_RX_ANT][NUM_SC],
  IN complex B[NUM_RX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_TX_ANT][NUM_SC]);

/** Multiply two matrices A and B
 * \param A input matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param B input matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param result output matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 */
void cmatmul_TxTx_TxTx(
  IN complex A[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex B[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_TX_ANT][NUM_SC]);

/** Add two matrices A and B
 * \param A input matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param B input matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param result output matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 */
void cmatadd_TxTx(
  IN complex A[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex B[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_TX_ANT][NUM_SC]);

/** Perfom Cholesky decomposition of a square matrix A
 * such that A = L*L^H with Cholesky–Banachiewicz algorithm
 * \param A input square matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param L output lower triangular matrix L. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 */
void ccholesky_TxTx(
  IN complex A[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex L[NUM_TX_ANT][NUM_TX_ANT][NUM_SC]);

/** Perfom Cholesky decomposition of a square matrix A
 * such that A = L*D*L^H with Cholesky–Banachiewicz algorithm,
 * where D is a diagonal matrix and
 * L is a lower triangular matrix with ones on the diagonal.
 * This function does not calculate use csqrt and sqrt functions.
 * \param A input square matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param L output lower triangular matrix L. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param D output diagonal matrix L. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 */
void ccholesky_nosqrt_TxTx(
  IN complex A[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex L[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex D[NUM_TX_ANT][NUM_TX_ANT][NUM_SC]);

/** Multiply a matrix A with a vector b
 * \param A input matrix. Shape [NUM_TX_ANT][NUM_RX_ANT][NUM_SC]
 * \param b input vector. Shape [NUM_RX_ANT][NUM_SC]
 * \param result output vector. Shape [NUM_TX_ANT][NUM_SC]
 */
void cmatvecmul_TxRx(
  IN complex A[NUM_TX_ANT][NUM_RX_ANT][NUM_SC],
  IN complex b[NUM_RX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_SC]);

/** Find z from L*z = b with a forward substitution
 * \param L lower triangular matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param b rhs vector. Shape [NUM_TX_ANTS][NUM_SC]
 * \param result output vector z. Shape [NUM_TX_ANT][NUM_SC]
 */
void cforwardsub_TxTx(
  IN complex L[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex b[NUM_TX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_SC]);

/** Find x from U*x = b with a backward substitution
 * \param U upper triangular matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \param b rhs vector. Shape [NUM_TX_ANTS][NUM_SC]
 * \param result output vector x. Shape [NUM_TX_ANT][NUM_SC]
 */
void cbackwardsub_TxTx(
  IN complex U[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex b[NUM_TX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_SC]);

/*
 * IO
 */

/** Load data from a binary file
 * \param file path to the file
 * \param size number of elements to read
 * \param out output array
 */
void load_data(
  IN const char *file,
  IN size_t size,
  OUT data_t *out);

/** Save data to a binary file
 * \param file path to the file
 * \param x_mmse data to save
 */
void save_data(
  IN const char* file,
  IN data_t x_mmse[NUM_TX_ANT][NUM_SC][2]);


#endif /* __MMSERV_H */
