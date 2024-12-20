#include "../include/mmserv.h"

#include "../../common/rivec/vector_defines.h"
#include "riscv_vector.h"

/* Profiling */
#ifdef DEBUG
#include <stdio.h>

#ifdef _WIN32
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

static uint64_t start, stop;
void start_timer()
{
  start = __rdtsc();
}
void stop_timer()
{
  stop = __rdtsc();
}
uint64_t get_timer()
{
  return stop - start;
}

#define TIME(msg, func, ...) \
  start_timer(); \
  func(__VA_ARGS__); \
  stop_timer(); \
  printf(msg, get_timer());

#else
#define TIME(msg, func, ...) func(__VA_ARGS__);

#endif

/*
 * Complex functions
 */

complex cmake(IN data_t re, IN data_t im)
{
  complex t;
  t.re = re;
  t.im = im;
  return t;
}
complex cmul(IN complex a, IN complex b)
{
  complex t;
  t.re = a.re * b.re - a.im * b.im;
  t.im = a.im * b.re + a.re * b.im;
  return t;
}
data_t cabs2(IN complex a)
{
  return a.re * a.re + a.im * a.im;
}
complex cadd(IN complex a, IN complex b)
{
  complex t;
  t.re = a.re + b.re;
  t.im = a.im + b.im;
  return t;
}
void cadd_acc(IN complex a, IN complex b)
{
  a.re += b.re;
  a.im += b.im;
}
complex csub(IN complex a, IN complex b)
{
  complex t;
  t.re = a.re - b.re;
  t.im = a.im - b.im;
  return t;
}
complex cconj(IN complex a)
{
  complex t;
  t.re = a.re;
  t.im = -a.im;
  return t;
}
complex cdiv(IN complex a, IN complex b)
{
  complex t;
  t.im = b.re * b.re + b.im * b.im;
  t.re = (a.re * b.re + a.im * b.im) / t.im;
  t.im = (a.im * b.re - a.re * b.im) / t.im;
  return t;
}

data_t sqrt(IN data_t x)
{
  data_t lo, hi, mid;
  data_t eps = 1e-6;
  if (x < eps)
    return 0;

  if (x < 1.){
    lo = x;
    hi = 1;
  }
  else{
    lo = 1;
    hi = x;
  }
  while (hi - lo > eps){
    mid = (lo + hi) / 2;
    if (mid * mid > x) hi = mid;
    else lo = mid;
  }
  return mid;
}
complex csqrt(IN complex a)
{
  if (a.im == 0){
    if (a.re < 0)
      return cmake(0, sqrt(-a.re));
    else
      return cmake(sqrt(a.re), 0);
  }
  data_t length = sqrt(cabs2(a));
  complex res;
  res.re = sqrt((length + a.re) / 2);
  res.im = sqrt((length - a.re) / 2);
  if (a.im < 0)
    res.re = -res.re;
  return res;
}

/*
 * Complex matrix operations
 */

void cmat_hermitian_transpose_RxTx(
  IN complex A[NUM_RX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex AH[NUM_TX_ANT][NUM_RX_ANT][NUM_SC])
{
  uint32_t i, j, k;
  for (i = 0; i < NUM_RX_ANT; ++i)
    for (j = 0; j < NUM_TX_ANT; ++j)
      for (k = 0; k < NUM_SC; ++k)
        AH[i][j][k] = cconj(A[j][i][k]);
}

void cmat_hermitian_transpose_TxTx(
  IN complex A[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex AH[NUM_TX_ANT][NUM_TX_ANT][NUM_SC])
{
  uint32_t i, j, k;
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < NUM_TX_ANT; ++j)
      for (k = 0; k < NUM_SC; ++k)
        AH[i][j][k] = cconj(A[j][i][k]);
}

void cmatmul_TxRx_RxTx(
  IN complex A[NUM_TX_ANT][NUM_RX_ANT][NUM_SC],
  IN complex B[NUM_RX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_TX_ANT][NUM_SC])
{
  uint32_t i, j, k, l;
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < NUM_TX_ANT; ++j)
      for (k = 0; k < NUM_SC; ++k)
        result[i][j][k] = cmake(0, 0);
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < NUM_TX_ANT; ++j)
      for (k = 0; k < NUM_RX_ANT; ++k)
        for (l = 0; l < NUM_SC; ++l)
          result[i][j][l] = cadd(result[i][j][l], cmul(A[i][k][l], B[k][j][l]));
}

void cmatmul_TxTx_TxTx(
  IN complex A[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex B[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_TX_ANT][NUM_SC])
{
  uint32_t i, j, k, l;
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < NUM_TX_ANT; ++j)
      for (k = 0; k < NUM_SC; ++k)
        result[i][j][k] = cmake(0, 0);
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < NUM_TX_ANT; ++j)
      for (k = 0; k < NUM_TX_ANT; ++k)
        for (l = 0; l < NUM_SC; ++l)
          result[i][j][l] = cadd(result[i][j][l], cmul(A[i][k][l], B[k][j][l]));
}

void cmatadd_TxTx(
  IN complex A[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex B[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_TX_ANT][NUM_SC])
{
  uint32_t i, j, k;
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < NUM_TX_ANT; ++j)
      for (k = 0; k < NUM_SC; ++k)
        result[i][j][k] = cadd(A[i][j][k], B[i][j][k]);
}

void ccholesky_TxTx(
  IN complex A[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex L[NUM_TX_ANT][NUM_TX_ANT][NUM_SC])
{
  uint32_t i, j, k, l;
  complex sum;
  /* Init upper triangle with zeros */
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = i+1; j < NUM_TX_ANT; ++j)
      for (k = 0; k < NUM_SC; ++k)
        L[i][j][k] = cmake(0, 0);

  /* Calculate diagonal and lower triangular entries */
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < i+1; ++j)
      for (k = 0; k < NUM_SC; ++k){
        sum = A[i][j][k];
        for (l = 0; l < j; ++l)
          sum = csub(sum, cmul(L[i][l][k], cconj(L[j][l][k])));
        if (i == j)
          L[i][j][k] = csqrt(sum);
        else
          L[i][j][k] = cdiv(sum, L[j][j][k]);
      }
}

void ccholesky_nosqrt_TxTx(
  IN complex A[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex L[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex D[NUM_TX_ANT][NUM_TX_ANT][NUM_SC])
{
  uint32_t i, j, k, l;
  complex sum, t;
  /* Init zeros in upper triangles of L and D */
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = i+1; j < NUM_TX_ANT; ++j)
      for (k = 0; k < NUM_SC; ++k)
        D[i][j][k] = L[i][j][k] = cmake(0, 0);
  /* Init zeros in lower triangle of D */
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < i; ++j)
      for (k = 0; k < NUM_SC; ++k)
        D[i][j][k] = cmake(0, 0);
  /* Init ones in the diagonal of the result */
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (k = 0; k < NUM_SC; ++k)
      L[i][i][k] = cmake(1, 0);

  /* Calculate the lower triangle of the result */
  for (i = 0; i < NUM_TX_ANT; ++i){
    /* L_ij = (A_ij - \sum_{l=0}^{j-1} L_il L_jl^* D_k) / D_j */
    for (j = 0; j < i; ++j){
      sum.re = 0;
      sum.im = 0;
      for (l = 0; l < j; ++l)
        for (k = 0; k < NUM_SC; ++k){
          /* TODO there is one sum for all scs what leads to contamination. Fix that */
          t = cmul(L[i][l][k], cconj(L[j][l][k]));
          t = cmul(t, D[l][l][k]);
          cadd_acc(sum, t);
        }
      for (k = 0; k < NUM_SC; ++k){
        L[i][j][k] = csub(A[i][j][k], sum);
        L[i][j][k] = cdiv(L[i][j][k], D[j][j][k]);
      }
    }
    /* D_i = A_ii - \sum_{j=0}^{i-1} L_ij L_ij^* D_j */
    sum.re = 0;
    sum.im = 0;
    for (j = 0; j < i; ++j)
      for (k = 0; k < NUM_SC; ++k){
        /* TODO there is one sum for all scs what leads to contamination. Fix that */
        t = cmul(L[i][j][k], cconj(L[i][j][k]));
        t = cmul(t, D[j][j][k]);
        cadd_acc(sum, t);
      }
    for (k = 0; k < NUM_SC; ++k)
      D[i][i][k] = csub(A[i][i][k], sum);
  }
}

void cmatvecmul_TxRx(
  IN complex A[NUM_TX_ANT][NUM_RX_ANT][NUM_SC],
  IN complex b[NUM_RX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_SC])
{
  uint32_t i, j, k;
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < NUM_SC; ++j)
      result[i][j] = cmake(0., 0.);
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < NUM_RX_ANT; ++j)
      for (k = 0; k < NUM_SC; ++k)
        result[i][k] = cadd(result[i][k], cmul(A[i][j][k], b[j][k]));
}

void cforwardsub_TxTx(
  IN complex L[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex b[NUM_TX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_SC])
{
  uint32_t i, j, k;
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < NUM_SC; ++j){
      result[i][j] = b[i][j];
      for (k = 0; k < i; ++k)
        result[i][j] = csub(result[i][j], cmul(L[i][k][j], result[k][j]));
      result[i][j] = cdiv(result[i][j], L[i][i][j]);
    }
}

void cbackwardsub_TxTx(
  IN complex U[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex b[NUM_TX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_SC])
{
  uint32_t i, j, k;
  for (i = NUM_TX_ANT; i != (uint32_t)-1; --i)
    for (j = 0; j < NUM_SC; ++j){
      result[i][j] = b[i][j];
      for (k = i+1; k < NUM_TX_ANT; ++k)
        result[i][j] = csub(result[i][j], cmul(U[i][k][j], result[k][j]));
      result[i][j] = cdiv(result[i][j], U[i][i][j]);
    }
}

/*
 * MMSE
 */

void mmse(
  IN complex H[NUM_RX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex y[NUM_RX_ANT][NUM_SC],
  IN complex R[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex x_MMSE[NUM_TX_ANT][NUM_SC])
{
  complex HH[NUM_TX_ANT][NUM_RX_ANT][NUM_SC];
  complex HH_H[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  complex L[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  complex HHy[NUM_TX_ANT][NUM_SC];
  complex z[NUM_TX_ANT][NUM_SC];
  complex LH[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];

  /* H^H */
  TIME(
    "Hermitian transpose (RxTx): %ld\n",
    cmat_hermitian_transpose_RxTx,
    H, HH);
  /* H^H*H */
  TIME(
    "Matmul (TxRx x RxTx): %ld\n",
    cmatmul_TxRx_RxTx,
    HH, H, HH_H);
  /* H^H*H + R */
  TIME(
    "Matadd (TxTx + TxTx): %ld\n",
    cmatadd_TxTx,
    HH_H, R, HH_H);
  /* L: (H^H*H + R) = L*L^H */
  TIME(
    "Cholesky (TxTx): %ld\n",
    ccholesky_TxTx,
    HH_H, L);
  /* z: L*z = H^H*y */
  TIME(
    "Matrix-vector multiplication (TxRx x Rx): %ld\n",
    cmatvecmul_TxRx,
    HH, y, HHy);
  TIME(
    "Forward substitution (TxTx): %ld\n",
    cforwardsub_TxTx,
    L, HHy, z);
  /* x_MMSE: L^H*x_MMSE = z */
  TIME(
    "Hermitian transpose (TxTx): %ld\n",
    cmat_hermitian_transpose_TxTx,
    L, LH);
  TIME(
    "Backward substitution (TxTx): %ld\n",
    cbackwardsub_TxTx,
    LH, z, x_MMSE);
}

void mmse_nosqrt(
  IN complex H[NUM_RX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex y[NUM_RX_ANT][NUM_SC],
  IN complex R[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex x_MMSE[NUM_TX_ANT][NUM_SC])
{
  /* H^H */
  complex HH[NUM_TX_ANT][NUM_RX_ANT][NUM_SC];
  complex HH_H[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  complex L[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  /* TODO D should be a vector to save memory */
  complex D[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  complex HHy[NUM_TX_ANT][NUM_SC];
  complex z[NUM_TX_ANT][NUM_SC];
  complex LH[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];

  TIME(
    "Hermitian transpose (RxTx): %ld\n",
    cmat_hermitian_transpose_RxTx,
    H, HH);
  TIME(
    "Matmul (TxRx x RxTx): %ld\n",
    cmatmul_TxRx_RxTx,
    HH, H, HH_H);
  TIME(
    "Matadd (TxTx + TxTx): %ld\n",
    cmatadd_TxTx,
    HH_H, R, HH_H);
  /* L, D: (H^H*H + R) = L*D*L^H */
  TIME(
    "Cholesky nosqrt (TxTx): %ld\n",
    ccholesky_nosqrt_TxTx,
    HH_H, L, D);
  /* z: L*z = H^H*y */
  TIME(
    "Matrix-vector multiplication (TxRx x Rx): %ld\n",
    cmatvecmul_TxRx,
    HH, y, HHy);
  TIME(
    "Forward substitution (TxTx): %ld\n",
    cforwardsub_TxTx,
    L, HHy, z);
  /* x_MMSE: L^H*x_MMSE = D^-1*z */
  start_timer();
  for (uint32_t i = 0; i != NUM_TX_ANT; ++i)
    for (uint32_t j = 0; j != NUM_SC; ++j)
      z[i][j] = cdiv(z[i][j], D[i][i][j]);
  stop_timer();
  printf("Division (TxTx): %ld\n", get_timer());
  TIME(
    "Hermitian transpose (TxTx): %ld\n",
    cmat_hermitian_transpose_TxTx,
    L, LH);
  TIME(
    "Backward substitution (TxTx): %ld\n",
    cbackwardsub_TxTx,
    LH, z, x_MMSE);
}

acc_t mse(
  IN complex x[NUM_TX_ANT][NUM_SC],
  IN complex x_MMSE[NUM_TX_ANT][NUM_SC])
{
  acc_t sum = 0;
  for (uint32_t i = 0; i != NUM_TX_ANT; ++i)
    for (uint32_t j = 0; j != NUM_SC; ++j)
      sum += cabs2(csub(x[i][j], x_MMSE[i][j]));
  return sum / (NUM_TX_ANT * NUM_SC);
}

acc_t mse(
  IN complex x[NUM_TX_ANT][NUM_SC],
  IN complex x_MMSE[NUM_TX_ANT][NUM_SC])
{
  acc_t sum = 0;
  for (uint32_t i = 0; i != NUM_TX_ANT; ++i)
    for (uint32_t j = 0; j != NUM_SC; ++j)
      sum += cabs2(csub(x[i][j], x_MMSE[i][j]));
  return sum / (NUM_TX_ANT * NUM_SC);
}
