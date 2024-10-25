#include "mmserv.h"

#include <stdlib.h> // for exit

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
acc_t cabs2(IN complex a)
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

data_t sqrt(IN acc_t x)
{
  if (x < 2) return x;
  acc_t lo = 1, hi = x, mid;
  while (100 * lo * lo < x) lo *= 10;
  while (hi * hi / 100 > x) hi /= 10;
  for (int i = 0; i < 100; i++) {
    mid = (lo + hi) / 2;
    if (mid * mid == x) return mid;
    if (mid * mid > x) hi = mid;
    else lo = mid;
  }
  return mid;
}
complex csqrt(IN complex a)
{
  data_t length = sqrt(cabs2(a));
  complex res;
  res.re = sqrt((length + a.re) / 2);
  res.im = sqrt((length - a.re) / 2);
  res.im *= (a.im > 0) - (a.im < 0);
  return res;
}

/*
 * Complex matrix operations
 */

void cmat_hermitian_transpose_RxTx(
  IN complex A[NUM_RX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex AH[NUM_TX_ANT][NUM_RX_ANT][NUM_SC])
{
  for (uint32_t i = 0; i < NUM_RX_ANT; ++i)
    for (uint32_t j = 0; j < NUM_TX_ANT; ++j)
      for (uint32_t k = 0; k < NUM_SC; ++k)
        AH[i][j][k] = cconj(A[j][i][k]);
}

void cmat_hermitian_transpose_TxTx(
  IN complex A[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex AH[NUM_TX_ANT][NUM_TX_ANT][NUM_SC])
{
  for (uint32_t i = 0; i < NUM_TX_ANT; ++i)
    for (uint32_t j = 0; j < NUM_TX_ANT; ++j)
      for (uint32_t k = 0; k < NUM_SC; ++k)
        AH[i][j][k] = cconj(A[j][i][k]);
}

void cmatmul_TxRx_RxTx(
  IN complex A[NUM_TX_ANT][NUM_RX_ANT][NUM_SC],
  IN complex B[NUM_RX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_TX_ANT][NUM_SC])
{
  for (uint32_t i = 0; i < NUM_TX_ANT; ++i)
    for (uint32_t j = 0; j < NUM_TX_ANT; ++j)
      for (uint32_t k = 0; k < NUM_RX_ANT; ++k)
        for (uint32_t l = 0; l < NUM_SC; ++l)
          result[i][j][l] = cadd(result[i][j][l], cmul(A[i][k][l], B[k][j][l]));
}

void cmatadd_TxTx(
  IN complex A[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex B[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_TX_ANT][NUM_SC])
{
  for (uint32_t i = 0; i < NUM_TX_ANT; ++i)
    for (uint32_t j = 0; j < NUM_TX_ANT; ++j)
      for (uint32_t k = 0; k < NUM_SC; ++k)
        result[i][j][k] = cadd(A[i][j][k], B[i][j][k]);
}

void ccholesky_TxTx(
  IN complex A[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex L[NUM_TX_ANT][NUM_TX_ANT][NUM_SC])
{
  /* Init upper triangle with zeros */
  for (uint32_t i = 0; i < NUM_TX_ANT; ++i)
    for (uint32_t j = i+1; j < NUM_TX_ANT; ++j)
      for (uint32_t k = 0; k < NUM_SC; ++k)
        L[i][j][k] = cmake(0, 0);
  /* Calculate diagonal and lower triangular entries */
  for (uint32_t i = 0; i < NUM_TX_ANT; ++i)
    for (uint32_t j = 0; j < i+1; ++j)
      for (uint32_t k = 0; k < NUM_SC; ++k){
        complex sum = {0, 0};
        for (uint32_t l = 0; l < j; ++l)
          cadd_acc(sum, cmul(L[i][l][k], cconj(L[j][l][k])));
        if (i == j)
          L[i][j][k] = csqrt(csub(A[i][i][k], sum));
        else
          L[i][j][k] = cdiv(csub(A[i][j][k], sum), L[j][j][k]);

        // TODO delete
        #if 0
        for (uint32_t i = 0; i < NUM_TX_ANT; ++i){
          for (uint32_t j = 0; j < NUM_TX_ANT; ++j)
            printf("(%6i, %6i) ", L[i][j][0].re, L[i][j][0].im);
          printf("\n");
        }
        printf("\n");
        #endif
      }
}

void ccholesky_nosqrt_TxTx(
  IN complex A[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex L[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex D[NUM_TX_ANT][NUM_TX_ANT][NUM_SC])
{
  /* Init zeros in upper triangles of L and D */
  for (uint32_t i = 0; i < NUM_TX_ANT; ++i)
    for (uint32_t j = i+1; j < NUM_TX_ANT; ++j)
      for (uint32_t k = 0; k < NUM_SC; ++k)
        D[i][j][k] = L[i][j][k] = cmake(0, 0);
  /* Init zeros in lower triangle of D */
  for (uint32_t i = 0; i < NUM_TX_ANT; ++i)
    for (uint32_t j = 0; j < i; ++j)
      for (uint32_t k = 0; k < NUM_SC; ++k)
        D[i][j][k] = cmake(0, 0);
  /* Init ones in the diagonal of the result */
  for (uint32_t i = 0; i < NUM_TX_ANT; ++i)
    for (uint32_t k = 0; k < NUM_SC; ++k)
      L[i][i][k] = cmake(1, 0);


  /* Calculate the lower triangle of the result */
  for (uint32_t i = 0; i < NUM_TX_ANT; ++i){
    /* L_ij = (A_ij - \sum_{l=0}^{j-1} L_il L_jl^* D_k) / D_j */
    for (uint32_t j = 0; j < i; ++j){
      complex sum = {0, 0};
      for (uint32_t l = 0; l < j; ++l)
        for (uint32_t k = 0; k < NUM_SC; ++k){
          complex t = cmul(L[i][l][k], cconj(L[j][l][k]));
          t = cmul(t, D[l][l][k]);
          cadd_acc(sum, t);
        }
      for (uint32_t k = 0; k < NUM_SC; ++k){
        L[i][j][k] = csub(A[i][j][k], sum);
        L[i][j][k] = cdiv(L[i][j][k], D[j][j][k]);
      }
    }
    /* D_i = A_ii - \sum_{j=0}^{i-1} L_ij L_ij^* D_j */
    complex sum = {0, 0};
    for (uint32_t j = 0; j < i; ++j)
      for (uint32_t k = 0; k < NUM_SC; ++k){
        complex t = cmul(L[i][j][k], cconj(L[i][j][k]));
        t = cmul(t, D[j][j][k]);
        cadd_acc(sum, t);
      }
    for (uint32_t k = 0; k < NUM_SC; ++k)
      D[i][i][k] = csub(A[i][i][k], sum);
  }
}

void cmatvecmul_TxRx(
  IN complex A[NUM_TX_ANT][NUM_RX_ANT][NUM_SC],
  IN complex b[NUM_RX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_SC])
{
  for (uint32_t i = 0; i < NUM_TX_ANT; ++i)
    for (uint32_t j = 0; j < NUM_RX_ANT; ++j)
      for (uint32_t k = 0; k < NUM_SC; ++k)
        result[i][k] = cadd(result[i][k], cmul(A[i][j][k], b[j][k]));
}

void cforwardsub_TxTx(
  IN complex L[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex b[NUM_TX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_SC])
{
  for (uint32_t i = 0; i < NUM_TX_ANT; ++i)
    for (uint32_t j = 0; j < NUM_SC; ++j){
      result[i][j] = b[i][j];
      for (uint32_t k = 0; k < i; ++k)
        result[i][j] = csub(result[i][j], cmul(L[i][k][j], result[k][j]));
      result[i][j] = cdiv(result[i][j], L[i][i][j]);
    }
}

void cbackwardsub_TxTx(
  IN complex U[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex b[NUM_TX_ANT][NUM_SC],
  OUT complex result[NUM_TX_ANT][NUM_SC])
{
  for (uint32_t i = NUM_TX_ANT; i > 0; --i)
    for (uint32_t j = 0; j < NUM_SC; ++j){
      result[i][j] = b[i][j];
      for (uint32_t k = i+1; k < NUM_TX_ANT; ++k)
        result[i][j] = csub(result[i][j], cmul(U[i][k][j], result[k][j]));
      result[i][j] = cdiv(result[i][j], U[i][i][j]);
    }
}

/*
 * IO
 */

void load_data(
  IN const char* file,
  IN size_t size,
  OUT uint16_t *out)
{
    FILE *fp = fopen(file, "rb");
    
    if (fp == NULL) {
        fprintf(stderr, "Error opening file %s.\n", file);
        exit(8);
    }
    if ((fread(out, sizeof(uint16_t), size, fp)) != size) {
        fprintf(stderr, "Error reading file %s.\n", file);
        exit(8);
    }

    fclose(fp);
}

void save_data(
  IN const char* file,
  IN uint16_t x_mmse[NUM_TX_ANT][NUM_SC][2])
{
    FILE *fp = fopen(file, "wb");
    
    if (fp == NULL) {
        fprintf(stderr, "Error opening file %s.\n", file);
        exit(8);
    }
    if ((fwrite(x_mmse, sizeof(uint16_t), NUM_TX_ANT*NUM_SC*2, fp)) != NUM_TX_ANT*NUM_SC*2) {
        fprintf(stderr, "Error writing file %s.\n", file);
        exit(8);
    }

    fclose(fp);
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
  /* H^H */
  complex HH[NUM_TX_ANT][NUM_RX_ANT][NUM_SC];
  cmat_hermitian_transpose_RxTx(H, HH);
  /* H^H*H */
  complex HH_H[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  cmatmul_TxRx_RxTx(HH, H, HH_H);
  /* H^H*H + R */
  cmatadd_TxTx(HH_H, R, HH_H);
  /* L: (H^H*H + R) = L*L^H */
  complex L[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  ccholesky_TxTx(HH_H, L);
  /* z: L*z = H^H*y */
  complex HHy[NUM_TX_ANT][NUM_SC];
  cmatvecmul_TxRx(HH, y, HHy);
  complex z[NUM_TX_ANT][NUM_SC];
  cforwardsub_TxTx(L, HHy, z);
  /* x_MMSE: L^H*x_MMSE = z */
  complex LH[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  cmat_hermitian_transpose_TxTx(L, LH);
  cbackwardsub_TxTx(LH, z, x_MMSE);
}

void mmse_nosqrt(
  IN complex H[NUM_RX_ANT][NUM_TX_ANT][NUM_SC],
  IN complex y[NUM_RX_ANT][NUM_SC],
  IN complex R[NUM_TX_ANT][NUM_TX_ANT][NUM_SC],
  OUT complex x_MMSE[NUM_TX_ANT][NUM_SC])
{
  /* H^H */
  complex HH[NUM_TX_ANT][NUM_RX_ANT][NUM_SC];
  cmat_hermitian_transpose_RxTx(H, HH);
  /* H^H*H */
  complex HH_H[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  cmatmul_TxRx_RxTx(HH, H, HH_H);
  /* H^H*H + R */
  cmatadd_TxTx(HH_H, R, HH_H);
  /* L, D: (H^H*H + R) = L*D*L^H */
  complex L[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  complex D[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  ccholesky_nosqrt_TxTx(HH_H, L, D);
  /* z: L*z = H^H*y */
  complex HHy[NUM_TX_ANT][NUM_SC];
  cmatvecmul_TxRx(HH, y, HHy);
  complex z[NUM_TX_ANT][NUM_SC];
  cforwardsub_TxTx(L, HHy, z);
  /* x_MMSE: L^H*x_MMSE = D^-1*z */
  complex LH[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  cmat_hermitian_transpose_TxTx(L, LH);
  cbackwardsub_TxTx(LH, z, x_MMSE);
  cmatvecmul_TxRx(D, x_MMSE, x_MMSE);
}
