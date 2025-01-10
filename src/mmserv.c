#include "../include/mmserv.h"

#include "../../common/rivec/vector_defines.h"
#include "riscv_vector.h"

#include <stddef.h> /* for size_t */

#ifdef DEBUG
#include "../../common/runtime.h"
#include "../../common/util.h"
#include "printf.h"

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

void cdiv(
  IN data_t a_re, IN data_t a_im,
  IN data_t b_re, IN data_t b_im,
  OUT data_t *res_re, OUT data_t *res_im)
{
  *res_im = b_re * b_re + b_im * b_im;
  *res_re = (a_re * b_re + a_im * b_im) / *res_im;
  *res_im = (a_im * b_re - a_re * b_im) / *res_im;
}

data_t sqrt(IN data_t x)
{
  data_t lo, hi, mid;
  data_t eps = 1e-4;
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
    mid = (lo + hi) / 2.f;
    if (mid * mid > x) hi = mid;
    else lo = mid;
  }
  return mid;
}
void csqrt(
  IN data_t re, IN data_t im,
  OUT data_t *res_re, OUT data_t *res_im)
{
  if (im == 0){
    if (re < 0) {
      *res_re = 0;
      *res_im = sqrt(-re);
      return;
    }
    else {
      *res_re = sqrt(re);
      *res_im = 0;
      return;
    }
  }
  data_t length = sqrt(re * re + im * im);
  *res_re = sqrt((length + re) / 2.f);
  *res_im = sqrt((length - re) / 2.f);
  if (im < 0.f)
    *res_re = -*res_re;
}

/*
 * Complex matrix operations
 */

void cmat_hermitian_transpose_RxTx(
  IN vcomplex *A,
  OUT vcomplex *AH)
{
  size_t i, j, k, off_ijk = 0, off_jik;
  for (i = 0; i < NUM_RX_ANT; ++i)
    for (j = 0; j < NUM_TX_ANT; ++j) {
      off_jik = j * NUM_RX_ANT * NUM_SC + i * NUM_SC;
      for (k = 0; k < NUM_SC; ++k) {
        AH->re[off_jik] = A->re[off_ijk];
        AH->im[off_jik] = -A->im[off_ijk];
        ++off_ijk;
        ++off_jik;
      }
    }
}

void cmat_hermitian_transpose_TxTx(
  IN vcomplex *A,
  OUT vcomplex *AH)
{
  size_t i, j, k, off_ijk = 0, off_jik;
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < NUM_TX_ANT; ++j) {
      off_jik = j * NUM_TX_ANT * NUM_SC + i * NUM_SC;
      for (k = 0; k < NUM_SC; ++k) {
        AH->re[off_jik] = A->re[off_ijk];
        AH->im[off_jik] = -A->im[off_ijk];
        ++off_ijk;
        ++off_jik;
      }
    }
}

void cmatmul_TxRx_RxTx(
  IN vcomplex *A,
  IN vcomplex *B,
  OUT vcomplex *result)
{
  size_t i, j, k, l;
  size_t off_ijk = 0, off_ilk, off_ljk;
  data_t A_re, A_im, B_re, B_im;

  for (i = 0; i < NUM_TX_ANT * NUM_TX_ANT * NUM_SC; ++i)
    result->re[i] = result->im[i] = 0.f;

  for (i = 0; i != NUM_TX_ANT; ++i) {
    for (j = 0; j != NUM_RX_ANT; ++j) {
      for (k = 0; k != NUM_SC; ++k) {
        off_ilk = i * NUM_RX_ANT * NUM_SC + k;
        off_ljk = j * NUM_SC + k;
        for (l = 0; l != NUM_RX_ANT; ++l) {
          A_re = A->re[off_ilk];
          A_im = A->im[off_ilk];
          B_re = B->re[off_ljk];
          B_im = B->im[off_ljk];
          result->re[off_ijk] += A_re * B_re - A_im * B_im;
          result->im[off_ijk] += A_re * B_im + A_im * B_re;
          off_ilk += NUM_SC;
          off_ljk += NUM_TX_ANT * NUM_SC;
        }
        ++off_ijk;
      }
    }
  }
}

void cmatmul_TxTx_TxTx(
  IN vcomplex *A,
  IN vcomplex *B,
  OUT vcomplex *result)
{
  size_t i, j, k, l;
  size_t off_ijk = 0, off_ilk, off_ljk;
  data_t A_re, A_im, B_re, B_im;

  for (i = 0; i < NUM_TX_ANT * NUM_TX_ANT * NUM_SC; ++i)
    result->re[i] = result->im[i] = 0.f;

  for (i = 0; i != NUM_TX_ANT; ++i) {
    for (j = 0; j != NUM_TX_ANT; ++j) {
      for (k = 0; k != NUM_SC; ++k) {
        off_ilk = i * NUM_TX_ANT * NUM_SC + k;
        off_ljk = j * NUM_SC + k;
        for (l = 0; l != NUM_TX_ANT; ++l) {
          A_re = A->re[off_ilk];
          A_im = A->im[off_ilk];
          B_re = B->re[off_ljk];
          B_im = B->im[off_ljk];
          result->re[off_ijk] += A_re * B_re - A_im * B_im;
          result->im[off_ijk] += A_re * B_im + A_im * B_re;
          off_ilk += NUM_SC;
          off_ljk += NUM_TX_ANT * NUM_SC;
        }
        ++off_ijk;
      }
    }
  }
}

void cmatadd_TxTx(
  IN vcomplex *A,
  IN vcomplex *B,
  OUT vcomplex *result)
{
  vfloat32m1_t vA, vB, vresult;
  size_t vl, sz = NUM_TX_ANT * NUM_TX_ANT * NUM_SC;
  size_t off = 0;
  while (sz > 0) {
    vl = vsetvl_e32m1(sz);

    /* real part */
    vA = vle32_v_f32m1(&A->re[off], vl);
    vB = vle32_v_f32m1(&B->re[off], vl);
    vresult = vfadd_vv_f32m1(vA, vB, vl);
    vse32_v_f32m1(&result->re[off], vresult, vl);

    /* imaginary part */
    vA = vle32_v_f32m1(&A->im[off], vl);
    vB = vle32_v_f32m1(&B->im[off], vl);
    vresult = vfadd_vv_f32m1(vA, vB, vl);
    vse32_v_f32m1(&result->im[off], vresult, vl);

    sz -= vl;
    off += vl;
  }
}

void ccholesky_TxTx(
  IN vcomplex *A,
  OUT vcomplex *L)
{
  size_t i, j, k, l;
  size_t off_ijk, off_ilk, off_jlk, off_jjk;
  data_t sum_re, sum_im;
  data_t t1_re, t1_im, t2_re, t2_im;

  /* Init upper triangle with zeros */
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = i+1; j < NUM_TX_ANT; ++j) {
      off_ijk = i * NUM_TX_ANT * NUM_SC + j * NUM_SC;
      for (k = 0; k < NUM_SC; ++k) {
        L->re[off_ijk] = L->im[off_ijk] = 0;
        ++off_ijk;
      }
    }

  /* Calculate diagonal and lower triangular entries */
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < i+1; ++j)
      for (k = 0; k < NUM_SC; ++k){
        off_ijk = i * NUM_TX_ANT * NUM_SC + j * NUM_SC + k;
        sum_re = A->re[off_ijk];
        sum_im = A->im[off_ijk];
        for (l = 0; l < j; ++l) {
          off_ilk = i * NUM_TX_ANT * NUM_SC + l * NUM_SC + k;
          off_jlk = j * NUM_TX_ANT * NUM_SC + l * NUM_SC + k;
          t1_re = L->re[off_ilk];
          t1_im = L->im[off_ilk];
          t2_re = L->re[off_jlk];
          t2_im = L->im[off_jlk];
          sum_re -= t1_re * t2_re + t1_im * t2_im;
          sum_im -= t1_im * t2_re - t1_re * t2_im;
        }
        if (i == j)
          csqrt(
            sum_re, sum_im,
            &L->re[off_ijk],
            &L->im[off_ijk]);
        else {
          off_jjk = j * NUM_TX_ANT * NUM_SC + j * NUM_SC + k;
          cdiv(
            sum_re, sum_im,
            L->re[off_jjk],
            L->im[off_jjk],
            &L->re[off_ijk],
            &L->im[off_ijk]);
        }
      }
}

void ccholesky_nosqrt_TxTx(
  IN vcomplex *A,
  OUT vcomplex *L,
  OUT vcomplex *D)
{
  size_t i, j, k, l;
  size_t off_iik, off_ijk, off_ilk, off_jlk;
  size_t off_ik, off_jk;
  data_t t1_re, t1_im, t2_re, t2_im;
  data_t sum_re, sum_im;

  /* Init zeros in upper triangles of L and D */
  for (off_ijk = 0; off_ijk < NUM_TX_ANT * NUM_TX_ANT * NUM_SC; ++off_ijk)
    L->re[off_ijk] = L->im[off_ijk] = D->re[off_ijk] = D->im[off_ijk] = 0;
  /* Init zeros in lower triangle of D */
  /* Init ones in the diagonal of the result */
  for (off_ik = 0; off_ik < NUM_TX_ANT * NUM_SC; ++off_ik) {
    L->im[off_ik] = D->re[off_ik] = D->im[off_ik] = 0;
    L->re[off_ik] = 1;
  }

  /* Calculate the lower triangle of the result */
  for (k = 0; k != NUM_SC; ++k) {
    for (i = 0; i != NUM_TX_ANT; ++i) {
      off_iik = i * NUM_TX_ANT * NUM_SC + i * NUM_SC + k;
      /* L_ij = (A_ij - \sum_{l=0}^{j-1} L_il L_jl^* D_k) / D_j */
      for (j = 0; j != i; ++j) {
        for (l = 0; l != j; ++l) {
          off_ilk = i * NUM_TX_ANT * NUM_SC + l * NUM_SC + k;
          off_jlk = j * NUM_TX_ANT * NUM_SC + l * NUM_SC + k;
          t1_re = L->re[off_ilk];
          t1_im = L->im[off_ilk];
          t2_re = L->re[off_jlk];
          t2_im = L->im[off_jlk];
          sum_re += t1_re * t2_re + t1_im * t2_im;
          sum_im += t1_im * t2_re - t1_re * t2_im;
        }
        off_ijk = i * NUM_TX_ANT * NUM_SC + j * NUM_SC + k;
        off_jk = j * NUM_SC + k;
        t1_re = A->re[off_ijk] - sum_re;
        t1_im = A->im[off_ijk] - sum_im;
        cdiv(
          t1_re, t1_im,
          D->re[off_jk], D->im[off_jk],
          &L->re[off_ijk],
          &L->im[off_ijk]);
      }
      /* D_i = A_ii - \sum_{j=0}^{i-1} L_ij L_ij^* D_j */
      sum_re = sum_im = 0;
      off_ik = i * NUM_SC + k;
      for (j = 0; j != i; ++j) {
        off_ijk = i * NUM_TX_ANT * NUM_SC + j * NUM_SC + k;
        off_jk = j * NUM_SC + k;
        t1_re = L->re[off_ijk];
        t1_im = L->im[off_ijk];
        t2_re = t1_re * t1_re + t1_im * t1_im;
        t2_im = t1_re * t1_im - t1_im * t1_re;
        t1_re = D->re[off_jk];
        t1_im = D->im[off_jk];
        sum_re += t2_re * t1_re - t2_im * t1_im;
        sum_im += t2_re * t1_im + t2_im * t1_re;
      }
      D->re[off_ik] = A->re[off_iik] - sum_re;
      D->im[off_ik] = A->im[off_iik] - sum_im;
    }
  }
}

void cmatvecmul_TxRx(
  IN vcomplex *A,
  IN vcomplex *b,
  OUT vcomplex *result)
{
  size_t i, j, k;
  size_t off_ik, off_ijk = 0, off_jk;
  size_t off_ik_bck;
  data_t A_re, A_im, b_re, b_im;

  for (i = 0; i < NUM_TX_ANT * NUM_SC; ++i)
    result->re[i] = result->im[i] = 0.f;

  for (i = 0; i < NUM_TX_ANT; ++i) {
    off_jk = 0;
    off_ik_bck = i * NUM_SC;
    for (j = 0; j < NUM_RX_ANT; ++j) {
      off_ik = off_ik_bck;
      for (k = 0; k < NUM_SC; ++k) {
        A_re = A->re[off_ijk];
        A_im = A->im[off_ijk];
        b_re = b->re[off_jk];
        b_im = b->im[off_jk];
        result->re[off_ik] += A_re * b_re - A_im * b_im;
        result->im[off_ik] += A_re * b_im + A_im * b_re;
        ++off_ik;
        ++off_ijk;
        ++off_jk;
      }
    }
  }
}

void cforwardsub_TxTx(
  IN vcomplex *L,
  IN vcomplex *b,
  OUT vcomplex *result)
{
  size_t i, j, k;
  size_t off_ij = 0, off_ikj, off_kj, off_iij;
  size_t off_ikj_bck;
  data_t L_re, L_im, result_re, result_im;

  for (i = 0; i < NUM_TX_ANT; ++i) {
    off_ij = i * NUM_SC;
    off_ikj_bck = i * NUM_TX_ANT * NUM_SC;
    off_iij = i * NUM_TX_ANT * NUM_SC + i * NUM_SC;
    for (j = 0; j < NUM_SC; ++j) {
      result->re[off_ij] = b->re[off_ij];
      result->im[off_ij] = b->im[off_ij];

      off_ikj = off_ikj_bck;
      off_kj = j;
      for (k = 0; k < i; ++k) {
        L_re = L->re[off_ikj];
        L_im = L->im[off_ikj];
        result_re = result->re[off_kj];
        result_im = result->im[off_kj];
        result->re[off_ij] -= L_re * result_re - L_im * result_im;
        result->im[off_ij] -= L_re * result_im + L_im * result_re;
        off_kj += NUM_SC;
        off_ikj += NUM_SC;
      }

      result_re = result->re[off_ij];
      result_im = result->im[off_ij];
      cdiv(
        result_re, result_im,
        L->re[off_iij], L->im[off_iij],
        &result->re[off_ij],
        &result->im[off_ij]);

      ++off_ij;
      ++off_ikj_bck;
      ++off_iij;
    }
  }
}

void cbackwardsub_TxTx(
  IN vcomplex *U,
  IN vcomplex *b,
  OUT vcomplex *result)
{
  /* TODO optimize offsets */

  size_t i, j, k;
  size_t off_ij, off_ikj, off_kj, off_iij;
  data_t U_re, U_im, result_re, result_im;

  for (i = NUM_TX_ANT - 1; i != (size_t)-1; --i) {
    for (j = 0; j < NUM_SC; ++j){
      off_ij = i * NUM_SC + j;
      off_iij = i * NUM_TX_ANT * NUM_SC + i * NUM_SC + j;
      result->re[off_ij] = b->re[off_ij];
      result->im[off_ij] = b->im[off_ij];

      for (k = i+1; k < NUM_TX_ANT; ++k){
        off_ikj = i * NUM_TX_ANT * NUM_SC + k * NUM_SC + j;
        off_kj = k * NUM_SC + j;
        U_re = U->re[off_ikj];
        U_im = U->im[off_ikj];
        result_re = result->re[off_kj];
        result_im = result->im[off_kj];
        result->re[off_ij] -= U_re * result_re - U_im * result_im;
        result->im[off_ij] -= U_re * result_im + U_im * result_re;
      }

      result_re = result->re[off_ij];
      result_im = result->im[off_ij];
      cdiv(
        result_re, result_im,
        U->re[off_iij], U->im[off_iij],
        &result->re[off_ij],
        &result->im[off_ij]);
    }
  }
}

/*
 * MMSE
 */

void mmse(
  IN vcomplex *H,
  IN vcomplex *y,
  IN vcomplex *R,
  OUT vcomplex *x_MMSE)
{
  vcomplex HH, HH_H, L, HHy, z, LH;

  data_t HH_re[NUM_TX_ANT][NUM_RX_ANT][NUM_SC];
  data_t HH_im[NUM_TX_ANT][NUM_RX_ANT][NUM_SC];
  data_t HH_H_re[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  data_t HH_H_im[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  data_t L_re[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  data_t L_im[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  data_t HHy_re[NUM_TX_ANT][NUM_SC];
  data_t HHy_im[NUM_TX_ANT][NUM_SC];
  data_t z_re[NUM_TX_ANT][NUM_SC];
  data_t z_im[NUM_TX_ANT][NUM_SC];
  data_t LH_re[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  data_t LH_im[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];

  HH.re = (data_t *)HH_re;
  HH.im = (data_t *)HH_im;
  HH_H.re = (data_t *)HH_H_re;
  HH_H.im = (data_t *)HH_H_im;
  L.re = (data_t *)L_re;
  L.im = (data_t *)L_im;
  HHy.re = (data_t *)HHy_re;
  HHy.im = (data_t *)HHy_im;
  z.re = (data_t *)z_re;
  z.im = (data_t *)z_im;
  LH.re = (data_t *)LH_re;
  LH.im = (data_t *)LH_im;

  /* H^H */
  TIME(
    "Hermitian transpose (RxTx): %ld\n",
    cmat_hermitian_transpose_RxTx,
    H, &HH);
  /* H^H*H */
  TIME(
    "Matmul (TxRx x RxTx): %ld\n",
    cmatmul_TxRx_RxTx,
    &HH, H, &HH_H);
  /* H^H*H + R */
  TIME(
    "Matadd (TxTx + TxTx): %ld\n",
    cmatadd_TxTx,
    &HH_H, R, &HH_H);
  /* L: (H^H*H + R) = L*L^H */
  TIME(
    "Cholesky (TxTx): %ld\n",
    ccholesky_TxTx,
    &HH_H, &L);
  /* z: L*z = H^H*y */
  TIME(
    "Matrix-vector multiplication (TxRx x Rx): %ld\n",
    cmatvecmul_TxRx,
    &HH, y, &HHy);
  TIME(
    "Forward substitution (TxTx): %ld\n",
    cforwardsub_TxTx,
    &L, &HHy, &z);
  /* x_MMSE: L^H*x_MMSE = z */
  TIME(
    "Hermitian transpose (TxTx): %ld\n",
    cmat_hermitian_transpose_TxTx,
    &L, &LH);
  TIME(
    "Backward substitution (TxTx): %ld\n",
    cbackwardsub_TxTx,
    &LH, &z, x_MMSE);
}

void mmse_nosqrt(
  IN vcomplex *H,
  IN vcomplex *y,
  IN vcomplex *R,
  OUT vcomplex *x_MMSE)
{
  vcomplex HH, HH_H, L, D, HHy, z, LH;

  data_t HH_re[NUM_TX_ANT][NUM_RX_ANT][NUM_SC];
  data_t HH_im[NUM_TX_ANT][NUM_RX_ANT][NUM_SC];
  data_t HH_H_re[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  data_t HH_H_im[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  data_t L_re[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  data_t L_im[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  data_t D_re[NUM_TX_ANT][NUM_SC];
  data_t D_im[NUM_TX_ANT][NUM_SC];
  data_t HHy_re[NUM_TX_ANT][NUM_SC];
  data_t HHy_im[NUM_TX_ANT][NUM_SC];
  data_t z_re[NUM_TX_ANT][NUM_SC];
  data_t z_im[NUM_TX_ANT][NUM_SC];
  data_t LH_re[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
  data_t LH_im[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];

  HH.re = (data_t *)HH_re;
  HH.im = (data_t *)HH_im;
  HH_H.re = (data_t *)HH_H_re;
  HH_H.im = (data_t *)HH_H_im;
  L.re = (data_t *)L_re;
  L.im = (data_t *)L_im;
  D.re = (data_t *)D_re;
  D.im = (data_t *)D_im;
  HHy.re = (data_t *)HHy_re;
  HHy.im = (data_t *)HHy_im;
  z.re = (data_t *)z_re;
  z.im = (data_t *)z_im;
  LH.re = (data_t *)LH_re;
  LH.im = (data_t *)LH_im;

  /* H^H */
  TIME(
    "Hermitian transpose (RxTx): %ld\n",
    cmat_hermitian_transpose_RxTx,
    H, &HH);
  TIME(
    "Matmul (TxRx x RxTx): %ld\n",
    cmatmul_TxRx_RxTx,
    &HH, H, &HH_H);
  TIME(
    "Matadd (TxTx + TxTx): %ld\n",
    cmatadd_TxTx,
    &HH_H, R, &HH_H);
  /* L, D: (H^H*H + R) = L*D*L^H */
  TIME(
    "Cholesky nosqrt (TxTx): %ld\n",
    ccholesky_nosqrt_TxTx,
    &HH_H, &L, &D);
  /* z: L*z = H^H*y */
  TIME(
    "Matrix-vector multiplication (TxRx x Rx): %ld\n",
    cmatvecmul_TxRx,
    &HH, y, &HHy);
  TIME(
    "Forward substitution (TxTx): %ld\n",
    cforwardsub_TxTx,
    &L, &HHy, &z);
  /* x_MMSE: L^H*x_MMSE = D^-1*z */
  start_timer();
  data_t z_t_re, z_t_im;
  for (size_t off_ij = 0; off_ij < NUM_TX_ANT * NUM_SC; ++off_ij) {
    z_t_re = z.re[off_ij];
    z_t_im = z.im[off_ij];
    cdiv(
      z_t_re, z_t_im,
      D.re[off_ij],
      D.im[off_ij],
      &z.re[off_ij],
      &z.im[off_ij]);
  }
  stop_timer();
  printf("Division (Tx): %ld\n", get_timer());
  TIME(
    "Hermitian transpose (TxTx): %ld\n",
    cmat_hermitian_transpose_TxTx,
    &L, &LH);
  TIME(
    "Backward substitution (TxTx): %ld\n",
    cbackwardsub_TxTx,
    &LH, &z, x_MMSE);
}

acc_t mse(
  IN vcomplex *x,
  IN vcomplex *x_MMSE)
{
  acc_t sum = 0.;
  size_t off_ij = 0;
  data_t sub1, sub2;
  for (; off_ij != NUM_TX_ANT * NUM_SC; ++off_ij) {
    sub1 = x->re[off_ij] - x_MMSE->re[off_ij];
    sub2 = x->im[off_ij] - x_MMSE->im[off_ij];
    sum += (sub1 * sub1 + sub2 * sub2) / (data_t)(NUM_TX_ANT * NUM_SC);
  }
  return sum;
}
