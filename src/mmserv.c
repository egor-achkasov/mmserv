#include "../include/define.h"

#include <stddef.h> /* for size_t */
#include <stdint.h> /* for uint64_t */

/*
 * Debug
 */

#ifdef DEBUG
#include "printf.h"

static uint64_t g_timer;
void start_timer()
{
  __asm__ volatile("rdcycle %0" : "=r"(g_timer));
}
void stop_timer()
{
  __asm__ volatile(
    "rdcycle t0\n"
    "sub %0, t0, %0"
    : "+r"(g_timer)
  );
}
uint64_t get_timer()
{
  return g_timer;
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
 * Global variables
 */

/* Externs */
extern vcomplex g_x, g_H, g_R, g_y;
extern vcomplex g_x_MMSE;

/* Raw data */
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

/* Same data but casted to vcomplex */
vcomplex g_HH = { .re = (data_t *)HH_re, .im = (data_t *)HH_im };
vcomplex g_HH_H = { .re = (data_t *)HH_H_re, .im = (data_t *)HH_H_im };
vcomplex g_L = { .re = (data_t *)L_re, .im = (data_t *)L_im };
vcomplex g_HHy = { .re = (data_t *)HHy_re, .im = (data_t *)HHy_im };
vcomplex g_z = { .re = (data_t *)z_re, .im = (data_t *)z_im };
vcomplex g_LH = { .re = (data_t *)LH_re, .im = (data_t *)LH_im };

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
  __asm__ volatile("fsqrt.s %0, %1" : "=f"(x) : "f"(x));
  return x;
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

/** Complex Gram matrix H^H*H and add complex matrix R
 * 
 * HH_H = H^H*H + R
 * 
 * \global g_H matrix of channel coefficients. Shape [NUM_RX_ANT][NUM_TX_ANT][NUM_SC]
 * \global g_R noise covariance matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \global g_HH_H output Gram matrix + R. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 */
void cmatgram_TxRx_cadd()
{
  size_t t1, t2, r;
  size_t sz, vl;
  size_t off_sc, off_A, off_AH;
  size_t off_HH_H_L, off_HH_H_U;
  data_t *A_re, *A_im, *AH_re, *AH_im;
  data_t *R_L_re, *R_L_im, *R_U_re, *R_U_im;
  data_t *HH_H_L_re, *HH_H_L_im, *HH_H_U_re, *HH_H_U_im;

  for (t1 = 0; t1 != NUM_TX_ANT; ++t1)
    for (t2 = t1; t2 != NUM_TX_ANT; ++t2) {
      off_sc = 0;
      off_HH_H_L = t1 * NUM_TX_ANT * NUM_SC + t2 * NUM_SC;
      off_HH_H_U = t2 * NUM_TX_ANT * NUM_SC + t1 * NUM_SC;
      HH_H_L_re = &g_HH_H.re[off_HH_H_L];
      HH_H_L_im = &g_HH_H.im[off_HH_H_L];
      HH_H_U_re = &g_HH_H.re[off_HH_H_U];
      HH_H_U_im = &g_HH_H.im[off_HH_H_U];
      sz = NUM_SC;

      while (sz > 0) {
        /* Initialize HH_H registers */
        /* v0 - HH_H real part */
        /* v1 - HH_H imaginary part */
        __asm__ volatile(
          "vsetvli %0, %1, e32, m1, ta, ma\n"
          "vmv.v.i v0, 0\n"
          "vmv.v.i v1, 0\n"
          : "=r"(vl) : "r"(sz));

        for (r = 0; r != NUM_RX_ANT; ++r) {
          off_A = r * NUM_TX_ANT * NUM_SC + t1 * NUM_SC + off_sc;
          off_AH = r * NUM_TX_ANT * NUM_SC + t2 * NUM_SC + off_sc;
          A_re = &g_H.re[off_A];
          A_im = &g_H.im[off_A];
          AH_re = &g_H.re[off_AH];
          AH_im = &g_H.im[off_AH];

          /* Calculate A^H*A */
          /* v2 - A real part */
          /* v3 - A imaginary part */
          /* v4 - A^H real part */
          /* v5 - A^H imaginary part */
          __asm__ volatile(
            "vle32.v v2, (%0)\n"
            "vle32.v v3, (%1)\n"
            "vle32.v v4, (%2)\n"
            "vle32.v v5, (%3)\n"
            /* real part */
            "vfmacc.vv v0, v2, v4\n"
            "vfmacc.vv v0, v3, v5\n"
            /* imaginary part */
            "vfmacc.vv v1, v3, v4\n"
            "vfnmsac.vv v1, v2, v5\n"
            :
            : "r"(A_re), "r"(A_im), "r"(AH_re), "r"(AH_im)
          );
        }

        /* Add R */
        /* v2 - R real part */
        /* v3 - R imaginary part */
        R_U_re = &g_R.re[off_HH_H_U];
        R_U_im = &g_R.im[off_HH_H_U];
        __asm__ volatile(
          "vle32.v v2, (%0)\n"
          "vle32.v v3, (%1)\n"
          "vfadd.vv v2, v2, v0\n"
          "vfadd.vv v3, v3, v1\n"
          "vse32.v v2, (%2)\n"
          "vse32.v v3, (%3)\n"
          :
          : "r"(R_U_re), "r"(R_U_im), "r"(HH_H_U_re), "r"(HH_H_U_im));
        if (t1 != t2) {
          R_L_re = &g_R.re[off_HH_H_L];
          R_L_im = &g_R.im[off_HH_H_L];
          __asm__ volatile(
            /* Lower triangle */
            "vle32.v v2, (%0)\n"
            "vle32.v v3, (%1)\n"
            "vfadd.vv v2, v2, v0\n"
            "vfsub.vv v3, v3, v1\n"
            "vse32.v v2, (%2)\n"
            "vse32.v v3, (%3)\n"
            :
            : "r"(R_L_re), "r"(R_L_im), "r"(HH_H_L_re), "r"(HH_H_L_im)
          );
        }

        sz -= vl;
        off_sc += vl;
      }
    }
}

/** Complex Cholesky decomposition L of a Hermitian positive-definite matrix HH_H
 *
 * HH_H = L*L^H
 *
 * \global g_HH_H matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \global g_L output lower triangular matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 */
void ccholesky_TxTx()
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
        g_L.re[off_ijk] = g_L.im[off_ijk] = 0;
        ++off_ijk;
      }
    }

  /* Calculate diagonal and lower triangular entries */
  for (i = 0; i < NUM_TX_ANT; ++i)
    for (j = 0; j < i+1; ++j)
      for (k = 0; k < NUM_SC; ++k){
        off_ijk = i * NUM_TX_ANT * NUM_SC + j * NUM_SC + k;
        sum_re = g_HH_H.re[off_ijk];
        sum_im = g_HH_H.im[off_ijk];
        for (l = 0; l < j; ++l) {
          off_ilk = i * NUM_TX_ANT * NUM_SC + l * NUM_SC + k;
          off_jlk = j * NUM_TX_ANT * NUM_SC + l * NUM_SC + k;
          t1_re = g_L.re[off_ilk];
          t1_im = g_L.im[off_ilk];
          t2_re = g_L.re[off_jlk];
          t2_im = g_L.im[off_jlk];
          sum_re -= t1_re * t2_re + t1_im * t2_im;
          sum_im -= t1_im * t2_re - t1_re * t2_im;
        }
        if (i == j)
          csqrt(
            sum_re, sum_im,
            &g_L.re[off_ijk],
            &g_L.im[off_ijk]);
        else {
          off_jjk = j * NUM_TX_ANT * NUM_SC + j * NUM_SC + k;
          cdiv(
            sum_re, sum_im,
            g_L.re[off_jjk],
            g_L.im[off_jjk],
            &g_L.re[off_ijk],
            &g_L.im[off_ijk]);
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

/** Complex matrix-vector multiplication HHy = HH*y
 *
 * \global g_HH matrix. Shape [NUM_TX_ANT][NUM_RX_ANT][NUM_SC]
 * \global g_y vector. Shape [NUM_RX_ANT][NUM_SC]
 * \global g_HHy output vector. Shape [NUM_TX_ANT][NUM_SC]
 */
void cmatvecmul_TxRx()
{
  size_t i, j;
  size_t off_HH, off_y, off_HHy, off_sc;
  size_t sz, vl;

  for (i = 0; i < NUM_TX_ANT; ++i) {
    off_HHy = i * NUM_SC;
    off_sc = 0;
    sz = NUM_SC;

    /* Initialize result registers */
    /* v0 - HHy real part */
    /* v1 - HHy imaginary part */
    __asm__ volatile(
      "vsetvli %0, %1, e32, m1, ta, ma\n"
      "vmv.v.i v0, 0\n"
      "vmv.v.i v1, 0\n"
      : "=r"(vl)
      : "r"(sz)
    );

    while (sz > 0) {
      for (j = 0; j < NUM_RX_ANT; ++j) {
        off_HH = i * NUM_RX_ANT * NUM_SC + j * NUM_SC + off_sc;
        off_y = j * NUM_SC + off_sc;
        __asm__ volatile(
          "vle32.v v2, (%0)\n"
          "vle32.v v3, (%1)\n"
          "vle32.v v4, (%2)\n"
          "vle32.v v5, (%3)\n"
          /* real part */
          "vfmacc.vv v0, v2, v4\n"
          "vfnmsac.vv v0, v3, v5\n"
          /* imaginary part */
          "vfmacc.vv v1, v3, v4\n"
          "vfmacc.vv v1, v2, v5\n"
          :
          : "r"(&g_HH.re[off_HH]), "r"(&g_HH.im[off_HH]),
            "r"(&g_y.re[off_y]), "r"(&g_y.im[off_y])
          );
      }

      /* Store result */
      __asm__ volatile(
        "vse32.v v0, (%0)\n"
        "vse32.v v1, (%1)\n"
        :
        : "r"(&g_HHy.re[off_HHy]), "r"(&g_HHy.im[off_HHy])
      );

      sz -= vl;
      off_HHy += vl;
      off_sc += vl;
    }
  }
}

/** Complex forward substitution L*z = HHy
 * 
 * z_i = (HHy_i - \sum_{k=0}^{i-1} L_{ik} z_k) / L_{ii}
 * 
 * \global g_L lower triangular matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \global g_HHy vector. Shape [NUM_TX_ANT][NUM_SC]
 * \global g_z output vector. Shape [NUM_TX_ANT][NUM_SC]
 */
void cforwardsub_TxTx()
{
  size_t i, j;
  size_t sz, vl;
  size_t off_sc;

  for (i = 0; i != NUM_TX_ANT; ++i) {
    off_sc = 0;
    sz = NUM_SC;

    while (sz > 0) {
      // printf("sz: %lu\n", sz);
      // printf("vl: %lu\n", vl);
      /* Initialize result registers as b */
      /* v0 - result real part */
      /* v1 - result imaginary part */
      __asm__ volatile (
        "vsetvli %0, %1, e32, m1, ta, ma\n"
        : "=r"(vl)
        : "r"(sz));
      __asm__ volatile (
        "vle32.v v0, (%0)\n"
        "vle32.v v1, (%1)\n"
        :
        : "r"(&g_HHy.re[i * NUM_SC + off_sc]),
          "r"(&g_HHy.im[i * NUM_SC + off_sc]));
      // printf("vl: %lu\n", vl);

      for (j = 0; j != i; ++j) {
        /* b - sum L_ij * z_j */
        /* v2 - L real part */
        /* v3 - L imaginary part */
        /* v4 - result_j real part */
        /* v5 - result_j imaginary part */
        __asm__ volatile (
          "vle32.v v2, (%0)\n"
          "vle32.v v3, (%1)\n"
          "vle32.v v4, (%2)\n"
          "vle32.v v5, (%3)\n"
          /* real part */
          "vfnmsac.vv v0, v2, v4\n"
          "vfmacc.vv v0, v3, v5\n"
          /* imaginary part */
          "vfnmsac.vv v1, v3, v4\n"
          "vfnmsac.vv v1, v2, v5\n"
          :
          : "r"(&g_L.re[i * NUM_TX_ANT * NUM_SC + j * NUM_SC + off_sc]),
            "r"(&g_L.im[i * NUM_TX_ANT * NUM_SC + j * NUM_SC + off_sc]),
            "r"(&g_z.re[j * NUM_SC + off_sc]),
            "r"(&g_z.im[j * NUM_SC + off_sc])
        );
      }

      /* Divide by L_ii */
      /* v2 - L_ii real part */
      /* v3 - L_ii imaginary part */
      __asm__ volatile (
        "vle32.v v2, (%0)\n"
        "vle32.v v3, (%1)\n"
        /* calculate L_ii_re^2 + L_ii_im^2 -> v4 */
        "vfmul.vv v4, v2, v2\n"
        "vfmacc.vv v4, v3, v3\n"
        /* real part */
        "vfmul.vv v5, v0, v2\n"
        "vfmacc.vv v5, v1, v3\n"
        "vdiv.vv v0, v5, v4\n"
        /* imaginary part */
        "vfmul.vv v6, v1, v2\n"
        "vfnmsac.vv v6, v0, v3\n"
        "vfdiv.vv v1, v6, v4\n"
        /* store result */
        "vse32.v v0, (%2)\n"
        "vse32.v v1, (%3)\n"
        :
        : "r"(&g_L.re[i * NUM_TX_ANT * NUM_SC + i * NUM_SC + off_sc]),
          "r"(&g_L.im[i * NUM_TX_ANT * NUM_SC + i * NUM_SC + off_sc]),
          "r"(&g_z.re[i * NUM_SC + off_sc]),
          "r"(&g_z.im[i * NUM_SC + off_sc])
      );

      off_sc += vl;
      sz -= vl;
    }
  }
}

/** Complex backward substitution L^H*x_MMSE = z
 * 
 * x_i = (z_i - \sum_{j=i+1}^{n-1} L_{ji} x_k) / L_{ii}
 * 
 * \global g_L lower triangular matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \global g_z rhs vector. Shape [NUM_TX_ANT][NUM_SC]
 * \global g_x_MMSE output vector. Shape [NUM_TX_ANT][NUM_SC]
 */
void cbackwardsub_TxTx()
{
  size_t i, j;
  size_t sz, vl;
  size_t off_sc;

  for (i = NUM_TX_ANT - 1; i != 0; --i) {
    off_sc = 0;
    sz = NUM_SC;

    while (sz > 0){
      /* Initialize result registers as z */
      /* v0 - result real part */
      /* v1 - result imaginary part */
      __asm__ volatile(
        "vsetvli %0, %1, e32, m1, ta, ma\n"
        : "=r"(vl)
        : "r"(sz));
      __asm__ volatile(
        "vle32.v v0, (%0)\n"
        "vle32.v v1, (%1)\n"
        :
        : "r"(&g_z.re[i * NUM_SC + off_sc]),
          "r"(&g_z.im[i * NUM_SC + off_sc]));

      for (j = i + 1; j < NUM_TX_ANT; ++j) {
        /* b - sum L_ji * z_j */
        /* v2 - L real part */
        /* v3 - L imaginary part */
        /* v4 - x_MMSE_j real part */
        /* v5 - x_MMSE_j imaginary part */
        __asm__ volatile(
          "vle32.v v2, (%0)\n"
          "vle32.v v3, (%1)\n"
          "vle32.v v4, (%2)\n"
          "vle32.v v5, (%3)\n"
          /* real part */
          "vfnmsac.vv v0, v2, v4\n"
          "vfmacc.vv v0, v3, v5\n"
          /* imaginary part */
          "vfnmsac.vv v1, v3, v4\n"
          "vfnmsac.vv v1, v2, v5\n"
          :
          : "r"(&g_L.re[j * NUM_TX_ANT * NUM_SC + i * NUM_SC + off_sc]),
            "r"(&g_L.im[j * NUM_TX_ANT * NUM_SC + i * NUM_SC + off_sc]),
            "r"(&g_x_MMSE.re[j * NUM_SC + off_sc]),
            "r"(&g_x_MMSE.im[j * NUM_SC + off_sc])
        );
      }

      /* Divide by L_ii */
      /* v2 - L_ii real part */
      /* v3 - L_ii imaginary part */
      __asm__ volatile (
        "vle32.v v2, (%0)\n"
        "vle32.v v3, (%1)\n"
        /* calculate L_ii_re^2 + L_ii_im^2 -> v4 */
        "vfmul.vv v4, v2, v2\n"
        "vfmacc.vv v4, v3, v3\n"
        /* real part */
        "vfmul.vv v5, v0, v2\n"
        "vfmacc.vv v5, v1, v3\n"
        "vdiv.vv v0, v5, v4\n"
        /* imaginary part */
        "vfmul.vv v6, v1, v2\n"
        "vfnmsac.vv v6, v0, v3\n"
        "vfdiv.vv v1, v6, v4\n"
        /* store HH_H */
        "vse32.v v0, (%2)\n"
        "vse32.v v1, (%3)\n"
      :
      : "r"(&g_L.re[i * NUM_TX_ANT * NUM_SC + i * NUM_SC + off_sc]),
        "r"(&g_L.im[i * NUM_TX_ANT * NUM_SC + i * NUM_SC + off_sc]),
        "r"(&g_x_MMSE.re[i * NUM_SC + off_sc]),
        "r"(&g_x_MMSE.im[i * NUM_SC + off_sc])
      );

      sz -= vl;
      off_sc += vl;
    }
  }
}

/*
 * MMSE
 */

void mmse()
{
  /* H^H*H + R */
  TIME(
    "Gram and add (RxTx x TxRx + TxTx): %ld\n",
    cmatgram_TxRx_cadd);
  /* L: (H^H*H + R) = L*L^H */
  TIME(
    "Cholesky (TxTx): %ld\n",
    ccholesky_TxTx);
  /* z: L*z = H^H*y */
  TIME(
    "Matrix-vector multiplication (TxRx x Rx): %ld\n",
    cmatvecmul_TxRx);
  TIME(
    "Forward substitution (TxTx): %ld\n",
    cforwardsub_TxTx);
  /* x_MMSE: L^H*x_MMSE = z */
  TIME(
    "Backward substitution (TxTx): %ld\n",
    cbackwardsub_TxTx);
}

acc_t mse()
{
  acc_t sum = 0.;
  size_t off_ij = 0;
  data_t sub1, sub2;
  for (; off_ij != NUM_TX_ANT * NUM_SC; ++off_ij) {
    sub1 = g_x.re[off_ij] - g_x_MMSE.re[off_ij];
    sub2 = g_x.im[off_ij] - g_x_MMSE.im[off_ij];
    sum += (sub1 * sub1 + sub2 * sub2) / (data_t)(NUM_TX_ANT * NUM_SC);
  }
  return sum;
}
