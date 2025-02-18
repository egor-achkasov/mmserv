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
extern size_t num_rx_cur, num_tx_cur, num_sc_cur;


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

  for (t1 = 0; t1 != num_tx_cur; ++t1)
    for (t2 = t1; t2 != num_tx_cur; ++t2) {
      off_sc = 0;
      off_HH_H_L = t1 * num_tx_cur * num_sc_cur + t2 * num_sc_cur;
      off_HH_H_U = t2 * num_tx_cur * num_sc_cur + t1 * num_sc_cur;
      HH_H_L_re = &g_HH_H.re[off_HH_H_L];
      HH_H_L_im = &g_HH_H.im[off_HH_H_L];
      HH_H_U_re = &g_HH_H.re[off_HH_H_U];
      HH_H_U_im = &g_HH_H.im[off_HH_H_U];
      sz = num_sc_cur;

      while (sz > 0) {
        /* Initialize HH_H registers */
        /* v0 - HH_H real part */
        /* v1 - HH_H imaginary part */
        __asm__ volatile(
          "vsetvli %0, %1, e32, m1, ta, ma\n"
          "vmv.v.i v0, 0\n"
          "vmv.v.i v1, 0\n"
          : "=r"(vl) : "r"(sz));

        for (r = 0; r != num_rx_cur; ++r) {
          off_A = r * num_tx_cur * num_sc_cur + t1 * num_sc_cur + off_sc;
          off_AH = r * num_tx_cur * num_sc_cur + t2 * num_sc_cur + off_sc;
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
 * L_ij = (HH_H_ij - \sum_{k=0}^{j-1} L_ik L_jk^*) / L_jj
 * L_ii = sqrt(HH_H_ii - \sum_{k=0}^{i-1} L_ik L_ik^*)
 *
 * \global g_HH_H matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 * \global g_L output lower triangular matrix. Shape [NUM_TX_ANT][NUM_TX_ANT][NUM_SC]
 */
void ccholesky_TxTx()
{
  size_t i, j, k;
  size_t sz, vl;
  size_t off_sc;

  /* Init float registers */
  register float f0 __asm__("f0") = 2.0f;

  for (i = 0; i < num_tx_cur; ++i)
    for (j = 0; j <= i; ++j) {
      sz = num_sc_cur;
      off_sc = 0;

      while (sz > 0) {
        /* Initialize L registers */
        /* v0 - L_ij real part = HH_H_ij.re */
        /* v1 - L_ij imaginary part  = HH_H_ij.im */
        __asm__ volatile(
          "vsetvli %0, %1, e32, m1, ta, ma\n"
          : "=r"(vl) : "r"(sz));
        __asm__ volatile(
          "vle32.v v0, (%0)\n"
          "vle32.v v1, (%1)\n"
        :
        : "r"(&g_HH_H.re[i * num_tx_cur * num_sc_cur + j * num_sc_cur + off_sc]),
          "r"(&g_HH_H.im[i * num_tx_cur * num_sc_cur + j * num_sc_cur + off_sc])
        );

        /* Calculate sum_{k=0}^{j-1} L_ik L_jk^* */
        /* v2 - sum real part */
        /* v3 - sum imaginary part */
        __asm__ volatile(
          "vmv.v.i v2, 0\n"
          "vmv.v.i v3, 0\n"
        );
        for (k = 0; k < j; ++k) {
          __asm__ volatile(
            "vle32.v v4, (%0)\n"
            "vle32.v v5, (%1)\n"
            "vle32.v v6, (%2)\n"
            "vle32.v v7, (%3)\n"
            /* real part */
            "vfmacc.vv v2, v4, v6\n"
            "vfmacc.vv v2, v5, v7\n"
            /* imaginary part */
            "vfmacc.vv v3, v5, v6\n"
            "vfnmsac.vv v3, v4, v7\n"
          :
          : "r"(&g_L.re[i * num_tx_cur * num_sc_cur + k * num_sc_cur + off_sc]),
            "r"(&g_L.im[i * num_tx_cur * num_sc_cur + k * num_sc_cur + off_sc]),
            "r"(&g_L.re[j * num_tx_cur * num_sc_cur + k * num_sc_cur + off_sc]),
            "r"(&g_L.im[j * num_tx_cur * num_sc_cur + k * num_sc_cur + off_sc])
          );
        }

        /* HH_H_ii - sum */
        __asm__ volatile(
          "vfsub.vv v0, v0, v2\n"
          "vfsub.vv v1, v1, v3\n"
        );

        if (i == j) {
          /* Calculate L_ii = sqrt(HH_H_ii - sum_{k=0}^{i-1} L_ik L_ik^*) */
          __asm__ volatile(
            /* Complex sqrt */

            /* v2 = r = sqrt(re^2 + im^2) */
            "vfmul.vv v2, v0, v0\n"
            "vfmacc.vv v2, v1, v1\n"
            "vfsqrt.v v2, v2\n"

            /* v3 - real part */
            "vfadd.vv v3, v2, v0\n" /* r + re */
            "vfdiv.vf v3, v3, f0\n" /* (r + re) / 2 */
            "vfsqrt.v v3, v3\n" /* sqrt((r + re) / 2) */
            /* v4 - imaginary part */
            "vfsub.vv v4, v2, v0\n" /* r - re */
            "vfdiv.vf v4, v4, f0\n" /* (r - re) / 2 */
            "vfsqrt.v v4, v4\n" /* sqrt((r - re) / 2) */
            "vfsgnj.vv v4, v4, v1\n" /* sgn(im) * sqrt((r - re) / 2) */

            /* TODO handle im == 0 */

            /* Move the result to v0 and v1 */
            "vmv.v.v v0, v3\n"
            "vmv.v.v v1, v4\n"
          );
        } else {
          /* Calculate L_ij = (HH_H_ij - sum) / L_jj */
          __asm__ volatile(
            /* L_jj */
            "vle32.v v2, (%0)\n"
            "vle32.v v3, (%1)\n"
            /* calculate L_jj_re^2 + L_jj_im^2 -> v4 */
            "vfmul.vv v4, v2, v2\n"
            "vfmacc.vv v4, v3, v3\n"
            /* real part */
            "vfmul.vv v5, v0, v2\n"
            "vfmacc.vv v5, v1, v3\n"
            /* imaginary part */
            "vfmul.vv v6, v1, v2\n"
            "vfnmsac.vv v6, v0, v3\n"
            /* divide and store at v0 and v1 */
            "vfdiv.vv v0, v5, v4\n"
            "vfdiv.vv v1, v6, v4\n"
          :
          : "r"(&g_L.re[j * num_tx_cur * num_sc_cur + j * num_sc_cur + off_sc]),
            "r"(&g_L.im[j * num_tx_cur * num_sc_cur + j * num_sc_cur + off_sc])
          );
        }

        /* Store result */
        __asm__ volatile(
          "vse32.v v0, (%0)\n"
          "vse32.v v1, (%1)\n"
          :
          : "r"(&g_L.re[i * num_tx_cur * num_sc_cur + j * num_sc_cur + off_sc]),
            "r"(&g_L.im[i * num_tx_cur * num_sc_cur + j * num_sc_cur + off_sc])
        );

        sz -= vl;
        off_sc += vl;
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

  for (i = 0; i < num_tx_cur; ++i) {
    off_HHy = i * num_sc_cur;
    off_sc = 0;
    sz = num_sc_cur;

    while (sz > 0) {
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

      for (j = 0; j < num_rx_cur; ++j) {
        off_HH = i * num_rx_cur * num_sc_cur + j * num_sc_cur + off_sc;
        off_y = j * num_sc_cur + off_sc;
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

  for (i = 0; i < num_tx_cur; ++i) {
    off_sc = 0;
    sz = num_sc_cur;

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
        : "r"(&g_HHy.re[i * num_sc_cur + off_sc]),
          "r"(&g_HHy.im[i * num_sc_cur + off_sc]));
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
          : "r"(&g_L.re[i * num_tx_cur * num_sc_cur + j * num_sc_cur + off_sc]),
            "r"(&g_L.im[i * num_tx_cur * num_sc_cur + j * num_sc_cur + off_sc]),
            "r"(&g_z.re[j * num_sc_cur + off_sc]),
            "r"(&g_z.im[j * num_sc_cur + off_sc])
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
        : "r"(&g_L.re[i * num_tx_cur * num_sc_cur + i * num_sc_cur + off_sc]),
          "r"(&g_L.im[i * num_tx_cur * num_sc_cur + i * num_sc_cur + off_sc]),
          "r"(&g_z.re[i * num_sc_cur + off_sc]),
          "r"(&g_z.im[i * num_sc_cur + off_sc])
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

  for (i = num_tx_cur - 1; i != (size_t)-1; --i) {
    off_sc = 0;
    sz = num_sc_cur;

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
        : "r"(&g_z.re[i * num_sc_cur + off_sc]),
          "r"(&g_z.im[i * num_sc_cur + off_sc]));

      for (j = i + 1; j < num_tx_cur; ++j) {
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
          : "r"(&g_L.re[j * num_tx_cur * num_sc_cur + i * num_sc_cur + off_sc]),
            "r"(&g_L.im[j * num_tx_cur * num_sc_cur + i * num_sc_cur + off_sc]),
            "r"(&g_x_MMSE.re[j * num_sc_cur + off_sc]),
            "r"(&g_x_MMSE.im[j * num_sc_cur + off_sc])
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
      : "r"(&g_L.re[i * num_tx_cur * num_sc_cur + i * num_sc_cur + off_sc]),
        "r"(&g_L.im[i * num_tx_cur * num_sc_cur + i * num_sc_cur + off_sc]),
        "r"(&g_x_MMSE.re[i * num_sc_cur + off_sc]),
        "r"(&g_x_MMSE.im[i * num_sc_cur + off_sc])
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
  size_t off = 0;
  data_t sub1, sub2;
  size_t sz = num_tx_cur * num_sc_cur, vl;
  register data_t num_tx_num_sc_reg __asm__("f0") = (data_t)(num_tx_cur * NUM_SC); 

  while (sz > 0) {
    __asm__ volatile (
      "vsetvli %0, %1, e32, m1, ta, ma\n"
      : "=r"(vl)
      : "r"(sz));
    __asm__ volatile (
      "vle32.v v0, (%1)\n"
      "vle32.v v1, (%2)\n"
      "vle32.v v2, (%3)\n"
      "vle32.v v3, (%4)\n"
      "vfsub.vv v0, v0, v2\n"
      "vfsub.vv v1, v1, v3\n"
      "vfmul.vv v0, v0, v0\n"
      "vfmul.vv v1, v1, v1\n"
      "vfadd.vv v0, v0, v1\n"
      "vfdiv.vf v0, v0, %5\n"
      "vfredusum.vs v0, v4, v5\n"
      "vmv.x.s %0, v0\n"
      : "+r"(sum)
      : "r"(&g_x.re[off]),
        "r"(&g_x.re[off]),
        "r"(&g_x_MMSE.re[off]),
        "r"(&g_x_MMSE.re[off]),
        "f"(num_tx_num_sc_reg));
    sz -= vl;
    off += vl;
  }

  return sum;
}
