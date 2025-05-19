#include "../include/common.h"

#include <stddef.h>

/** Complex Cholesky decomposition of a Hermitian positive-definite matrix G
 *
 * LL (floating point solution): 
 * G = L*L^H
 * L_ij = (G_ij - \sum_{k=0}^{j-1} L_ik L_jk^*) / L_jj
 * L_ii = sqrt(G_ii - \sum_{k=0}^{i-1} L_ik L_ik^*)
 * 
 * LDL (fixed point solution):
 * G = L*D*L^H
 * L_ij = (G_ij - \sum_{k=0}^{j-1} L_ik D_k L_jk^*) / D_j
 * D_i = G_ii - \sum_{k=0}^{i-1} L_ik D_k L_ik^*
 *
 * \global g_G matrix. Shape [NUM_TX][NUM_TX][NUM_SC]
 * \global g_L output lower triangular matrix. Shape [NUM_TX][NUM_TX][NUM_SC]
 * \global g_D output diagonal matrix (if DATA_TYPE_fixed is defined). Shape [NUM_TX][NUM_SC]
 */
void ccholesky()
{
#if defined(ARCH_x86) || defined(ARCH_rv)
  size_t i, j, k, s;
  size_t off_ij, off_jj, off_ii;
  size_t off_ik, off_jk;
  size_t off_i, off_j, off_k;
  data_t tmp; /* Temporary variable for sqrt */
  acc_t sum_re, sum_im;
  for (i = 0; i < NUM_TX; ++i) {
    for (j = 0; j <= i; ++j) {
      for (s = 0; s < NUM_SC; ++s) {
        off_ij = i * NUM_TX * NUM_SC + j * NUM_SC + s;
        sum_im = sum_re = 0;

        /* Calculate the sum */
        for (k = 0; k < j; ++k) {
          off_ik = i * NUM_TX * NUM_SC + k * NUM_SC + s;
          off_jk = j * NUM_TX * NUM_SC + k * NUM_SC + s;
#if defined(DATA_TYPE_float)
          sum_re += g_L.re[off_ik] * g_L.re[off_jk]
                  - g_L.im[off_ik] * g_L.im[off_jk];
#elif defined(DATA_TYPE_fixed)
          sum_re += (g_L.re[off_ik] * g_L.re[off_jk]
                  - g_L.im[off_ik] * g_L.im[off_jk])
                  * g_D[k * NUM_SC + s];
#else
#error "Unknown data type"
#endif
          sum_im += g_L.re[off_ik] * g_L.im[off_jk]
                  + g_L.im[off_ik] * g_L.re[off_jk];
        }

        if (i == j) {
          off_ii = i * NUM_TX * NUM_SC + i * NUM_SC + s;
#if defined(DATA_TYPE_float)

#if defined(ARCH_x86)
          __asm__ volatile (
            "flds %1\n"
            "fsubs %2\n"
            "fsqrt\n"
            "fstps %0\n"
            : "=m" (g_L.re[off_ii])
            : "m" (g_G.re[off_ij]), "m" (sum_re)
          );
#elif defined(ARCH_rv)
          __asm__ volatile (
            "fsub.s %0, %1, %2\n"   /* tmp = g_G.re[off_ij] - sum_re */
            "fsqrt.s %0, %0\n"      /* tmp = sqrtf(tmp) */
            : "=&f"(tmp) : "f"(g_G.re[off_ij]), "f"(sum_re)
          );
          g_L.re[off_ii] = tmp;
          g_L.im[off_ii] = 0;
#else
#error "Unknown architecture"
#endif

#elif defined(DATA_TYPE_fixed)
          /* Calculate D_i = G_ii - sum */
          g_D[i * NUM_SC + s] = g_G.re[off_ii] - (data_t)(sum_re >> FP_Q);
#else
#error "Unknown data type"
#endif
        } else { /* i != j */
#if defined(DATA_TYPE_float)
          /* Calculate L_ij = (G_ij - sum) / L_jj */
          off_jj = j * NUM_TX * NUM_SC + j * NUM_SC + s;
          g_L.re[off_ij] = (g_G.re[off_ij] - sum_re) / g_L.re[off_jj];
          g_L.im[off_ij] = (g_G.im[off_ij] - sum_im) / g_L.re[off_jj];
#elif defined(DATA_TYPE_fixed)
          /* Calculate L_ij = (G_ij - sum) / D_j */
          off_j = j * NUM_SC + s;
          /* real */
          sum_re = ((acc_t)g_G.re[off_ij] << FP_Q) - sum_re;
          /* TODO roubding? */
          g_L.re[off_ij] = (data_t)(sum_re / (acc_t)g_D[off_j]);
          /* imaginary */
          sum_im = ((acc_t)g_G.im[off_ij] << FP_Q) - sum_im;
          /* TODO roubding? */
          g_L.im[off_ij] = (data_t)(sum_im / (acc_t)g_D[off_j]);
#else
#error "Unknown data type"
#endif
        }
      }
    }
  }
#elif defined(ARCH_rvv)
  size_t i, j, k;
  size_t sz, vl;
  size_t off_sc;

  /* Init float registers */
  register float f0 __asm__("f0") = 2.0f;

  for (i = 0; i < NUM_TX; ++i)
    for (j = 0; j <= i; ++j) {
      sz = NUM_SC;
      off_sc = 0;

      while (sz > 0) {
        /* Initialize L registers */
        /* v0 - L_ij real part = G_ij.re */
        /* v1 - L_ij imaginary part  = G_ij.im */
        __asm__ volatile(
          "vsetvli %0, %1, e32, m1, ta, ma\n"
          : "=r"(vl) : "r"(sz));
        __asm__ volatile(
          "vle32.v v0, (%0)\n"
          "vle32.v v1, (%1)\n"
        :
        : "r"(&g_G.re[i * NUM_TX * NUM_SC + j * NUM_SC + off_sc]),
          "r"(&g_G.im[i * NUM_TX * NUM_SC + j * NUM_SC + off_sc])
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
          : "r"(&g_L.re[i * NUM_TX * NUM_SC + k * NUM_SC + off_sc]),
            "r"(&g_L.im[i * NUM_TX * NUM_SC + k * NUM_SC + off_sc]),
            "r"(&g_L.re[j * NUM_TX * NUM_SC + k * NUM_SC + off_sc]),
            "r"(&g_L.im[j * NUM_TX * NUM_SC + k * NUM_SC + off_sc])
          );
        }

        /* G_ii - sum */
        __asm__ volatile(
          "vfsub.vv v0, v0, v2\n"
          "vfsub.vv v1, v1, v3\n"
        );

        if (i == j) {
          /* Calculate L_ii = sqrt(G_ii - sum_{k=0}^{i-1} L_ik L_ik^*) */
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
          /* Calculate L_ij = (G_ij - sum) / L_jj */
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
          : "r"(&g_L.re[j * NUM_TX * NUM_SC + j * NUM_SC + off_sc]),
            "r"(&g_L.im[j * NUM_TX * NUM_SC + j * NUM_SC + off_sc])
          );
        }

        /* Store result */
        __asm__ volatile(
          "vse32.v v0, (%0)\n"
          "vse32.v v1, (%1)\n"
          :
          : "r"(&g_L.re[i * NUM_TX * NUM_SC + j * NUM_SC + off_sc]),
            "r"(&g_L.im[i * NUM_TX * NUM_SC + j * NUM_SC + off_sc])
        );

        sz -= vl;
        off_sc += vl;
      }
    }
#else
#error "Unknown architecture"
#endif
}
