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
  size_t i, j, k, s;
  size_t off_ij, off_jj, off_ii;
  size_t off_ik, off_jk;
  size_t off_i, off_j, off_k;

#if defined(ARCH_x86) || defined(ARCH_rv)
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
          sum_re = ((acc_t)g_G.re[off_ij] << FP_Q) - sum_re;
          g_L.re[off_ij] = (data_t)(sum_re / (acc_t)g_D[off_j]);
          sum_im = ((acc_t)g_G.im[off_ij] << FP_Q) - sum_im;
          g_L.im[off_ij] = (data_t)(sum_im / (acc_t)g_D[off_j]);
#else
#error "Unknown data type"
#endif
        }
      }
    }
  }
#elif defined(ARCH_rvv)
  size_t sz, vl;

  for (i = 0; i < NUM_TX; ++i)
    for (j = 0; j <= i; ++j) {
      sz = NUM_SC;
      s = 0;

      while (sz > 0) {
        /* Calculate
         * sum_{k=0}^{j-1} L_ik L_jk^* or
         * sum_{k=0}^{j-1} L_ik D_k L_jk^*
         * 
         * v0 - sum real part
         * v4 - sum imaginary part */
        __asm__ volatile(
          "vsetvli %0, %1, e32, m4, ta, ma\n"
          "vmv.v.i v0, 0\n"
          "vmv.v.i v4, 0\n"
          : "=r"(vl) : "r"(sz)
        );

        for (k = 0; k < j; ++k) {
          off_ik = i * NUM_TX * NUM_SC + k * NUM_SC + s;
          off_jk = j * NUM_TX * NUM_SC + k * NUM_SC + s;
          __asm__ volatile(
            "vle32.v v8, (%0)\n"
            "vle32.v v12, (%1)\n"
            "vle32.v v16, (%2)\n"
            "vle32.v v20, (%3)\n"
#if defined(DATA_TYPE_float)
            /* real part */
            "vfmacc.vv v0, v8, v16\n"
            "vfmacc.vv v0, v12, v20\n"
            /* imaginary part */
            "vfmacc.vv v4, v12, v16\n"
            "vfnmsac.vv v4, v8, v20\n"
#elif defined(DATA_TYPE_fixed)
            "vle32.v v24, (%4)\n"
            /* real part */
            "vsmul.vv v28, v8, v24\n"
            "vsmul.vv v28, v28, v16\n"
            "vsadd.vv v0, v0, v28\n"
            "vsmul.vv v28, v12, v24\n"
            "vsmul.vv v28, v28, v20\n"
            "vsadd.vv v0, v0, v28\n"
            /* imaginary part */
            "vsmul.vv v28, v12, v16\n"
            "vsmul.vv v28, v28, v20\n"
            "vsadd.vv v4, v4, v28\n"
            "vsmul.vv v28, v8, v24\n"
            "vsmul.vv v28, v28, v20\n"
            "vssub.vv v4, v4, v28\n"
#else
#error "Unknown data type"
#endif
          :
          : "r"(&g_L.re[off_ik]), "r"(&g_L.im[off_ik]),
            "r"(&g_L.re[off_jk]), "r"(&g_L.im[off_jk])
#if defined(DATA_TYPE_fixed)
            , "r"(&g_D[k * NUM_SC + s])
#endif
          );
        }

        if (i == j) {
          off_ii = i * NUM_TX * NUM_SC + i * NUM_SC + s;
          /* L_ii = sqrt(G_ii - sum) or D_i = (G_ii - sum) */
          /* G_ii imaginary part is 0, so we can ignore it */
          __asm__ volatile(
            "vle32.v v8, (%0)\n"
#if defined(DATA_TYPE_float)
            "vfsub.vv v0, v8, v0\n"
            "vfsqrt.v v0, v0\n"
            "vfneg.v v4, v4\n"
#elif defined(DATA_TYPE_fixed)
            "vssub.vv v0, v8, v0\n"
            "vneg.v v4, v4\n"
#else
#error "Unknown data type"
#endif
            "vse32.v v0, (%1)\n"
            : : "r"(&g_G.re[off_ii]),
#if defined(DATA_TYPE_float)
            "r"(&g_L.re[off_ii])
#elif defined(DATA_TYPE_fixed)
            "r"(&g_D[i * NUM_SC + s])
#else
#error "Unknown data type"
#endif
          );
        } else { /* i != j */
          off_ij = i * NUM_TX * NUM_SC + j * NUM_SC + s;
          off_jj = j * NUM_TX * NUM_SC + j * NUM_SC + s;
          off_j = j * NUM_SC + s;
          /* Calculate L_ij = (G_ij - sum) / L_jj or L_ij = (G_ij - sum) / D_j */
          /* L_jj is always real, so we can ignore the imaginary part */
          __asm__ volatile(
            "vle32.v v8, (%0)\n"  /* G_ij.re */
            "vle32.v v12, (%1)\n" /* G_ij.im */
            "vle32.v v16, (%2)\n" /* L_jj or D_j */
#if defined(DATA_TYPE_float)
            "vfsub.vv v0, v8, v0\n"
            "vfsub.vv v4, v12, v4\n"
            "vfdiv.vv v0, v0, v16\n"
            "vfdiv.vv v4, v4, v16\n"
#elif defined(DATA_TYPE_fixed)
            "vssub.vv v0, v8, v0\n"
            "vssub.vv v4, v12, v4\n"
            "vdiv.vv v0, v0, v16\n"
            "vdiv.vv v4, v4, v16\n"
            "vsll.vi v0, v0, 1\n"
            "vsll.vi v4, v4, 1\n"
#else
#error "Unknown data type"
#endif
            /* Store result */
            "vse32.v v0, (%3)\n"
            "vse32.v v4, (%4)\n"
            :
            : "r"(&g_G.re[off_ij]), "r"(&g_G.im[off_ij]),
              "r"(&g_L.re[off_jj]),
              "r"(&g_L.re[off_ij]), "r"(&g_L.im[off_ij])
          );
        }

        sz -= vl;
        s += vl;
      }
    }
#else
#error "Unknown architecture"
#endif
}
