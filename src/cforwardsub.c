#include "../include/common.h"

#include <stddef.h>

/** Complex forward substitution L*z = HHy
 * 
 * z_t = (HHy_t - \sum_{tt=0}^{t-1} L_{t tt} z_tt) / L_{t t} (for LL / float solution)
 * z_t = (HHy_t - \sum_{tt=0}^{t-1} L_{t tt} z_tt) (for LDL / fixed solution)
 * 
 * \global g_L lower triangular matrix. Shape [NUM_TX][NUM_TX][NUM_SC]
 * \global g_HHy vector. Shape [NUM_TX][NUM_SC]
 * \global g_z output vector. Shape [NUM_TX][NUM_SC]
 */
void cforwardsub()
{
#if defined(ARCH_x86) || defined(ARCH_rv)
  size_t t, tt, s;
  size_t off_L, off_HHy, off_z;
  acc_t sum_re, sum_im;
  for (t = 0; t < NUM_TX; ++t) {
    for (s = 0; s < NUM_SC; ++s) {
      sum_re = sum_im = 0;
      for (tt = 0; tt < t; ++tt) {
        off_L = t * NUM_TX * NUM_SC + tt * NUM_SC + s;
        off_z = tt * NUM_SC + s;
        sum_re += (acc_t)g_L.re[off_L] * (acc_t)g_z.re[off_z]
                - (acc_t)g_L.im[off_L] * (acc_t)g_z.im[off_z];
        sum_im += (acc_t)g_L.re[off_L] * (acc_t)g_z.im[off_z]
                + (acc_t)g_L.im[off_L] * (acc_t)g_z.re[off_z];
      }
      off_HHy = t * NUM_SC + s;
      off_z = t * NUM_SC + s;
#if defined(DATA_TYPE_float)
      off_L = t * NUM_TX * NUM_SC + t * NUM_SC + s;
      g_z.re[off_z] = (g_HHy.re[off_HHy] - (data_t)sum_re) / g_L.re[off_L];
      g_z.im[off_z] = (g_HHy.im[off_HHy] - (data_t)sum_im) / g_L.re[off_L];
#elif defined(DATA_TYPE_fixed)
      g_z.re[off_z] = g_HHy.re[off_HHy] - (data_t)(sum_re >> FP_Q);
      g_z.im[off_z] = g_HHy.im[off_HHy] - (data_t)(sum_im >> FP_Q);
#else
#error "Unknown data type"
#endif
    }
  }
#elif defined(ARCH_rvv)
  size_t i, j;
  size_t sz, vl;
  size_t off_sc;

  for (i = 0; i < NUM_TX; ++i) {
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
          : "r"(&g_L.re[i * NUM_TX * NUM_SC + j * NUM_SC + off_sc]),
            "r"(&g_L.im[i * NUM_TX * NUM_SC + j * NUM_SC + off_sc]),
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
        : "r"(&g_L.re[i * NUM_TX * NUM_SC + i * NUM_SC + off_sc]),
          "r"(&g_L.im[i * NUM_TX * NUM_SC + i * NUM_SC + off_sc]),
          "r"(&g_z.re[i * NUM_SC + off_sc]),
          "r"(&g_z.im[i * NUM_SC + off_sc])
      );

      off_sc += vl;
      sz -= vl;
    }
  }
#else
#error "Unknown architecture"
#endif
}
