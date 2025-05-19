#include "../include/common.h"

#include <stddef.h>

/** Complex backward substitution L^H*x_MMSE = z
 * 
 * x_MMSE_t = (z_t - \sum_{tt=t+1}^{NUM_TX-1} L_{tt t} x_tt) / L_{t t} (for LL / float solution)
 * x_MMSE_t = (z_t / D_t - \sum_{tt=t+1}^{NUM_TX-1} L_{tt t} x_tt) (for LDL / fixed solution)
 * 
 * \global g_L lower triangular matrix. Shape [NUM_TX][NUM_TX][NUM_SC]
 * \global g_D diagonal matrix (only if DATA_TYPE_fixed is defined). Shape [NUM_TX][NUM_SC]
 * \global g_z rhs vector. Shape [NUM_TX][NUM_SC]
 * \global g_x_MMSE output vector. Shape [NUM_TX][NUM_SC]
 */
void cbackwardsub()
{
#if defined(ARCH_x86) || defined(ARCH_rv)
  size_t t, tt, s;
  size_t off_L, off_z, off_x_MMSE, off_D;
  acc_t sum_re, sum_im;
  for (t = NUM_TX - 1; t != (size_t)-1; --t) {
    for (s = 0; s < NUM_SC; ++s) {
      sum_re = sum_im = 0;
      for (tt = t + 1; tt < NUM_TX; ++tt) {
        off_L = tt * NUM_TX * NUM_SC + t * NUM_SC + s;
        off_x_MMSE = tt * NUM_SC + s;
        sum_re += (acc_t)g_L.re[off_L] * (acc_t)g_x_MMSE.re[off_x_MMSE]
                - (acc_t)g_L.im[off_L] * (acc_t)g_x_MMSE.im[off_x_MMSE];
        sum_im += (acc_t)g_L.re[off_L] * (acc_t)g_x_MMSE.im[off_x_MMSE]
                + (acc_t)g_L.im[off_L] * (acc_t)g_x_MMSE.re[off_x_MMSE];
      }
      off_z = t * NUM_SC + s;
      off_x_MMSE = t * NUM_SC + s;
#if defined(DATA_TYPE_float)
      off_L = t * NUM_TX * NUM_SC + t * NUM_SC + s;
      g_x_MMSE.re[off_x_MMSE] = (g_z.re[off_z] - (data_t)sum_re) / g_L.re[off_L];
      g_x_MMSE.im[off_x_MMSE] = (g_z.im[off_z] - (data_t)sum_im) / g_L.re[off_L];
#elif defined(DATA_TYPE_fixed)
      off_D = t * NUM_SC + s;
      g_x_MMSE.re[off_x_MMSE] = (data_t)((acc_t)(g_z.re[off_z] << FP_Q) / g_D[off_D])
                              - (data_t)(sum_re >> FP_Q);
      g_x_MMSE.im[off_x_MMSE] = (data_t)((acc_t)(g_z.im[off_z] << FP_Q) / g_D[off_D])
                              - (data_t)(sum_im >> FP_Q);
#else
#error "Unknown data type"
#endif
    }
  }
#elif defined(ARCH_rvv)
  size_t i, j;
  size_t sz, vl;
  size_t off_sc;

  for (i = NUM_TX - 1; i != (size_t)-1; --i) {
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

      for (j = i + 1; j < NUM_TX; ++j) {
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
          : "r"(&g_L.re[j * NUM_TX * NUM_SC + i * NUM_SC + off_sc]),
            "r"(&g_L.im[j * NUM_TX * NUM_SC + i * NUM_SC + off_sc]),
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
      : "r"(&g_L.re[i * NUM_TX * NUM_SC + i * NUM_SC + off_sc]),
        "r"(&g_L.im[i * NUM_TX * NUM_SC + i * NUM_SC + off_sc]),
        "r"(&g_x_MMSE.re[i * NUM_SC + off_sc]),
        "r"(&g_x_MMSE.im[i * NUM_SC + off_sc])
      );

      sz -= vl;
      off_sc += vl;
    }
  }
#else
#error "Unknown architecture"
#endif
}
