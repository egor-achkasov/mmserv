#include "../include/common.h"

#include <stddef.h>

/** Complex backward substitution L^H*x_MMSE = z
 * 
 * x_{MMSE}_t = (z_t - \sum_{tt=t+1}^{NUM_TX-1} L_{tt t} x_{MMSE}_{tt}) / L_{t t} (for LL / float solution)
 * x_{MMSE}_t = (z_t / D_t - \sum_{tt=t+1}^{NUM_TX-1} L_{tt t} x_{MMSE}_{tt}) (for LDL / fixed solution)
 * 
 * \global g_L lower triangular matrix. Shape [NUM_TX][NUM_TX][NUM_SC]
 * \global g_D diagonal matrix (only if DATA_TYPE_fixed is defined). Shape [NUM_TX][NUM_SC]
 * \global g_z rhs vector. Shape [NUM_TX][NUM_SC]
 * \global g_x_MMSE output vector. Shape [NUM_TX][NUM_SC]
 */
void cbackwardsub()
{
  size_t t, tt, s;
  size_t off_L, off_z, off_x_MMSE, off_D;

#if defined(ARCH_x86) || defined(ARCH_rv)
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
  size_t sz, vl;

  for (t = NUM_TX - 1; t != (size_t)-1; --t) {
    sz = NUM_SC;
    s = 0;

    while (sz > 0){
      /* Initialize x_MMSE as z */
      /* v0 - x_MMSE real part */
      /* v4 - x_MMSE imaginary part */
      off_z = t * NUM_SC + s;
      __asm__ volatile(
        "vsetvli %0, %1, e32, m4, ta, ma\n"
        "vle32.v v0, (%2)\n"
        "vle32.v v4, (%3)\n"
        : "=r"(vl)
        : "r"(sz),
          "r"(&g_z.re[off_z]), "r"(&g_z.im[off_z])
      );

#if defined(DATA_TYPE_fixed)
      /* Divide by D_t */
      __asm__ volatile (
        "vle32.v v8, (%0)\n"
        "vfdiv.vv v0, v0, v8\n"
        "vfdiv.vv v4, v4, v8\n"
        :
        : "r"(&g_D[t * NUM_SC + s])
      );
#endif

      for (tt = tt + 1; tt < NUM_TX; ++tt) {
        /* z - sum L_{tt t} * x_{MMSE}_{tt} or
         * z / D_t - sum L_{tt t} * x_{MMSE}_{tt} */
        off_L = tt * NUM_TX * NUM_SC + t * NUM_SC + s;
        off_x_MMSE = tt * NUM_SC + s;
        __asm__ volatile(
          "vle32.v v8, (%0)\n"
          "vle32.v v12, (%1)\n"
          "vle32.v v16, (%2)\n"
          "vle32.v v20, (%3)\n"
#if defined(DATA_TYPE_float)
          /* real part */
          "vfnmsac.vv v0, v8, v16\n"
          "vfmacc.vv v0, v12, v20\n"
          /* imaginary part */
          "vfnmsac.vv v4, v12, v16\n"
          "vfnmsac.vv v4, v8, v20\n"
#elif defined(DATA_TYPE_fixed)
          /* real part */
          "vsmul.vv v24, v8, v16\n"
          "vssub.vv v0, v0, v24\n"
          "vsmul.vv v24, v12, v20\n"
          "vsadd.vv v0, v0, v24\n"
          /* imaginary part */
          "vsmul.vv v24, v12, v16\n"
          "vssub.vv v4, v4, v24\n"
          "vsmul.vv v24, v8, v20\n"
          "vssub.vv v4, v4, v24\n"
#else
#error "Unknown data type"
#endif
          :
          : "r"(&g_L.re[off_L]), "r"(&g_L.im[off_L]),
            "r"(&g_x_MMSE.re[off_x_MMSE]), "r"(&g_x_MMSE.im[off_x_MMSE])
        );
      }

#if defined(DATA_TYPE_float)
      /* Divide by L_ii */
      off_L = t * NUM_TX * NUM_SC + t * NUM_SC + s;
      __asm__ volatile (
        "vle32.v v8, (%0)\n"
        "vfdiv.vv v0, v0, v8\n"
        "vfdiv.vv v4, v4, v8\n"
        :
        : "r"(&g_L.re[off_L])
      );
#endif

      /* Store result */
      off_x_MMSE = t * NUM_SC + s;
      __asm__ volatile(
        "vse32.v v0, (%0)\n"
        "vse32.v v4, (%1)\n"
        :
        : "r"(&g_x_MMSE.re[off_x_MMSE]), "r"(&g_x_MMSE.im[off_x_MMSE])
      );

      sz -= vl;
      s += vl;
    }
  }

#else
#error "Unknown architecture"
#endif
}
