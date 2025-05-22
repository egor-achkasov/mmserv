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
  size_t t, tt, s;
  size_t off_L, off_HHy, off_z;

#if defined(ARCH_x86) || defined(ARCH_rv)
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
  size_t sz, vl;

  for (t = 0; t < NUM_TX; ++t) {
    sz = NUM_SC;
    s = 0;

    while (sz > 0) {
      off_HHy = t * NUM_SC + s;

      /* Initialize z as HHy */
      /* v0 - z real part */
      /* v4 - z imaginary part */
      __asm__ volatile (
        "vsetvli %0, %1, e32, m4, ta, ma\n"
        "vle32.v v0, (%2)\n"
        "vle32.v v4, (%3)\n"
        : "=r"(vl)
        : "r"(sz),
          "r"(&g_HHy.re[off_HHy]), "r"(&g_HHy.im[off_HHy])
      );

      for (tt = 0; tt != t; ++tt) {
        off_L = t * NUM_TX * NUM_SC + tt * NUM_SC + s;
        off_z = tt * NUM_SC + s;

        /* HHy_t - sum L_{t tt} * z_{tt} */
        __asm__ volatile (
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
            "r"(&g_z.re[off_z]), "r"(&g_z.im[off_z])
        );
      }

#if defined(DATA_TYPE_float)
      /* Divide by L_ii */
      /* L_ii imaginary part is zero, so we ignore it */
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
      off_z = t * NUM_SC + s;
      __asm__ volatile (
        "vse32.v v0, (%0)\n"
        "vse32.v v4, (%1)\n"
        :
        : "r"(&g_z.re[off_z]), "r"(&g_z.im[off_z])
      );

      s += vl;
      sz -= vl;
    }
  }

#else
#error "Unknown architecture"
#endif
}
