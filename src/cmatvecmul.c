#include "../include/common.h"

#include <stddef.h>

/** Complex matrix-vector multiplication HHy = HH*y
 * 
 * HHy_t = \sum_{r=0}^{NUM_RX-1} H_{rt}^* * y_r
 *
 * \global g_H matrix. Shape [NUM_RX][NUM_TX][NUM_SC]
 * \global g_y vector. Shape [NUM_RX][NUM_SC]
 * \global g_HHy output vector. Shape [NUM_TX][NUM_SC]
 */
void cmatvecmul()
{
  size_t t, r, s;
  size_t off_H, off_y, off_HHy;

#if defined(ARCH_x86) || defined(ARCH_rv)
  acc_t sum_re, sum_im;

  for (t = 0; t < NUM_TX; ++t) {
    for (s = 0; s < NUM_SC; ++s) {
      sum_re = sum_im = 0;
      for (r = 0; r < NUM_RX; ++r) {
        off_H = r * NUM_TX * NUM_SC + t * NUM_SC + s;
        off_y = r * NUM_SC + s;
        sum_re += (acc_t)g_H.re[off_H] * (acc_t)g_y.re[off_y]
                - (acc_t)g_H.im[off_H] * (acc_t)g_y.im[off_y];
        sum_im += (acc_t)g_H.re[off_H] * (acc_t)g_y.im[off_y]
                + (acc_t)g_H.im[off_H] * (acc_t)g_y.re[off_y];
      }
      off_HHy = t * NUM_SC + s;
#if defined(DATA_TYPE_float)
      g_HHy.re[off_HHy] = (data_t)sum_re;
      g_HHy.im[off_HHy] = (data_t)sum_im;
#elif defined(DATA_TYPE_fixed)
      g_HHy.re[off_HHy] = (data_t)(sum_re >> FP_Q);
      g_HHy.im[off_HHy] = (data_t)(sum_im >> FP_Q);
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
      /* Initialize HHy with 0 */
      /* v0 - HHy real part */
      /* v4 - HHy imaginary part */
      __asm__ volatile(
        "vsetvli %0, %1, e32, m4, ta, ma\n"
        "vmv.v.i v0, 0\n"
        "vmv.v.i v4, 0\n"
        : "=r"(vl)
        : "r"(sz)
      );

      for (r = 0; r < NUM_RX; ++r) {
        off_H = r * NUM_TX * NUM_SC + t * NUM_SC + s;
        off_y = r * NUM_SC + s;
        __asm__ volatile(
          "vle32.v v8, (%0)\n"
          "vle32.v v12, (%1)\n"
          "vle32.v v16, (%2)\n"
          "vle32.v v20, (%3)\n"
#if defined(DATA_TYPE_float)
          /* real part */
          "vfmacc.vv v0, v8, v16\n"
          "vfnmsac.vv v0, v12, v20\n"
          /* imaginary part */
          "vfmacc.vv v4, v12, v16\n"
          "vfmacc.vv v4, v8, v20\n"
#elif defined(DATA_TYPE_fixed)
          /* real part */
          "vsmul.vv v24, v8, v16\n"
          "vsadd.vv v0, v0, v24\n"
          "vsmul.vv v24, v12, v20\n"
          "vssub.vv v0, v0, v24\n"
          /* imaginary part */
          "vsmul.vv v24, v12, v16\n"
          "vsadd.vv v4, v4, v24\n"
          "vsmul.vv v24, v8, v20\n"
          "vsadd.vv v4, v4, v24\n"
#else
#error "Unknown data type"
#endif
          :
          : "r"(&g_H.re[off_H]), "r"(&g_H.im[off_H]),
            "r"(&g_y.re[off_y]), "r"(&g_y.im[off_y])
        );
      }

      /* Store result */
      off_HHy = t * NUM_SC + s;
      __asm__ volatile(
        "vse32.v v0, (%0)\n"
        "vse32.v v4, (%1)\n"
        :
        : "r"(&g_HHy.re[off_HHy]), "r"(&g_HHy.im[off_HHy])
      );

      sz -= vl;
      s += vl;
    }
  }
#else
#error "Unknown architecture"
#endif
}
