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
#if defined(ARCH_x86) || defined(ARCH_rv)
  size_t t, r, s;
  size_t off_H, off_y, off_HHy;
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
  size_t i, j;
  size_t off_HH, off_y, off_HHy, off_sc;
  size_t sz, vl;

  for (i = 0; i < NUM_TX; ++i) {
    off_HHy = i * NUM_SC;
    off_sc = 0;
    sz = NUM_SC;

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

      for (j = 0; j < NUM_RX; ++j) {
        off_HH = i * NUM_RX * NUM_SC + j * NUM_SC + off_sc;
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
#else
#error "Unknown architecture"
#endif
}
