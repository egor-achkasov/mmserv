#include "../include/common.h"

#include <stddef.h>

/** Complex Gram matrix H^H*H and add complex matrix R
 * 
 * G = H^H*H + R
 * G_{t1t2} = \sum_{r=0}^{NUM_RX - 1} (H_{rt1}^* H_{rt2}) + R_{t1t2}
 * 
 * \global g_H matrix of channel coefficients. Shape [NUM_RX][NUM_TX][NUM_SC]
 * \global g_R noise covariance matrix. Shape [NUM_TX][NUM_TX][NUM_SC]
 * \global g_G output Gram matrix + R. Shape [NUM_TX][NUM_TX][NUM_SC]
 */
void cmatgram()
{
  size_t r, t1, t2, s;
  size_t off_G, off_GH, off_H1, off_H2;

#if defined(ARCH_x86) || defined(ARCH_rv)
  acc_t sum_re, sum_im;

  for (t1 = 0; t1 < NUM_TX; ++t1) {
    for (t2 = 0; t2 <= t1; ++t2) { /* t2 <= t1 since G is hermitian and R is symmetric */
      for (s = 0; s < NUM_SC; ++s) {

        /* Calculate the sum */
        sum_re = sum_im = 0;
        for (r = 0; r < NUM_RX; ++r) {
          off_H1 = r * NUM_TX * NUM_SC + t1 * NUM_SC + s;
          off_H2 = r * NUM_TX * NUM_SC + t2 * NUM_SC + s;
          sum_re += (acc_t)g_H.re[off_H1] * (acc_t)g_H.re[off_H2]
                  + (acc_t)g_H.im[off_H1] * (acc_t)g_H.im[off_H2];
          sum_im += (acc_t)g_H.re[off_H1] * (acc_t)g_H.im[off_H2]
                  - (acc_t)g_H.im[off_H1] * (acc_t)g_H.re[off_H2];
        }

        /* Add R */
        off_G = t1 * NUM_TX * NUM_SC + t2 * NUM_SC + s;
#if defined(DATA_TYPE_float)
        g_G.re[off_G] = (data_t)sum_re + g_R.re[off_G];
        g_G.im[off_G] = (data_t)sum_im + g_R.im[off_G];
#elif defined(DATA_TYPE_fixed)
        g_G.re[off_G] = (data_t)(sum_re >> FP_Q) + g_R.re[off_G];
        g_G.im[off_G] = (data_t)(sum_im >> FP_Q) + g_R.im[off_G];
#else
#error "Unknown data type"
#endif

        /* Fill the upper triangle */
        /* G_{t2t1} = G_{t1t2}^* */
        if (t1 != t2) {
          off_GH = t2 * NUM_TX * NUM_SC + t1 * NUM_SC + s;
          g_G.re[off_GH] = g_G.re[off_G];
          g_G.im[off_GH] = -g_G.im[off_G];
        }
      } 
    }
  }

#elif defined(ARCH_rvv)
  size_t sz, vl;

  for (t1 = 0; t1 != NUM_TX; ++t1)
    for (t2 = t1; t2 != t1; ++t2) {
      sz = NUM_SC;
      s = 0;

      while (sz > 0) {
        /* Initialize G registers */
        /* v0 - G real part */
        /* v4 - G imaginary part */
        __asm__ volatile(
          "vsetvli %0, %1, e32, m4, ta, ma\n"
          "vmv.v.i v0, 0\n"
          "vmv.v.i v4, 0\n"
          : "=r"(vl) : "r"(sz)
        );

        for (r = 0; r != NUM_RX; ++r) {
          off_H1 = r * NUM_TX * NUM_SC + t1 * NUM_SC + s;
          off_H2 = r * NUM_TX * NUM_SC + t2 * NUM_SC + s;

          /* Calculate H^H*H */
          /* v8  - H_{rt1} real part */
          /* v12 - H_{rt2} imaginary part */
          /* v16 - H_{rt2} real part */
          /* v20 - H_{rt2} imaginary part */
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
            /* real part */
            "vsmul.vv v24, v8, v16\n"
            "vsadd.vv v0, v0, v24\n"
            "vsmul.vv v24, v12, v20\n"
            "vsadd.vv v0, v0, v24\n"
            /* imaginary part */
            "vsmul.vv v24, v12, v16\n"
            "vsadd.vv v4, v4, v24\n"
            "vsmul.vv v24, v8, v20\n"
            "vssub.vv v4, v4, v24\n"
#else
#error "Unknown data type"
#endif
            :
            : "r"(g_H.re[off_H1]), "r"(g_H.im[off_H1]),
              "r"(g_H.re[off_H2]), "r"(g_H.im[off_H2])
          );
        }

        /* Add R */
        off_G = t1 * NUM_TX * NUM_SC + t2 * NUM_SC + s;
        __asm__ volatile(
          "vle32.v v8, (%0)\n"
          "vle32.v v12, (%1)\n"
#if defined(DATA_TYPE_float)
          "vfadd.vv v0, v0, v8\n"
          "vfadd.vv v4, v4, v12\n"
#elif defined(DATA_TYPE_fixed)
          "vsadd.vv v0, v0, v8\n"
          "vsadd.vv v4, v4, v12\n"
#else
#error "Unknown data type"
#endif
          "vse32.v v0, (%2)\n"
          "vse32.v v4, (%3)\n"
          :
          : "r"(&g_R.re[off_G]), "r"(&g_R.im[off_G]),
            "r"(&g_G.re[off_G]), "r"(&g_G.im[off_G])
        );

        /* Fill the upper triangle */
        /* G_{t2t1} = G_{t1t2}^* */
        if (t1 != t2) {
          off_GH = t2 * NUM_TX * NUM_SC + t1 * NUM_SC + s;
          __asm__ volatile(
#if defined(DATA_TYPE_float)
            "vfneg.v v4, v4\n"
#elif defined(DATA_TYPE_fixed)
            "vneg.v v4, v4\n"
#else
#error "Unknown data type"
#endif
            "vse32.v v0, (%0)\n"
            "vse32.v v4, (%1)\n"
            :
            : "r"(&g_G.re[off_GH]), "r"(&g_G.im[off_GH])
          );
        }

        sz -= vl;
        s += vl;
      }
    }
#endif
}
