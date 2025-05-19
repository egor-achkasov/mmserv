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
#if defined(ARCH_x86) || defined(ARCH_rv)
  size_t r, t1, t2, s;
  size_t off_G, off_R, off_H1, off_H2;
  for (t1 = 0; t1 < NUM_TX; ++t1) {
    for (t2 = 0; t2 < NUM_TX; ++t2) {
      for (s = 0; s < NUM_SC; ++s) {
        off_R = off_G = t1 * NUM_TX * NUM_SC + t2 * NUM_SC + s;
        g_G.re[off_G] = g_R.re[off_R];
        g_G.im[off_G] = g_R.im[off_R];

        for (r = 0; r < NUM_RX; ++r) {
          off_H1 = r * NUM_TX * NUM_SC + t1 * NUM_SC + s;
          off_H2 = r * NUM_TX * NUM_SC + t2 * NUM_SC + s;
          g_G.re[off_G] += g_H.re[off_H1] * g_H.re[off_H2]
                         + g_H.im[off_H1] * g_H.im[off_H2];
          g_G.im[off_G] += g_H.re[off_H1] * g_H.im[off_H2]
                         - g_H.im[off_H1] * g_H.re[off_H2];
        }
      } 
    }
  }
#elif defined(ARCH_rvv)
  size_t t1, t2, r;
  size_t sz, vl;
  size_t off_sc, off_A, off_AH;
  size_t off_G_L, off_G_U;
  data_t *A_re, *A_im, *AH_re, *AH_im;
  data_t *R_L_re, *R_L_im, *R_U_re, *R_U_im;
  data_t *G_L_re, *G_L_im, *G_U_re, *G_U_im;

  for (t1 = 0; t1 != NUM_TX; ++t1)
    for (t2 = t1; t2 != NUM_TX; ++t2) {
      off_sc = 0;
      off_G_L = t1 * NUM_TX * NUM_SC + t2 * NUM_SC;
      off_G_U = t2 * NUM_TX * NUM_SC + t1 * NUM_SC;
      G_L_re = &g_G.re[off_G_L];
      G_L_im = &g_G.im[off_G_L];
      G_U_re = &g_G.re[off_G_U];
      G_U_im = &g_G.im[off_G_U];
      sz = NUM_SC;

      while (sz > 0) {
        /* Initialize G registers */
        /* v0 - G real part */
        /* v1 - G imaginary part */
        __asm__ volatile(
          "vsetvli %0, %1, e32, m1, ta, ma\n"
          "vmv.v.i v0, 0\n"
          "vmv.v.i v1, 0\n"
          : "=r"(vl) : "r"(sz));

        for (r = 0; r != NUM_RX; ++r) {
          off_A = r * NUM_TX * NUM_SC + t1 * NUM_SC + off_sc;
          off_AH = r * NUM_TX * NUM_SC + t2 * NUM_SC + off_sc;
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
        R_U_re = &g_R.re[off_G_U];
        R_U_im = &g_R.im[off_G_U];
        __asm__ volatile(
          "vle32.v v2, (%0)\n"
          "vle32.v v3, (%1)\n"
          "vfadd.vv v2, v2, v0\n"
          "vfadd.vv v3, v3, v1\n"
          "vse32.v v2, (%2)\n"
          "vse32.v v3, (%3)\n"
          :
          : "r"(R_U_re), "r"(R_U_im), "r"(G_U_re), "r"(G_U_im));
        if (t1 != t2) {
          R_L_re = &g_R.re[off_G_L];
          R_L_im = &g_R.im[off_G_L];
          __asm__ volatile(
            /* Lower triangle */
            "vle32.v v2, (%0)\n"
            "vle32.v v3, (%1)\n"
            "vfadd.vv v2, v2, v0\n"
            "vfsub.vv v3, v3, v1\n"
            "vse32.v v2, (%2)\n"
            "vse32.v v3, (%3)\n"
            :
            : "r"(R_L_re), "r"(R_L_im), "r"(G_L_re), "r"(G_L_im)
          );
        }

        sz -= vl;
        off_sc += vl;
      }
    }
#endif
}
