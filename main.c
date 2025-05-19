#include "include/common.h"

#include <stddef.h> /* for size_t */


/*
 * Defines
 */

/* Got from https://elm-chan.org/junk/32bit/binclude.html */
/* Import a binary file */
#define IMPORT_BIN(sect, file, sym) __asm__ (\
    ".section " #sect "\n"                  /* Change section */\
    ".balign 4\n"                           /* Word alignment */\
    ".global " #sym "\n"                    /* Export the object address */\
    #sym ":\n"                              /* Define the object label */\
    ".incbin \"" file "\"\n"                /* Import the file */\
    ".global _sizeof_" #sym "\n"            /* Export the object size */\
    ".set _sizeof_" #sym ", . - " #sym "\n" /* Define the object size */\
    ".balign 4\n"                           /* Word alignment */\
    ".section \".text\"\n")                 /* Restore section */


/*
 * Global variables
 */

/* Import data from binary files */
IMPORT_BIN(.rodata, "data/x_re.bin", x_re);
IMPORT_BIN(.rodata, "data/x_im.bin", x_im);
IMPORT_BIN(.rodata, "data/H_re.bin", H_re);
IMPORT_BIN(.rodata, "data/H_im.bin", H_im);
IMPORT_BIN(.rodata, "data/R_re.bin", R_re);
IMPORT_BIN(.rodata, "data/R_im.bin", R_im);
IMPORT_BIN(.rodata, "data/y_re.bin", y_re);
IMPORT_BIN(.rodata, "data/y_im.bin", y_im);

/* Allocate space for MMSE raw data */
data_t G_re[NUM_TX][NUM_TX][NUM_SC];
data_t G_im[NUM_TX][NUM_TX][NUM_SC];
data_t L_re[NUM_TX][NUM_TX][NUM_SC];
data_t L_im[NUM_TX][NUM_TX][NUM_SC];
data_t g_D[NUM_TX][NUM_SC]; /* no imaginary part in D */
data_t HHy_re[NUM_TX][NUM_SC];
data_t HHy_im[NUM_TX][NUM_SC];
data_t z_re[NUM_TX][NUM_SC];
data_t z_im[NUM_TX][NUM_SC];
data_t x_MMSE_re[NUM_TX][NUM_SC];
data_t x_MMSE_im[NUM_TX][NUM_SC];

/* Initialize data */
vcomplex g_x = { .re = (data_t *)x_re, .im = (data_t *)x_im };
vcomplex g_H = { .re = (data_t *)H_re, .im = (data_t *)H_im };
vcomplex g_R = { .re = (data_t *)R_re, .im = (data_t *)R_im };
vcomplex g_y = { .re = (data_t *)y_re, .im = (data_t *)y_im };
vcomplex g_G = { .re = (data_t *)G_re, .im = (data_t *)G_im };
vcomplex g_L = { .re = (data_t *)L_re, .im = (data_t *)L_im };
vcomplex g_HHy = { .re = (data_t *)HHy_re, .im = (data_t *)HHy_im };
vcomplex g_z = { .re = (data_t *)z_re, .im = (data_t *)z_im };
vcomplex g_x_MMSE = { .re = (data_t *)x_MMSE_re, .im = (data_t *)x_MMSE_im };


/*
 * Read cycles macro and function
 */

/** Read and return the cycle counter value */
uint64_t readcycle() {
#if defined(ARCH_rvv) || defined(ARCH_rv)
  uint64_t val;
  __asm__ volatile("rdcycle %0" : "=r"(val));
  return val;
#elif defined(ARCH_x86)
  unsigned int hi, lo;
  __asm__ volatile("rdtsc" : "=a"(lo), "=d"(hi));
  return ((uint64_t)hi << 32) | lo;
#else
#error "Unknown architecture"
#endif
}

/** Read cycles for a function call and store the result in out */
#define FUNC_CYCLES(func, out) \
  do { \
    uint64_t start = readcycle(); \
    func; \
    uint64_t end = readcycle(); \
    out = end - start; \
  } while (0)


/*
 * Printf function
 */

#if defined(PLATFORM_ara)
#include "../../common/printf.h"
#elif defined(PLATFORM_linux)
#include <stdio.h>
#elif defined(PLATFORM_baremetal)
/* TODO */
#else
#error "Unknown platform"
#endif


/*
 * Main
 */

int main() {
  uint64_t start, end;
  uint64_t t1,t2,t3,t4,t5;

  FUNC_CYCLES(cmatgram(), t1);
  FUNC_CYCLES(ccholesky(), t2);
  FUNC_CYCLES(cmatvecmul(), t3);
  FUNC_CYCLES(cforwardsub(), t4);
  FUNC_CYCLES(cbackwardsub(), t5);
  printf("%lu,%lu,%lu,%lu,%lu,", t1, t2, t3, t4, t5);

  return 0;
}
