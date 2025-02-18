#include "include/define.h"

#include "../common/printf.h"

/* extern functions */
extern void cmatgram_TxRx_cadd();
extern void ccholesky_TxTx();
extern void cmatvecmul_TxRx();
extern void cforwardsub_TxTx();
extern void cbackwardsub_TxTx();
extern acc_t mse();

/* Raw data */
/* Transmitted signal */
extern data_t x_re[NUM_TX_ANT][NUM_SC];
extern data_t x_im[NUM_TX_ANT][NUM_SC];
/* Channel */
extern data_t H_re[NUM_RX_ANT][NUM_TX_ANT][NUM_SC];
extern data_t H_im[NUM_RX_ANT][NUM_TX_ANT][NUM_SC];
/* Noise covariance matrix */
extern data_t R_re[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
extern data_t R_im[NUM_TX_ANT][NUM_TX_ANT][NUM_SC];
/* Received signal */
extern data_t y_re[NUM_RX_ANT][NUM_SC];
extern data_t y_im[NUM_RX_ANT][NUM_SC];

/* Same data but casted to vcomplex */
vcomplex g_x = { .re = (data_t *)x_re, .im = (data_t *)x_im };
vcomplex g_H = { .re = (data_t *)H_re, .im = (data_t *)H_im };
vcomplex g_R = { .re = (data_t *)R_re, .im = (data_t *)R_im };
vcomplex g_y = { .re = (data_t *)y_re, .im = (data_t *)y_im };

/* MMSE approximation will be stored in this global*/
data_t x_MMSE_re[NUM_TX_ANT][NUM_SC];
data_t x_MMSE_im[NUM_TX_ANT][NUM_SC];
vcomplex g_x_MMSE = { .re = (data_t *)x_MMSE_re, .im = (data_t *)x_MMSE_im };

extern vcomplex g_HH, h_y, g_HHy;
size_t num_rx_cur=1, num_tx_cur=1, num_sc_cur=1;

int main() {
  unsigned long long start, end;

  size_t t1,t2,t3,t4,t5;
  for (num_rx_cur = 1; num_rx_cur <= NUM_RX_ANT; ++num_rx_cur)
    for (num_tx_cur = 1; num_tx_cur <= NUM_TX_ANT; ++num_tx_cur)
      for (num_sc_cur = 1; num_sc_cur <= NUM_SC; num_sc_cur += 1) {
        __asm__ volatile("rdcycle %0" : "=r"(start));
        cmatgram_TxRx_cadd();
        __asm__ volatile("rdcycle %0" : "=r"(end));
        t1 = end - start;
        __asm__ volatile("rdcycle %0" : "=r"(start));
        ccholesky_TxTx();
        __asm__ volatile("rdcycle %0" : "=r"(end));
        t2 = end - start;
        __asm__ volatile("rdcycle %0" : "=r"(start));
        cmatvecmul_TxRx();
        __asm__ volatile("rdcycle %0" : "=r"(end));
        t3 = end - start;
        __asm__ volatile("rdcycle %0" : "=r"(start));
        cforwardsub_TxTx();
        __asm__ volatile("rdcycle %0" : "=r"(end));
        t4 = end - start;
        __asm__ volatile("rdcycle %0" : "=r"(start));
        cbackwardsub_TxTx();
        __asm__ volatile("rdcycle %0" : "=r"(end));
        t5 = end - start;
        printf("%llu,%llu,%llu,%llu,%llu,", t1, t2, t3, t4, t5);
      }

  return 0;
}
