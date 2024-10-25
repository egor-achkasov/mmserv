#!/usr/bin/env python3

from util import read_defines, deinterleave

import numpy as np
from matplotlib import pyplot as plt

# read `out/x_mmse.bin` and plot the complex signal data samples

NUM_RX_ANT, NUM_TX_ANT, NUM_SC = read_defines()

with open("out/x_mmse.bin", "rb") as f:
    x_mmse = deinterleave(np.fromfile(f, dtype=np.int16), (NUM_TX_ANT, NUM_SC))
with open("data/x.bin", "rb") as f:
    x = deinterleave(np.fromfile(f, dtype=np.int16), (NUM_TX_ANT, NUM_SC))

for tx in range(NUM_TX_ANT):
    plt.scatter(x_mmse[tx].real,
                x_mmse[tx].imag,
                label=f"MMSE Tx {tx}",
                marker="x",
                color=plt.cm.hsv(tx / NUM_TX_ANT))
    plt.scatter(x[tx].real,
                x[tx].imag,
                label=f"Orig Tx {tx}",
                marker="o",
                color=plt.cm.hsv(tx / NUM_TX_ANT))
plt.axhline(0, color='black')
plt.axvline(0, color='black')
plt.xlim(-1.1, 1.1)
plt.ylim(-1.1, 1.1)
plt.title("Approximated MMSE Signal vs Original Signal Samples")
plt.ylabel("Imaginary")
plt.xlabel("Real")
plt.legend()
plt.show()
