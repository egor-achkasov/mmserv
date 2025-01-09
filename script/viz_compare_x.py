#!/usr/bin/env python3

from util import read_defines

import numpy as np
from matplotlib import pyplot as plt

# read `out/x_mmse.bin` and plot the complex signal data samples

NUM_RX_ANT, NUM_TX_ANT, NUM_SC = read_defines()
x = np.fromfile("data/x_re.bin", dtype=np.float32) \
    + 1j * np.fromfile("data/x_im.bin", dtype=np.float32)
x = x.reshape((NUM_TX_ANT, NUM_SC))
x_mmse = np.fromfile("out/x_mmse_re.bin", dtype=np.float32) \
    + 1j * np.fromfile("out/x_mmse_im.bin", dtype=np.float32)
x_mmse = x_mmse.reshape((NUM_TX_ANT, NUM_SC))
x_mmse_python = np.fromfile("out/x_mmse_python_re.bin", dtype=np.float32) \
    + 1j * np.fromfile("out/x_mmse_python_im.bin", dtype=np.float32)
x_mmse_python = x_mmse_python.reshape((NUM_TX_ANT, NUM_SC))

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
    plt.scatter(x_mmse_python[tx].real,
                x_mmse_python[tx].imag,
                label=f"MMSE Python Tx {tx}",
                marker="^",
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
