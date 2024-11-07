#!/usr/bin/python

from util import read_defines, deinterleave, interleave

import numpy as np
from os import path, makedirs

NUM_RX_ANT, NUM_TX_ANT, NUM_SC = read_defines()

# Read the data from the files
x = deinterleave(np.fromfile("data/x.bin", dtype=np.float32), (NUM_TX_ANT, NUM_SC))
H = deinterleave(np.fromfile("data/H.bin", dtype=np.float32), (NUM_RX_ANT, NUM_TX_ANT, NUM_SC))
R = deinterleave(np.fromfile("data/R.bin", dtype=np.float32), (NUM_TX_ANT, NUM_TX_ANT, NUM_SC))
y = deinterleave(np.fromfile("data/y.bin", dtype=np.float32), (NUM_RX_ANT, NUM_SC))

x_mmse = np.empty((NUM_TX_ANT, NUM_SC), np.complex64)

for sc in range(NUM_SC):
    HH = H[..., sc].conj().T
    L = np.linalg.cholesky(HH @ H[..., sc] + R[..., sc])
    HHy = np.einsum("ij,j->i", HH, y[:, sc])
    print("\nHHy:\n", HHy)
    print("\nL:\n", L)

    z = np.empty((NUM_TX_ANT,), dtype=np.complex64)
    for i in range(NUM_TX_ANT):
        z[i] = (HHy[i] - np.sum([L[i, j]*z[j] for j in range(i)])) / L[i, i]
    print("\nz:\n", z)

    LH = L.conj().T
    for i in range(NUM_TX_ANT-1, -1, -1):
        x_mmse[i, sc] = (z[i] - np.sum([LH[i, j]*x_mmse[j, sc] for j in range(i+1, NUM_TX_ANT)], 0)) / LH[i, i]
    print("\nx_mmse:\n", x_mmse)

# Output the data
x_mmse_data = interleave(x_mmse)
out_dir = path.join(path.dirname(__file__), "..", "out")
if not path.exists(out_dir):
    makedirs(out_dir)
with open(path.join(out_dir, "x_mmse_python.bin"), "wb") as f:
    f.write(x_mmse_data.tobytes())
