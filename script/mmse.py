#!/usr/bin/python

""" Python implementation of the MMSE receiver with Cholesky decomposition """

from util import read_defines, load_xHRy

import numpy as np
from os import path, makedirs

# Read the data from the files
NUM_RX_ANT, NUM_TX_ANT, NUM_SC = read_defines()
x, H, R, y = load_xHRy(NUM_RX_ANT, NUM_TX_ANT, NUM_SC)
x_mmse = np.empty((NUM_TX_ANT, NUM_SC), np.complex64)

for sc in range(NUM_SC):
    HH = H[..., sc].conj().T
    L = np.linalg.cholesky(HH @ H[..., sc] + R[..., sc])
    HHy = np.einsum("ij,j->i", HH, y[:, sc])

    z = np.empty((NUM_TX_ANT,), dtype=np.complex64)
    for i in range(NUM_TX_ANT):
        z[i] = (HHy[i] - np.sum([L[i, j]*z[j] for j in range(i)])) / L[i, i]

    LH = L.conj().T
    for i in range(NUM_TX_ANT-1, -1, -1):
        x_mmse[i, sc] = (z[i] - np.sum([LH[i, j]*x_mmse[j, sc] for j in range(i+1, NUM_TX_ANT)], 0)) / LH[i, i]

# Output the data
out_dir = path.join(path.dirname(__file__), "..", "out")
if not path.exists(out_dir):
    makedirs(out_dir)
with open(path.join(out_dir, "x_mmse_python_re.bin"), "wb") as f:
    f.write(x_mmse.real.astype(np.float32).tobytes())
with open(path.join(out_dir, "x_mmse_python_im.bin"), "wb") as f:
    f.write(x_mmse.imag.astype(np.float32).tobytes())
