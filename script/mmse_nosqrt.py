#!/usr/bin/python

from util import read_defines, load_xHRy

import numpy as np
#from scipy.linalg import ldl
from os import path, makedirs

# Read the data from the files
NUM_RX_ANT, NUM_TX_ANT, NUM_SC = read_defines()
x, H, R, y = load_xHRy(NUM_RX_ANT, NUM_TX_ANT, NUM_SC)
x_mmse = np.empty((NUM_TX_ANT, NUM_SC), np.complex64)


def ldl(A):
    n = A.shape[0]
    D = np.zeros_like(A)
    L = np.eye(n, dtype=A.dtype)
    for i in range(n):
        for j in range(i):
            s = 0.+0j
            for l in range(j):
                s += L[i][l] * L[j][l].conj() * D[l][l]
            L[i][j] = (A[i][j] - s) / D[j][j]
        s = 0.+0.j
        for j in range(i):
            s += L[i][j] * L[i][j].conj() * D[j][j]
        D[i][i] = A[i][i] - s
    return L, D


for sc in range(NUM_SC):
    HH = H[..., sc].conj().T
    L, D = ldl(HH @ H[..., sc] + R[..., sc])
    HHy = np.einsum("ij,j->i", HH, y[:, sc])

    z = np.empty((NUM_TX_ANT,), dtype=np.complex64)
    for i in range(NUM_TX_ANT):
        z[i] = (HHy[i] - np.sum([L[i, j]*z[j] for j in range(i)]))

    LH = L.conj().T
    for i in range(NUM_TX_ANT-1, -1, -1):
        x_mmse[i, sc] = z[i] / D[i, i] - np.sum([L[i, j] * x[j, sc] for j in range(i+1, NUM_TX_ANT)])

# Output the data
out_dir = path.join(path.dirname(__file__), "..", "out")
if not path.exists(out_dir):
    makedirs(out_dir)
with open(path.join(out_dir, "x_mmse_nosqrt_python_re.bin"), "wb") as f:
    f.write(x_mmse.real.astype(np.float32).tobytes())
with open(path.join(out_dir, "x_mmse_nosqrt_python_im.bin"), "wb") as f:
    f.write(x_mmse.imag.astype(np.float32).tobytes())
