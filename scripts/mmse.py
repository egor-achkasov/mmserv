#!/usr/bin/python

from util import read_defines, deinterleave, interleave

import numpy as np
from os import path, makedirs

NUM_RX_ANT, NUM_TX_ANT, NUM_SC = read_defines()

# Read the data from the files
x = deinterleave(np.fromfile("data/x.bin", dtype=np.int16), (NUM_TX_ANT, NUM_SC))
H = deinterleave(np.fromfile("data/H.bin", dtype=np.int16), (NUM_RX_ANT, NUM_TX_ANT, NUM_SC))
R = deinterleave(np.fromfile("data/R.bin", dtype=np.int16), (NUM_TX_ANT, NUM_TX_ANT, NUM_SC))
y = deinterleave(np.fromfile("data/y.bin", dtype=np.int16), (NUM_RX_ANT, NUM_SC))

# Move the SC dimension to the beginning
x = np.moveaxis(x, -1, 0)
H = np.moveaxis(H, -1, 0)
R = np.moveaxis(R, -1, 0)
y = np.moveaxis(y, -1, 0)

HH = H.conj().transpose(0, 2, 1)
L = np.linalg.cholesky(HH @ H + R)
HHy = np.einsum("...ij,...j->...i", HH, y)

z = np.empty((NUM_SC, NUM_TX_ANT), dtype=np.complex64)
for i in range(NUM_TX_ANT):
    z[:, i] = (HHy[:, i] - np.sum([L[:, i, j]*z[:, j] for j in range(i)], 0)) / L[:, i, i]

LH = L.conj().transpose(0, 2, 1)
x_mmse = np.empty((NUM_SC, NUM_TX_ANT), dtype=np.complex64)
for i in range(NUM_TX_ANT-1, -1, -1):
    x_mmse[:, i] = (z[:, i] - np.sum([LH[:, i, j]*x_mmse[:, j] for j in range(i+1, NUM_TX_ANT)], 0)) / LH[:, i, i]

# Move the SC dimension back to the end
x_mmse = np.moveaxis(x_mmse, 0, -1)
x = np.moveaxis(x, 0, -1)
H = np.moveaxis(H, 0, -1)
R = np.moveaxis(R, 0, -1)
y = np.moveaxis(y, 0, -1)

# Cast the complex data to int16
x_mmse_data = interleave(x_mmse)
out_dir = path.join(path.dirname(__file__), "..", "out")
if not path.exists(out_dir):
    makedirs(out_dir)
with open(path.join(out_dir, "x_mmse_python.bin"), "wb") as f:
    f.write(x_mmse_data.tobytes())
