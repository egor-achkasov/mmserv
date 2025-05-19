#!/usr/bin/python

import numpy as np
from numpy.random import random, normal
from sys import argv
from os import path, makedirs

if len(argv) != 4:
    print("Usage: python gen_data.py <NUM_RX> <NUM_TX> <NUM_SC>")
    exit(1)
NUM_RX = int(argv[1])  # number of receive antennas
NUM_TX = int(argv[2])  # number of transmit antennas
NUM_SC = int(argv[3])  # number of subcarriers

NOISE_STD_DEVIATION = np.sqrt(.5) / 100  # noise standard deviation

# Transmitter signal
x = random((NUM_TX, NUM_SC)) \
    + 1.j * random((NUM_TX, NUM_SC))
x = (x - 0.5 - 0.5j) * 2  # scale it from [0, 1] to [-1, 1]
# Channel
H = random((NUM_RX, NUM_TX, NUM_SC)) \
    + 1.j * random((NUM_RX, NUM_TX, NUM_SC))
H = (H - 0.5 - 0.5j) * 2  # scale it from [0, 1] to [-1, 1]
# Noise
n = normal(0, NOISE_STD_DEVIATION, (NUM_RX, NUM_SC)) \
    + 1.j * normal(0, NOISE_STD_DEVIATION, (NUM_RX, NUM_SC))
# Received signal
y = np.einsum("ijk,jk->ik", H, x) + n
# Noise covariance matrix
R = np.eye(NUM_TX, NUM_TX, dtype=np.complex64) * NOISE_STD_DEVIATION**2

data_tuple = (
    x.real.astype(np.float32), x.imag.astype(np.float32),
    H.real.astype(np.float32), H.imag.astype(np.float32),
    R.real.astype(np.float32), R.imag.astype(np.float32),
    y.real.astype(np.float32), y.imag.astype(np.float32)
)
data_filenames = (
    "data/x_re.bin", "data/x_im.bin",
    "data/H_re.bin", "data/H_im.bin",
    "data/R_re.bin", "data/R_im.bin",
    "data/y_re.bin", "data/y_im.bin"
)

# Create "data" directory if it does not exist
data_dir = path.join(path.dirname(__file__), "..", "data")
if not path.exists(data_dir):
    makedirs(data_dir)

# Write data to bin files
for data, filename in zip(data_tuple, data_filenames):
    with open(filename, "wb") as f:
        f.write(data.tobytes())
