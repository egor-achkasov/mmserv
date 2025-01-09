#!/usr/bin/python

from util import read_defines

import numpy as np
from numpy.random import random, normal
from sys import argv
from os import path, makedirs

if "--help" in argv:
    print("Usage: python scripts/gen_data.py [--txt] [--bin] [--s]")
    print("  --txt: write the data/*.txt files")
    print("  --bin: write the data/*.bin files")
    print("  --s: print the .S file to stdout")
    exit(0)

# Change these flags to control the output
WRITE_DATA_TXT = "--txt" in argv
WRITE_DATA_BIN = "--bin" in argv
WRITE_DATA_S = "--s" in argv or len(argv) == 1

NUM_RX_ANT, NUM_TX_ANT, NUM_SC = read_defines()
NOISE_STD_DEVIATION = np.sqrt(.5) / 100  # noise standard deviation

# Transmitter signal
x = random((NUM_TX_ANT, NUM_SC)) \
    + 1.j * random((NUM_TX_ANT, NUM_SC))
x = (x - 0.5 - 0.5j) * 2  # scale it from [0, 1] to [-1, 1]
# Channel
H = random((NUM_RX_ANT, NUM_TX_ANT, NUM_SC)) \
    + 1.j * random((NUM_RX_ANT, NUM_TX_ANT, NUM_SC))
H = (H - 0.5 - 0.5j) * 2  # scale it from [0, 1] to [-1, 1]
# Noise
n = normal(0, NOISE_STD_DEVIATION, (NUM_RX_ANT, NUM_SC)) \
    + 1.j * normal(0, NOISE_STD_DEVIATION, (NUM_RX_ANT, NUM_SC))
# Received signal
y = np.einsum("ijk,jk->ik", H, x) + n
# Noise covariance matrix
R = np.eye(NUM_TX_ANT, NUM_TX_ANT, dtype=np.complex64) * NOISE_STD_DEVIATION**2

data_tuple = (
    x.real.astype(np.float32), x.imag.astype(np.float32),
    H.real.astype(np.float32), H.imag.astype(np.float32),
    R.real.astype(np.float32), R.imag.astype(np.float32),
    y.real.astype(np.float32), y.imag.astype(np.float32)
)


class Section:
    def __init__(self, name, source, align, sizeof, length, duplicate=False):
        self.name = name
        self.source = source
        self.align = align
        self.sizeof = sizeof
        self.length = length
        self.duplicate = duplicate


sections = [
    Section("x_re", "data/x_re.txt", "3", 32, x.size),
    Section("x_im", "data/x_im.txt", "3", 32, x.size),
    Section("H_re", "data/H_re.txt", "3", 32, H.size),
    Section("H_im", "data/H_im.txt", "3", 32, H.size),
    Section("R_re", "data/R_re.txt", "3", 32, R.size),
    Section("R_im", "data/R_im.txt", "3", 32, R.size),
    Section("y_re", "data/y_re.txt", "3", 32, y.size),
    Section("y_im", "data/y_im.txt", "3", 32, y.size),
]

# Create "data" directory if it does not exist
if WRITE_DATA_BIN or WRITE_DATA_TXT:
    data_dir = path.join(path.dirname(__file__), "..", "data")
    if not path.exists(data_dir):
        makedirs(data_dir)

if WRITE_DATA_TXT:
    for data, sec in zip(data_tuple, sections):
        with open(sec.source, "w") as f:
            for sample in data:
                f.write(f"{sample}\n")

if WRITE_DATA_BIN:
    for data, sec in zip(data_tuple, sections):
        with open(sec.source.replace(".txt", ".bin"), "wb") as f:
            f.write(data.tobytes())

if WRITE_DATA_S:
    print(".section .data,\"aw\",@progbits")
    for data, sec in zip(data_tuple, sections):
        print(f".global {sec.name}")
        print(f"{sec.name}:")
        for sample in data:
            bs = sample.tobytes()
            for i in range(0, len(bs), 4):
                s = ""
                for n in range(4):
                    s += "%02x" % bs[i+3-n]
                print("    .word 0x%s" % s)
