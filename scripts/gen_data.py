#!/usr/bin/python

from util import read_defines, interleave

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
WRITE_DATA_S = "--s" in argv

NUM_RX_ANT, NUM_TX_ANT, NUM_SC = read_defines()
NOISE_STD_DEVIATION = np.sqrt(.5) / 10  # noise standard deviation

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

# Cast the complex data to int16
x_data = interleave(x)
H_data = interleave(H)
R_data = interleave(R)
y_data = interleave(y)
data_tuple = (x_data, H_data, R_data, y_data)


class Section:
    def __init__(self, name, source, align, sizeof, length, duplicate=False):
        self.name = name
        self.source = source
        self.align = align
        self.sizeof = sizeof
        self.length = length
        self.duplicate = duplicate


sections = [
    Section("x", "data/x.txt", "3", 16, x.size * 2),
    Section("H", "data/H.txt", "3", 16, H.size * 2),
    Section("R", "data/R.txt", "3", 16, R.size * 2),
    Section("y", "data/y.txt", "3", 16, y.size * 2),
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
            print(f"    .hword 0x{sample.view(np.uint16):04x} // {sample}")
