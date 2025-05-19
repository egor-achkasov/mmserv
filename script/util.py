from os import path
import numpy as np


def load_xHRy(NUM_RX, NUM_TX, NUM_SC) -> tuple:
    """Load x, H, R and y from the data files.
    Assumes that the following files are present in the data directory:
    x_re.bin, x_im.bin, H_re.bin, H_im.bin, R_re.bin, R_im.bin, y_re.bin, y_im.bin
    and they contain float32 data of the appropriate sizes.

    Returns:
        tuple: x, H, R, y
    """
    x = np.fromfile("data/x_re.bin", dtype=np.float32)
    x = x + 1j*np.fromfile("data/x_im.bin", dtype=np.float32)
    x = x.reshape((NUM_TX, NUM_SC))
    H = np.fromfile("data/H_re.bin", dtype=np.float32)
    H = H + 1j*np.fromfile("data/H_im.bin", dtype=np.float32)
    H = H.reshape((NUM_RX, NUM_TX, NUM_SC))
    R = np.fromfile("data/R_re.bin", dtype=np.float32)
    R = R + 1j*np.fromfile("data/R_im.bin", dtype=np.float32)
    R = R.reshape((NUM_TX, NUM_TX, NUM_SC))
    y = np.fromfile("data/y_re.bin", dtype=np.float32)
    y = y + 1j*np.fromfile("data/y_im.bin", dtype=np.float32)
    y = y.reshape((NUM_RX, NUM_SC))

    return x, H, R, y
