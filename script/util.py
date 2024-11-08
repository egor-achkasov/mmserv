from os import path
import numpy as np


def read_defines():
    """Read the defines from the define.h file

    Returns:
        int: Number of receive antennas
        int: Number of transmit antennas
        int: Number of subcarriers
    """
    with open(path.join(path.dirname(__file__), "..", "include", "define.h"), "r") as f:
        lines = f.read().split("\n")
    for line in lines:
        if line.startswith("#define NUM_RX_ANT "):
            NUM_RX_ANT = int(line[19:])
        if line.startswith("#define NUM_TX_ANT "):
            NUM_TX_ANT = int(line[19:])
        if line.startswith("#define NUM_SC "):
            NUM_SC = int(line[15:])

    # Assert that all the defines are read
    assert NUM_RX_ANT
    assert NUM_TX_ANT
    assert NUM_SC

    return NUM_RX_ANT, NUM_TX_ANT, NUM_SC


def interleave(data: np.ndarray) -> np.ndarray:
    """Cast a np.complex64 array to np.float32 array.

    Args:
        data (np.ndarray): The complex data to be casted. Dtype: np.complex64

    Returns:
        np.ndarray: The interleaved data. Dtype: np.float32
    """
    res = np.empty(2 * data.size, dtype=np.float32)
    res[0::2] = data.flatten().real.astype(np.float32)
    res[1::2] = data.flatten().imag.astype(np.float32)
    return res


def deinterleave(data: np.ndarray, shape: tuple) -> np.ndarray:
    """Cast a np.float32 array to np.complex64 array.

    Args:
        data (np.ndarray): The interleaved data to be casted. Dtype: np.float32
        shape (tuple): The shape of the deinterleaved data.

    Returns:
        np.ndarray: The deinterleaved data. Dtype: np.complex64
    """
    res = np.empty(shape, np.complex64)
    res.real = data[0::2].reshape(shape)
    res.imag = data[1::2].reshape(shape)
    return res
