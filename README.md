# MMSE on RISC-V

Forked from [gfdmrv](https://gitlab.vodafone-chair.org/viktor.razilov/gfdmrv) and [lte-benchmark-code](https://gitlab.vodafone-chair.org/viktor.razilov/lte-benchmark-code).

## Configuration

Set up defines in `inc/define.h` before using the tool.
- NUM_RX_ANT - number of recieving antennas
- NUM_TX_ANT - number of transmitting antennas
- NUM_SC - number of subcarriers 

## Data generation

Generate data with `scripts/gen_data.py`.

```
$ python scripts/gen_data.py --help
```

There are 4 variables that are generated:
- x - transmitted signal
- H - channel signal
- R - noise correlation matrix
- y - received signal

The data is generated as complex floating point. It is interleaved into int16 before output. Interleaved data contains real values on even positions and imaginary values on odd positions. Interleaving and deinterleaving functions can be found in `scripts/util.py`.

- txt files are line separated entries
- bin files are densly concatenated values
- S file an asm data file with multiple sections for each variable. Note that running `python scripts/gen_data.py --s` outputs the file into stdout, so consider using `> data.S`.  

## Data visualization

The program outputs the approximated `x` into `out/x_mmse.bin`. Run `python scripts/viz_compare_x.py` to plot the actual and approximated signal samples, assuming that `data/x.bin` exists.

It is possible to compare the C implementation with a simple python numpy one. Running `python scripts/mmse.py` creates `out/x_mmse_python.bin` assuming that `data/x.bin` exists. Plot it together with the C approximation by `cp out/x_mmse_python.bin bin/x.bin`. Compare it with the original signal by `cp out/x_mmse_python.bin out/x_mmse.bin`.
