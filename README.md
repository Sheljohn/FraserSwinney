Fraser-Swinney lag estimate
===========================

This is a Mex wrapper for the original script `minfo.c` written by Eric Weeks (see references).
The script estimates the mutual-information between an input time-series and its lagged versions for several lags.
This is typically used to select the optimal lag for phase-space reconstruction. The method is published in:

> Independent coordinates for strange attractors from mutual information<br>
> A. Fraser, H. Swinney, PRA 1986 (doi: 10.1103/PhysRevA.33.1134)

## Compiling

Simply run
```
    mex -largeArrayDims -O -lut milag.cpp
```
from the command-line in Matlab. Note the linking to `-lut` for keyboard-interruption support.

## Usage

The script takes one to three arguments:

- `X` -- an `Ntimes` x `Nchannels` matrix of time-series in each column
- `lag_max` (optional) -- the maximum lag to consider (see `DEFAULT_LAG_MAX`)
- `nbins` (optional) -- number of bins

The computation can be interrupted from the Matlab command-line using `Ctrl+C`.

## Notes

It is actually unclear to me how `nbins` affects the results, and I cannot see how this corresponds to
the number of bins (which assumes that mutual information is estimated from a histogram, which is really not immediate
from the code). I would advise not to use this input.

The original author (Eric Weeks) released his code in the public domain, but Dmytro Lituiev copyrighted his version despite a minimal contribution.
I'm releasing this in the public domain, please don't hesitate to drop me a line for copyright concerns.

## Output

Returns a `(lag_max+1)` x `Nchannels` matrix with the MI estimates for each lag and channel.

## References

1. [Andrew M. Fraser and Harry L. Swinney. "Independent coordinates for strange attractors from mutual information", Phys. Rev. A 33 (1986) 1134-1140.]( http://dx.doi.org/10.1103%2fPhysRevA.33.1134 )
2. [Eric Weeks' page](http://www.physics.emory.edu/~weeks/software/minfo.html) with the description of the modified algorithm

## Contributors

- Jonathan Hadida (2016), University of Oxford
- Dmytro S. Lituiev (2013), University of Zurich
- Eric Weeks C script (1997)
