<!-- # halos
 C++ code used in https://arxiv.org/abs/2504.20043 to compute the UV luminosity function fits and in https://arxiv.org/abs/2601.06023 to compute the weak lensing magnifications and Hubble diagram reconstruction. -->

# halos

C++ code for the stochastic modelling of dark matter halos, with two
applications: weak lensing magnification distributions and Hubble diagram
reconstruction ([arXiv:2601.06023](https://arxiv.org/abs/2601.06023)), and
UV luminosity function fits in CDM and beyond-CDM cosmologies
([arXiv:2504.20043](https://arxiv.org/abs/2504.20043)).

## Build

Requires a C++17 compiler and [GSL](https://www.gnu.org/software/gsl/).

## Lensing model settings

The defaults in `lensing.h` are the extended model of
[arXiv:XXXX.XXXXX](https://arxiv.org/abs/XXXX.XXXXX),[add later] which adds subhalo
substructure, the correlated clustering field with its conditional weak arm,
filament bias, and the real-space top-hat σ8 anchor to the model of
[arXiv:2601.06023](https://arxiv.org/abs/2601.06023). These form a coupled set,
so overriding only some of them throws. The commented block in
`main_lensing.cpp` recovers the original physics as a whole.

## Emulator

Running the Monte Carlo for every point of a cosmological parameter scan is
expensive. **[FLUMEN](https://github.com/baltabaygal/flumen)** is a conditional
normalizing-flow emulator of the magnification PDF trained on this model,
returning the distribution for a source redshift and six cosmological
parameters in milliseconds.

## License

MIT