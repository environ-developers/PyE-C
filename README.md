# PyE-C - Environ Python wrapping tools
    `PyE-C` turns Environ into a Python engine for integration with other multiscale methods

## Installation
Set the Environ installation directory so that PyE-C can latch onto these libraries
 + envdir=${ENVIRON}
 + pip install `PyE-C`

## FAQ
 - Some Intel MPI/MKL errors occur. What do I do?
	+ Try `export LD_PRELOAD=/opt/intel/mkl/lib/intel64/libmkl_rt.so`
