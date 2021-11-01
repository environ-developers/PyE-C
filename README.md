# PyE-C - Environ Python wrapping tools
    `PyE-C` turns Environ into a Python engine for integration with other multiscale methods

## Installation
+ Compile Environ (see compilation instructions in Environ's README file)
+ Set the path to the Environ directory (if not `PyE-C/Environ`) with
  ```
  export envdir=absolute-path-to-Environ
  ```
+ Install PyE-C with
  ```
  pip install .
  ```

## FAQ
 - Some Intel MPI/MKL errors occur. What do I do?
	+ Try
      ```
      export LD_PRELOAD=/opt/intel/mkl/lib/intel64/libmkl_rt.so
      ```
    + exact path may vary
