# Metabolite inference

## Compiling the repo

```bash
git clone https://github.com/anticipated-chemistry-of-life/markov-random-field
cd markov-random-field
```

```bash
micromamba create -n acol_env "cmake>=3.30" compilers zlib fmt armadillo "libopenblas=*=*openmp*" -c conda-forge
micromamba activate acol_env
export CC=$CONDA_PREFIX/bin/gcc
export CXX=$CONDA_PREFIX/bin/g++
export FC=$CONDA_PREFIX/bin/gfortran
cmake -S . -B build -DLOTUS=ON
```

Or on Mac :

```bash
micromamba create -n acol_env "cmake>=3.30" compilers zlib fmt armadillo "libopenblas=*=*openmp*" -c conda-forge
micromamba activate acol_env
export CC=$CONDA_PREFIX/bin/clang
export CXX=$CONDA_PREFIX/bin/clang++
export FC=$CONDA_PREFIX/bin/gfortran
cmake -S . -B build -DLOTUS=ON
```

Once this is done, you can compile the files using :

```bash
cmake --build build -j8
```

## Running the test suite

In the build directory:

```bash
cd build
make -j8 acol_unitTests
./acol_unitTests
```
