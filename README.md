# Metabolite inference

## Compiling the repo

### Using `xmake`

You can install `xmake` using the following command :

```bash
curl -fsSL https://xmake.io/shget.text | bash
```

Then clone the repo

```bash
git clone https://github.com/anticipated-chemistry-of-life/markov-random-field
cd markov-random-field
```

Finally:

```bash
xmake f -m release ; xmake
```

If you want to compile in debug mode :

```bash
xmake f -m debug ; xmake
```

### Using `cmake`

```bash
git clone https://github.com/anticipated-chemistry-of-life/markov-random-field
cd markov-random-field
```

#### Linux

```bash
micromamba create -n acol_env "cmake>=3.30" cxx-compiler zlib fmt armadillo "libopenblas=*=*openmp*" -c conda-forge
micromamba activate acol_env
export CC=$CONDA_PREFIX/bin/gcc
export CXX=$CONDA_PREFIX/bin/g++
cmake -S . -B build -DLOTUS=ON
```

#### MacOS

```bash
micromamba create -n acol_env "cmake>=3.30" cxx-compiler zlib fmt armadillo "libopenblas=*=*openmp*" -c conda-forge
micromamba activate acol_env
export CC=$CONDA_PREFIX/bin/clang
export CXX=$CONDA_PREFIX/bin/clang++
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
