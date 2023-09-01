# Molecular Dynamics Code for IMTEK-HPC-MD class

## Compiling in the command line

The command line (terminal) may look daunting at first, but it has the advantage
of being the same across all UNIX platforms, and does not depend on a specific
IDE. The standard CMake workflow is to create a `build/` directory which will
contain all the build files. To do that, and compile your project, run:

```bash
cd <your repository>

# Create build directory
mkdir build
cd build

# Configure & compile
cmake -DCMAKE_BUILD_TYPE=Debug ..
make

# Run executable and tests
./milestones/01/01
make test
```

If there are no errors then you are all set! Note that the flag
`-DCMAKE_BUILD_TYPE=Debug` should be changed to
`-DCMAKE_BUILD_TYPE=Release` when you run a production simulation, i.e. a
simulation with more than a few hundred atoms. This turns on aggressive compiler
optimizations, which results in speedup. However, when writing the code and
looking for bugs, `Debug` should be used instead.

Try compiling and running tests with both compilation configurations.
