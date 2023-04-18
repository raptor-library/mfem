# Fluids Solver

### Downloading and configuring spack (only needs to be done once)
Download spack from github and checkout stable release
```bash
git clone -c feature.manyFiles=true https://github.com/spack/spack.git
```

Add spack setup-env.sh to your .bashrc from the spack repo directory
```bash
echo -e "\nsource `pwd`/share/spack/setup-env.sh\n" >> ~/.bashrc
```

Source your ~/.bashrc
```bash
source ~/.bashrc
```

Configure spack to find compilers

```bash
spack compiler find
```

If you need load a feature currently in development on an MFEM branch edit the mfem spack package entry to change to custom branch
```bash
spack edit mfem
```
and change the contents to reflect the following
```bash
-    version('develop', branch='master')
+    version('develop', branch= NAMEOFBRANCH)
```

Use spack to install mfem and required dependencies
```bash
spack install mfem@develop~amgx~conduit~cuda~debug+examples~gnutls~gslib+lapack+libceed~libunwind+metis+miniapps~mpfr+mpi~netcdf~occa~openmp~petsc~pumi~raja~shared+static~strumpack+suite-sparse+sundials~superlu-dist~threadsafe~umpire+zlib
```

### Google Test installation with spack
Use spack to install googletest

```bash
spack install googletest
```

### Artistic Style (Code Beautifier) installation with spack
```bash
spack install astyle
```

### Installling from source
Create build directory
```bash
spack load mfem googletest astyle
mkdir build
cd build
```

Use CMake to build and install
```bash
cmake ../
make
make install
```
Run unit tests from build directory
```bash
./test/Test
```

### Code Beautifying
Run artistic style on source code to auto-format it.
```bash
astyle --options=./config/mfem.astylerc <source>
```
