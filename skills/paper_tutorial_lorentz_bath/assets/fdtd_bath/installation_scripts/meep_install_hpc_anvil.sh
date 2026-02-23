#!/bin/bash
###############################################################################
#  meep_install_hpc_anvil.sh  –  Build MEEP + prerequisites in $HOME on Anvil
#  Author: Tao E. Li <taoeli@udel.edu>
###############################################################################
#  The Anvil HPC system is available to U.S. researchers through the Advanced 
#  Cyberinfrastructure Coordination Ecosystem: Services & Support (ACCESS)
#  program, which is supported by the National Science Foundation (NSF).
#  To obtain an Anvil HPC account, see:
#  https://allocations.access-ci.org/get-your-first-project 
###############################################################################

set -eu -o pipefail

### 1. Choose install prefix & build workspace ################################
export PREFIX=$HOME/install/meep-stack         
mkdir -p "$PREFIX"

### 2. Load Anvil module environment ##########################################
module purge
module load modtree/cpu  
# core toolchain
module load gcc/11.2.0 openmpi/4.1.6
# libraries that Meep needs (all are in the CPU stack)
module load fftw gsl hdf5 libpng gmp pcre zlib libffi
module load swig cmake    # dev tools
module load python/3.9.5
module load openblas/0.3.17 netlib-lapack/3.8.0 

### 3. Helpful build flags #####################################################
export MPICC=$(which mpicc)
export MPICXX=$(which mpicxx)
export LDFLAGS="-L$PREFIX/lib -L$PREFIX/lib64 -Wl,-rpath,$PREFIX/lib:$PREFIX/lib64"
export CPPFLAGS="-I$PREFIX/include"
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH"
export LD_LIBRARY_PATH="$PREFIX/lib:$PREFIX/lib64:$LD_LIBRARY_PATH"
export PATH="$PREFIX/bin:$PATH"

# A helper for parallel builds
NPROC=$(nproc)

### 4. Build the NanoComp stack ################################################
cd "$PREFIX"

# ---- Libtool (libltdl) ------------------------------------------------------
wget https://ftp.gnu.org/gnu/libtool/libtool-2.4.7.tar.gz
tar xf libtool-2.4.7.tar.gz && cd libtool-2.4.7

./configure --prefix="$PREFIX"           \
            --enable-ltdl-install        # installs libltdl.* + ltdl.h
make -j $NPROC && make install
cd ..

export CPPFLAGS="-I$PREFIX/include $CPPFLAGS"
export LDFLAGS="-L$PREFIX/lib $LDFLAGS"
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH"
export LD_LIBRARY_PATH="$PREFIX/lib:$LD_LIBRARY_PATH"

# ---- libunistring -----------------------------------------------------------
wget https://ftp.gnu.org/gnu/libunistring/libunistring-1.2.tar.gz
tar xf libunistring-1.2.tar.gz
cd libunistring-1.2

./configure --prefix="$PREFIX"   
make -j "$NPROC"
make install
cd ..

# ---- Boehm GC (bdw-gc) ------------------------------------------------------
wget https://github.com/ivmai/bdwgc/releases/download/v8.2.4/gc-8.2.4.tar.gz
tar xf gc-8.2.4.tar.gz
cd gc-8.2.4
./configure --prefix="$PREFIX" \
            --enable-cplusplus
make -j "$NPROC"
make install
cd ..

# ---- Guile -----------------------------------------------------------------
wget https://ftp.gnu.org/gnu/guile/guile-2.2.7.tar.gz
tar xvf guile-2.2.7.tar.gz
cd guile-2.2.7/
./configure --prefix="$PREFIX"
make -j $NPROC && make install
cd .. 

# ---- Harminv ---------------------------------------------------------------
git clone https://github.com/NanoComp/harminv.git
cd harminv
sh autogen.sh --prefix="$PREFIX" --enable-shared
make -j $NPROC && make install
cd ..

# ---- libctl -----------------------------------------------------------------
git clone https://github.com/NanoComp/libctl.git
cd libctl
sh autogen.sh --prefix="$PREFIX" --enable-shared
make -j $NPROC && make install
cd ..

# ---- h5utils ----------------------------------------------------------------
git clone https://github.com/NanoComp/h5utils.git
cd h5utils
sh autogen.sh --prefix="$PREFIX" CC="$MPICC" LDFLAGS="$LDFLAGS" CPPFLAGS="$CPPFLAGS"
make -j $NPROC && make install
cd ..

# ---- MPB --------------------------------------------------------------------
export LIBCTL_DIR="$PREFIX/share/libctl"   # convenience
git clone https://github.com/NanoComp/mpb.git
cd mpb
sh autogen.sh \
     --prefix="$PREFIX" \
     --enable-shared \
     --with-libctl="$LIBCTL_DIR" \
     CC="$MPICC" \
     CPPFLAGS="$CPPFLAGS" \
     LDFLAGS="$LDFLAGS" \
     --with-hermitian-eps
make -j $NPROC && make install
cd ..

# ---- libGDSII ---------------------------------------------------------------
git clone https://github.com/HomerReid/libGDSII.git
cd libGDSII
sh autogen.sh --prefix="$PREFIX"
make -j $NPROC && make install
cd ..

### 5. Python environment (user) ##############################################
pip install --upgrade pip wheel numpy==2.0.2 scipy matplotlib==3.3.4 mpi4py h5py

### 6. Build our modified MEEP ################################################
git clone https://github.com/TaoELi/meep.git
cd meep
sh autogen.sh --prefix="$PREFIX" \
       --enable-shared --with-mpi --disable-openmp \
       PYTHON="python" \
       MPICC="$MPICC" MPICXX="$MPICXX" \
       CPPFLAGS="$CPPFLAGS" LDFLAGS="$LDFLAGS"
make -j $NPROC && make install
cd ..

### 7. Add to ~/.bashrc ######################################################
grep -q "## MEEP-stack" ~/.bashrc 2>/dev/null || cat >> ~/.bashrc <<EOF
## MEEP-stack (Anvil) ##########################################################
module purge
module load modtree/cpu gcc/11.2.0 openmpi/4.1.6
module load fftw gsl hdf5 libpng gmp pcre zlib libffi
module load swig cmake    # dev tools
module load python/3.9.5                          # any ≥3.9 works
module load openblas/0.3.17 netlib-lapack/3.8.0

export PREFIX=$PREFIX
export PATH=\$PREFIX/bin:\$PATH
export LD_LIBRARY_PATH=\$PREFIX/lib:\$PREFIX/lib64:\$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=\$PREFIX/lib/pkgconfig:\$PKG_CONFIG_PATH
# make sure python can see where the meep is installed
export PYTHONPATH=\$PREFIX/lib/python3.9/site-packages:\$PYTHONPATH
EOF

echo "-------------------------------------------------------------------"
echo "MEEP stack installed in \$PREFIX"
echo "Open a new shell or 'source ~/.bashrc' before using."
echo "Run with:  mpirun -n <cores> python your_script.py"


