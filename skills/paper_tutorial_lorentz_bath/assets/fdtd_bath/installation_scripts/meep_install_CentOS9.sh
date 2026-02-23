#!/bin/bash
###############################################################################
#  meep-install_CentOS9.sh  –  Build MEEP + prerequisites in $HOME with sudo
#  Author: Tao E. Li <taoeli@udel.edu>
###############################################################################

# works for a clean CentOS9 system; see also https://meep.readthedocs.io/en/latest/Build_From_Source/
RPATH_FLAGS="-Wl,-rpath,/usr/local/lib:/usr/local/lib/openmpi"
MY_LDFLAGS="-L/usr/local/lib -L/usr/local/lib/openmpi ${RPATH_FLAGS}"
MY_CPPFLAGS="-I/usr/local/include -I/usr/local/include/openmpi"

sudo dnf install epel-release epel-next-release

sudo dnf -y install   \
    bison             \
    byacc             \
    cscope            \
    ctags             \
    cvs               \
    diffstat          \
    plasma-oxygen            \
    flex              \
    gcc               \
    gcc-c++           \
    gcc-gfortran      \
    gettext           \
    git               \
    indent            \
    intltool          \
    libtool           \
    patch             \
    patchutils        \
    rcs               \
    redhat-rpm-config \
    rpm-build         \
    subversion        \
    systemtap         \
    wget


# enable PRM fusion free
sudo dnf install \
  https://download1.rpmfusion.org/free/el/rpmfusion-free-release-$(rpm -E %rhel).noarch.rpm \
  https://download1.rpmfusion.org/nonfree/el/rpmfusion-nonfree-release-$(rpm -E %rhel).noarch.rpm \
  -y

sudo dnf -y install    \
    openblas-devel     \
    fftw3-devel        \
    libpng-devel       \
    gsl-devel          \
    gmp-devel          \
    pcre-devel         \
    libtool-ltdl-devel \
    libunistring-devel \
    libffi-devel       \
    gc-devel           \
    zlib-devel         \
    openssl-devel      \
    sqlite-devel       \
    bzip2-devel        \
    ffmpeg             \
    ffmpeg-devel




mkdir -p ~/install

# install swig
cd ~/install
wget https://github.com/swig/swig/archive/rel-3.0.12.tar.gz
tar xvf rel-3.0.12.tar.gz
cd swig-rel-3.0.12
./autogen.sh
./configure
make -j
sudo make -j install

# install guile
cd ~/install
wget https://ftp.gnu.org/gnu/guile/guile-2.2.7.tar.gz
tar xvf guile-2.2.7.tar.gz
cd guile-2.2.7/
./configure
make -j
sudo make -j install

# install python3.9
sudo dnf -y install python3.9 python3.9-pip

# install openmpi
cd ~/install
wget https://download.open-mpi.org/release/open-mpi/v2.1/openmpi-2.1.1.tar.gz
tar xvf openmpi-2.1.1.tar.gz
cd openmpi-2.1.1/
./configure
make -j all
sudo make -j install

# install hdf5
sudo dnf -y install hdf5 hdf5-devel hdf5-openmpi hdf5-openmpi-devel

# install harminv
cd ~/install
git clone https://github.com/NanoComp/harminv.git
cd harminv/
sh autogen.sh --enable-shared
make -j
sudo make -j install

# install libctl
cd ~/install
git clone https://github.com/NanoComp/libctl.git
cd libctl/
sh autogen.sh  --enable-shared
make -j
sudo make -j install

# install h5utils
cd ~/install
git clone https://github.com/NanoComp/h5utils.git
cd h5utils/
sh autogen.sh CC=/usr/local/bin/mpicc LDFLAGS="${MY_LDFLAGS}" CPPFLAGS="${MY_CPPFLAGS}"
make -j
sudo make -j install

# install MPB
cd ~/install
git clone https://github.com/NanoComp/mpb.git
cd mpb/
sh autogen.sh --enable-shared CC=/usr/local/bin/mpicc LDFLAGS="${MY_LDFLAGS}" CPPFLAGS="${MY_CPPFLAGS}" --with-hermitian-eps
make -j
sudo make -j install

# install libGDSII
cd ~/install
git clone https://github.com/HomerReid/libGDSII.git
cd libGDSII/
sh autogen.sh
sudo make -j install

# install mpi4py and h5py
python3 -m pip install --no-cache-dir mpi4py
python3 -m pip install --no-cache-dir h5py

# install modified meep code
cd ~/install
git clone git@github.com:TaoELi/meep.git
cd meep/
# please disable openmp
sh autogen.sh --enable-shared --with-mpi PYTHON=python3 MPICC=/usr/local/bin/mpicc MPICXX=/usr/local/bin/mpic++ LDFLAGS="${MY_LDFLAGS}" CPPFLAGS="${MY_CPPFLAGS}"
make -j
sudo make install

# finally write PYTHONPATH to bashrc
echo 'export PYTHONPATH="$HOME/install/meep/python/meep:$PYTHONPATH"' >> ~/.bashrc
