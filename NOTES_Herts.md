
module load openmpi-5.0.5

mkdir -p casacore_data
cd casacore_data
wget https://www.astron.nl/iers/WSRT_Measures.ztar
tar -zxvf WSRT_Measures.ztar

cd ..

python3.9 -m venv venv
source venv/bin/activate
export VE_FOLDER=$PWD/venv

mkdir -p pip-cache pip-tmp
export PIP_CACHE_DIR="$PWD/pip-cache"
export TMPDIR="$PWD/pip-tmp"
pip install --upgrade pip

pip install numpy==1.22.4
export CFLAGS="-I$(python - <<'EOF'
import numpy
print(numpy.get_include())
EOF
)"
pip install sharedarray==3.2.1

mkdir sources
cd $VE_FOLDER/../sources
git clone https://github.com/sailfishos-mirror/readline.git
cd readline/
./configure --prefix=$VE_FOLDER
make -j 20
make install


cd $VE_FOLDER/../sources
wget https://archives.boost.io/release/1.74.0/source/boost_1_74_0.tar.gz
tar -zxvf boost_1_74_0.tar.gz 
cd boost_1_74_0/
./bootstrap.sh --prefix=$VE_FOLDER --with-python=/usr/bin/python3.9  variant=release
./b2 --with-python variant=release threading=multi link=shared runtime-link=shared install

cd $VE_FOLDER/../sources
git clone https://github.com/healpy/cfitsio.git
cd cfitsio/
./configure --prefix=$VE_FOLDER
make -j 20
make install

cd $VE_FOLDER/../sources
git clone https://github.com/casacore/casacore.git
cd casacore/
git checkout v3.4.0
mkdir cbuild
cd cbuild
cmake -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV -DBUILD_PYTHON=OFF -DBUILD_PYTHON3=ON -DREADLINE_LIBRARY=$VE_FOLDER/lib/libreadline.so.8.3 -DNCURSES_LIBRARY=/usr/lib64/libncurses.so.6 -DBUILD_TESTING=OFF -DDATA_DIR=/home/$USER/casacore_data -DBOOST_ROOT=$VIRTUAL_ENV ..
make -j 20
make install

cd $VE_FOLDER/../sources
git clone -b v3.5.2 https://github.com/casacore/python-casacore.git
CASACORE_DATA=/home/$USER/casacore_data CMAKE_ARGS="-DCASACORE_ROOT_DIR=$VIRTUAL_ENV -DCFITSIO_INCLUDE_DIR=$VIRTUAL_ENV/include -DCFITSIO_LIBRARY=$VIRTUAL_ENV/lib" pip install -e ./python-casacore


cd $VE_FOLDER/../sources
git clone -b NoBoost https://github.com/cyriltasse/LOFARBeam.git
mkdir -p LOFARBeam/build
cd LOFARBeam/build
cmake -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV ..
make -j 20
make install

cd $VE_FOLDER/../sources
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.3/src/hdf5-1.14.3.tar.gz
tar -xzf hdf5-1.14.3.tar.gz
cd hdf5-1.14.3
./configure --prefix=$VIRTUAL_ENV --enable-cxx
make -j 20
make install

cd $VE_FOLDER/../sources
git clone -b v0.5.8 https://git.astron.nl/RD/EveryBeam.git
cd EveryBeam/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV -DBUILD_WITH_PYTHON=On -DCASACORE_ROOT_DIR=$VIRTUAL_ENV -DCFITSIO_INCLUDE_DIR=$VIRTUAL_ENV/include -DCFITSIO_LIBRARY=$VIRTUAL_ENV/lib/libcfitsio.so ..
make -j 20

cd $VE_FOLDER/../sources
git clone https://github.com/aroffringa/dysco.git
cd dysco
git checkout 3fd7a5fd17f3d09db89ad7827c9bdc4febf66eff
mkdir build
cd build
cmake -DCASACORE_ROOT_DIR=$VIRTUAL_ENV  -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV ../
make -j $J
make install

cd $VE_FOLDER/../sources
export DDFACET_BRANCH=MassiveMerge_PR_MergeSSD3_NancepMPI_APP
export KMS_BRANCH=APP_Predict_Compress_PolSmooth_HybridSM_OpFit_MultiField_MPI_MultiChain
export DDFPIPE_BRANCH=Hackaton_mpipool_test_NancepMPI_Herts

git clone -b $DDFACET_BRANCH https://github.com/cyriltasse/DDFacet.git
pip install -e ./DDFacet[mpi-support]

git clone -b $KMS_BRANCH https://github.com/cyriltasse/killMS.git
pip install -e ./killMS

git clone -b $DDFPIPE_BRANCH https://github.com/mhardcastle/ddf-pipeline.git
pip install future sshtunnel mysqldb pyregion

git clone https://github.com/mhardcastle/lotss-query.git

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/beegfs/car/mjh/cyril-version/venv/lib
export PYTHONPATH=`pwd`/sources/ddf-pipeline/scripts:`pwd`/sources/ddf-pipeline/utils:`pwd`/sources/lotss-query


mpirun -np 2 -x DDF_PIPELINE_CATALOGS -x DDF_LOCAL_DEV -x PATH -x VIRTUAL_ENV -x LD_LIBRARY_PATH -x PYTHONPATH -host smp5:1,smp4:1 DDF.py -h

