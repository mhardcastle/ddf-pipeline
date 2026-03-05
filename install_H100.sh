#!/bin/sh
module purge

module load arch/h100
module load anaconda-py3/2024.06
module load gcc/11.4.1
module load openmpi/4.1.6
module load cmake/3.21.3
module load gsl/2.7.1

conda create -n ddf_mpi python=3.10
conda activate ddf_mpi

PYVER=$($CONDA_PREFIX/bin/python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")

mkdir -p casacore_data
cd casacore_data
wget https://www.astron.nl/iers/WSRT_Measures.ztar
tar -zxvf WSRT_Measures.ztar
export CASA_DATA=$PWD
cd ..

pip install --upgrade pip
pip install numpy==1.22.4
export CFLAGS="-I$(python - <<'EOF'
import numpy
print(numpy.get_include())
EOF
)"
pip install sharedarray==3.2.1
pip install scikit-build-core
conda install -c conda-forge casacore==3.4.0 python-casacore pybind11 setuptools==65.6.3

mkdir -p $CONDA_PREFIX/lib/casa/data/ephemerides
mkdir -p $CONDA_PREFIX/lib/casa/data/geodetic

cp -r $CASA_DATA/ephemerides/* $CONDA_PREFIX/lib/casa/data/ephemerides/
cp -r $CASA_DATA/geodetic/* $CONDA_PREFIX/lib/casa/data/geodetic/

mkdir sources && cd sources
git clone -b NoBoost  https://github.com/cyriltasse/LOFARBeam.git
    cd LOFARBeam
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
      -DPYTHON_EXECUTABLE=$CONDA_PREFIX/bin/python \
      -DPYTHON_INSTALL_DIR=$CONDA_PREFIX/lib/python${PYVER}/site-packages \
      ..
    make -j 10
    make install
    cd ../..

git clone https://github.com/aroffringa/dysco.git
    cd dysco
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
      -DPYTHON_EXECUTABLE=$CONDA_PREFIX/bin/python \
      -DPYTHON_INSTALL_DIR=$CONDA_PREFIX/lib/python${PYVER}/site-packages \
      ..
    make -j 10
    make install
    cd ../..

export DDFACET_BRANCH=MassiveMerge_PR_MergeSSD3_NancepMPI_APP
export KMS_BRANCH=APP_Predict_Compress_PolSmooth_HybridSM_OpFit_MultiField_MPI_MultiChain
export DDFPIPE_BRANCH=Hackaton_mpipool_test_NancepMPI_Herts

git clone -b $DDFACET_BRANCH https://github.com/cyriltasse/DDFacet.git
pip install --no-build-isolation -e ./DDFacet[mpi-support]

git clone -b $KMS_BRANCH https://github.com/cyriltasse/killMS.git
pip install -e ./killMS

git clone -b $DDFPIPE_BRANCH https://github.com/besnardb/ddf-pipeline.git

pip uninstall mpi4py -y
pip install --no-cache-dir --no-binary=mpi4py mpi4py==4.1.1

pip install future
pip install sshtunnel
pip install pymysql
pip install pyregion
pip install ipython