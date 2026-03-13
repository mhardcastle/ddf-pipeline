#!/bin/sh
module purge

module load arch/h100
module load anaconda-py3/2024.06
module load gcc/11.4.1
module load openmpi/4.1.6
module load cmake/3.21.3
module load gsl/2.7.1

# Clear cache before installation
conda clean --all
pip cache purge

# Changing the Default Package Directory, avoid Disk quota exceeded
#rm -r /lustre/fswork/projects/rech/doz/udd71uc/packages
#mkdir -p /lustre/fswork/projects/rech/doz/udd71uc/packages
export DATA=/lustre/fswork/projects/rech/doz/udd71uc/packages
echo $DATA
conda config --add envs_dirs $DATA/.conda/envs
conda config --add pkgs_dirs $DATA/.conda/pkgs

conda create -n ddf_env-0903 python=3.9
conda activate ddf_env-0903
echo $CONDA_PREFIX

rm -r VE_MPI
mkdir VE_MPI && cd VE_MPI
git clone https://github.com/cyriltasse/MPIpeline.git

rm -r casacore_data
mkdir -p casacore_data

cd casacore_data
wget https://www.astron.nl/iers/WSRT_Measures.ztar
tar -zxvf WSRT_Measures.ztar
export CASA_DATA=$PWD
cd ..

mkdir -p $SCRATCH/tmp/pip-cache $SCRATCH/tmp/pip-tmp
export PIP_CACHE_DIR=$SCRATCH/tmp/pip-cache
export TMPDIR=$SCRATCH/tmp/pip-tmp
pip install --upgrade pip

pip install numpy==1.22.4
export CFLAGS="-I$(python - <<'EOF'
import numpy
print(numpy.get_include())
EOF
)"
pip install sharedarray==3.2.1

conda install -c conda-forge casacore==3.4.0 python-casacore setuptools==65.6.3 pybind11
pip install scikit-build-core

mkdir -p $CONDA_PREFIX/lib/casa/data/ephemerides
mkdir -p $CONDA_PREFIX/lib/casa/data/geodetic

cp -r $CASA_DATA/ephemerides/* $CONDA_PREFIX/lib/casa/data/ephemerides/
cp -r $CASA_DATA/geodetic/* $CONDA_PREFIX/lib/casa/data/geodetic/

rm -r sources
mkdir sources && cd sources

git clone --branch NoBoost https://github.com/cyriltasse/LOFARBeam.git
mkdir -p LOFARBeam/build
cd LOFARBeam/build
cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
      -DPYTHON_EXECUTABLE=$CONDA_PREFIX/bin/python \
      -DPYTHON_INSTALL_DIR=$CONDA_PREFIX/lib/python3.9/site-packages \
      ..
make -j 10 
make install
cd ../..

export DDFACET_BRANCH=MassiveMerge_PR_MergeSSD3_NancepMPI_APP
export KMS_BRANCH=APP_Predict_Compress_PolSmooth_HybridSM_OpFit_MultiField_MPI_MultiChain
export DDFPIPE_BRANCH=Hackaton_mpipool_test_NancepMPI_Herts
export DDFACET_COMMIT=403050b69c51e77ac9b7db69e6a69201db785a27
export KMS_COMMIT=d9de06d57bd5f083f41ccf3dae85588bbc7b33df
export DDFPIPE_COMMIT=109e389d0fb523d282c321ceb0cfe884bfe98004

git clone --branch $DDFACET_BRANCH https://github.com/cyriltasse/DDFacet.git
pip install --no-build-isolation -e ./DDFacet[mpi-support]

git clone --branch $KMS_BRANCH https://github.com/cyriltasse/killMS.git
pip install -e ./killMS

git clone --branch $DDFPIPE_BRANCH https://github.com/besnardb/ddf-pipeline.git

cd ddf-pipeline
git checkout $DDFPIPE_COMMIT
cd ..

cd DDFacet
git checkout $DDFACET_COMMIT
cd ..
cd killMS
git checkout $KMS_COMMIT
cd ..

cd ddf-pipeline/examples/

## Configuration file
wget -O lotss_mpi.cfg https://sdrive.cnrs.fr/public.php/dav/files/omnAgeZL5cPHAAn
wget -O Simul-Hackaton.cfg https://sdrive.cnrs.fr/public.php/dav/files/qgptXzGa6BCE56J
cd ../..

pip uninstall mpi4py -y
python -m pip install --no-cache-dir --no-binary=mpi4py mpi4py==4.1.1

pip install future
pip install sshtunnel
pip install pymysql
pip install pyregion
pip install ipython

## Dysco
cd /lustre/fsn1/projects/rech/doz/udd71uc/Dev/ddf-mpi/VE_MPI
git clone https://github.com/aroffringa/dysco.git
cd dysco
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
  -DPYTHON_EXECUTABLE=$CONDA_PREFIX/bin/python \
  -DPYTHON_INSTALL_DIR=$CONDA_PREFIX/lib/python3.9/site-packages \
  ..
make -j 10
make install

## Catalogs
cd /lustre/fsn1/projects/rech/doz/udd71uc/Dev/ddf-mpi/VE_MPI
wget -O CATALOGS.zip https://sdrive.cnrs.fr/public.php/dav/files/qDa8owWokGEaGBK
unzip CATALOGS.zip
rm -r CATALOGS.zip

