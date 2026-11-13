cd
python3.12 -m venv VE_P3.12_Merge
cd VE_P3.12_Merge

source bin/activate
#cp $HOME/VE_P3.12_Merge/init.sh .

cd $VIRTUAL_ENV
git clone -b HackathonRennes_June26 git@github.com:cyriltasse/DDFacet.git
git clone -b HackathonRennes_June26 git@github.com:cyriltasse/killMS.git
git clone -b HackathonRennes_June26 git@github.com:mhardcastle/ddf-pipeline.git
pip install -e ./DDFacet[mpi-support]
pip install -e ./killMS
pip install -e ./ddf-pipeline

# install LOFAR beam
wget https://github.com/bennahugo/LOFARBeam/archive/refs/tags/DDF_KMS_22.04.tar.gz
tar zxvf DDF_KMS_22.04.tar.gz
cd LOFARBeam-DDF_KMS_22.04
mkdir -p build/gnucxx11_opt
cd build/gnucxx11_opt
cmake -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV       -DPYTHON_PACKAGES_DIR=$VIRTUAL_ENV/lib/python3.10/site-packages       ../..
make -j4

# install oskar
pip install 'git+https://github.com/OxfordSKA/OSKAR.git@master#egg=oskarpy&subdirectory=python'


# install oskar
cd $VIRTUAL_ENV
git clone https://github.com/OxfordSKA/OSKAR.git
cd OSKAR/oskar/
mkdir build
cd build/
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV/OSKAR/install -DFIND_CUDA=OFF -DFIND_OPENCL=OFF
make -j$(nproc)
make install
emacs ~/VE_P3.12_Merge/init.sh &
cd ../python
pip install .


cd $VIRTUAL_ENV
