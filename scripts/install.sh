#!/bin/bash

export WD=`pwd`
git clone https://github.com/cyriltasse/SkyModel.git
git clone https://github.com/saopicc/killMS.git
cd killMS
git checkout lofar-stable
git tag -a 'v1.0' -m 'Dummy tag'
cd Predict
make
cd ../Array/Dot
make
cd ../../Gridder
make
cd $WD
git clone https://github.com/saopicc/DDFacet.git
cd DDFacet
git checkout lofar-stable
python setup.py build
cd $WD
sed -e "s|INSTALLDIR|$WD|" ddf-pipeline/misc/DDF.sh > init.sh
ln -s killMS killMS2
