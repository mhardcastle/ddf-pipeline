cd $DDF_DIR/DDFacet
git pull
python setup.py build
cd $DDF_DIR/killMS
git pull
cd Predict
make
cd ../Array/Dot
make
cd ../../Gridder
make
cd 
cd $DDF_DIR/DynSpecMS
git pull
