#!/bin/bash
# Install additional dependencies

# Install the catalogs
wget https://sdrive.cnrs.fr/public.php/dav/files/awWCRsaoeax8t4K/CATALOGS.zip
unzip CATALOGS.zip

# Install the simulated dataset
wget -O Sim.tgz https://sdrive.cnrs.fr/public.php/dav/files/PNmRiW3zc7XmEwN
tar -xvf Sim.tgz

# Install the ddf pipeline configuration file working with the simulated dataset
wget -O Simul-Hackaton.cfg https://sdrive.cnrs.fr/public.php/dav/files/qgptXzGa6BCE56J
mv Simul-Hackaton.cfg examples/
