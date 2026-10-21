# Note
To run ddf-pipeline standalone (outside container):

1. export DDF_LOCAL_DEV=1
2. python3 -m venv myenv
3. source myenv/bin/activate
4. pip install -r requirements.txt
5. export DDF_PIPELINE_CATALOGS=~/CATALOGS
7. go into your dataset folder
8. python3 scripts/pipeline.py examples/mpi.cfg