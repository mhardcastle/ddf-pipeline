export DDF_PIPELINE_CATALOGS=/home/$USER/CATALOGS
export DDF_LOCAL_DEV=1
export  OPENBLAS_NUM_THREADS=1
export OPENBLAS_MAX_THREADS=1

export PIPELINE_EXE="$(command -v pipeline.py)"
export CLEAN_EXE="$(command -v cleanup.sh)"
export OMPI_MCA_btl_tcp_if_include="$(python -c "import subprocess,re;print(','.join(i for i in re.findall(r'^([A-Za-z0-9]+):', subprocess.check_output(['ifconfig'], text=True), re.M) if i!='lo' and not i.startswith(('docker','br-'))))")"
echo $OMPI_MCA_btl_tcp_if_include

HOSTS="$(python -c "print(''.join(f'{h}:1,' for h in sorted({l.split(':',1)[0] for l in open('big-mslist.txt')})))")"
HOSTS=${HOSTS%,}   # remove trailing comma

NODES=$(echo "$HOSTS" | tr ',' '\n' | wc -l)

if [ "$1" = "Clean" ]; then
export EXEC=$CLEAN_EXE
else
export EXEC="python -m mpi4py.futures $PIPELINE_EXE $1"
fi


EXP="mpirun -np $NODES -x OPENBLAS_NUM_THREADS -x OPENBLAS_MAX_THREADS -x DDF_PIPELINE_CATALOGS -x DDF_LOCAL_DEV -x PATH -x VIRTUAL_ENV -x LD_LIBRARY_PATH -x PYTHONPATH -wdir $PWD -host $HOSTS $EXEC"
echo "====================="
echo running: $EXP
echo "====================="

$EXP
