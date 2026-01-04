export DDF_PIPELINE_CATALOGS=/home/$USER/CATALOGS
export DDF_LOCAL_DEV=1
export PIPELINE_EXE="$(command -v pipeline.py)"
export CLEAN_EXE="$(command -v cleanup.sh)"
export OMPI_MCA_btl_tcp_if_include="$(ifconfig | awk '/^[a-zA-Z0-9]/ {print $1}' | sed 's/://' | grep -Ev '^(lo|docker|br-)' | paste -sd, -)"
HOSTS=$(cut -d: -f1 big-mslist.txt | sort -u | awk '{printf "%s:1,", $0}')
HOSTS=${HOSTS%,}   # remove trailing comma

NODES=$(echo "$HOSTS" | tr ',' '\n' | wc -l)

if [ "$1" = "Clean" ]; then
export EXEC=$CLEAN_EXE
else
export EXEC="python -m mpi4py.futures $PIPELINE_EXE $1"
fi

EXP="mpirun -np $NODES -x DDF_PIPELINE_CATALOGS -x DDF_LOCAL_DEV -x PATH -x VIRTUAL_ENV -x LD_LIBRARY_PATH -x PYTHONPATH -wdir $PWD -host $HOSTS $EXEC"
echo $EXP

$EXP
