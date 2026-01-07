import mpi_manager
from auxcodes import run,run_serial

SetMS=mpi_manager.MSSet("mslist.txt")
FullSetMS=mpi_manager.MSSet("big-mslist.txt")
MPI_Manager=mpi_manager.mpi_manager({},SetMS, FullSetMS)
runcommand="""/home/tasse/VE_MPI/venv/bin/python -c "from mpi4py import MPI; print(MPI.Get_processor_name())" """
run(runcommand,
    mpiManager=MPI_Manager,
    mpi_disabled_in_serial_call=False)

