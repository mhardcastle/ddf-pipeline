import mpi_manager
from auxcodes import run,run_serial
from mpi4py import MPI
import mpi4py
import sys
import os

def test_bcast():
    i=10
    print(MPI.Get_processor_name())
    print(mpi4py.__version__,mpi4py.__file__)
    print("i=",i)
    os.system("mpirun --version")
    jsum=MPI.COMM_WORLD.allreduce(i, MPI.SUM)
    j= MPI.COMM_WORLD.bcast(i, root=0)
    print("(rank %i) j=%i sumj=%i"%(MPI.COMM_WORLD.rank,j,jsum))


def test_init(MPI_Manager):
    runcommand="""python -c "from mpi4py import MPI; print(MPI.COMM_WORLD.rank,MPI.Get_processor_name())" """
    run(runcommand,
        mpiManager=MPI_Manager,
        mpi_disabled_in_serial_call=False)

def test_run_bcast(MPI_Manager):
    runcommand="""python %s bcast """%(__file__)
    run(runcommand,
        mpiManager=MPI_Manager,
        mpi_disabled_in_serial_call=False)
    
if __name__=="__main__":
    #test_init()
    # test_bcast()
    if len(sys.argv)==1:
        SetMS=mpi_manager.MSSet("mslist.txt")
        FullSetMS=mpi_manager.MSSet("big-mslist.txt")
        MPI_Manager=mpi_manager.mpi_manager({},SetMS, FullSetMS)
        test_init(MPI_Manager)
        # test_run_bcast(MPI_Manager)
    elif sys.argv[1]=="bcast":
        print("   process test bcast")
        test_bcast()
