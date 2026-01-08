import mpi_manager
from auxcodes import run,run_serial
from mpi4py import MPI
import sys

def test_bcast():
    i=10
    print("i=",i)
    j= MPI.COMM_WORLD.bcast(i, root=0)
    print("(rank %i) j=%i"%(MPI.COMM_WORLD.rank,j))


def test_init():
    SetMS=mpi_manager.MSSet("mslist.txt")
    FullSetMS=mpi_manager.MSSet("big-mslist.txt")
    MPI_Manager=mpi_manager.mpi_manager({},SetMS, FullSetMS)
    runcommand="""python -c "from mpi4py import MPI; print(MPI.Get_processor_name())" """
    run(runcommand,
        mpiManager=MPI_Manager,
        mpi_disabled_in_serial_call=False)

def test_run_bcast():
    SetMS=mpi_manager.MSSet("mslist.txt")
    FullSetMS=mpi_manager.MSSet("big-mslist.txt")
    MPI_Manager=mpi_manager.mpi_manager({},SetMS, FullSetMS)
    runcommand="""python %s bcast """%(__file__)
    run(runcommand,
        mpiManager=MPI_Manager,
        mpi_disabled_in_serial_call=False)
    
if __name__=="__main__":
    #test_init()
    print(__file__)
    if len(sys.argv)==1:
        # test_init()
        # test_bcast()
        test_run_bcast()
    elif sys.argv[1]=="bcast":
        print("   process test bcast")
        test_bcast()
