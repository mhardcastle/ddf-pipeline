import mpi_manager
from auxcodes import run,run_serial
from mpi4py import MPI
import mpi4py
import sys
import os
import time

def test_bcast():
    i=10
    # os.system("mpirun --version")
    print(mpi4py.__version__,mpi4py.__file__)
    print("%s (rank %i, size=%i) i=%i"%(MPI.Get_processor_name(),MPI.COMM_WORLD.rank,MPI.COMM_WORLD.size,i))
    print("Barrier")
    # MPI.COMM_WORLD.Barrier()
    print("Barrier ok")
    time.sleep(10)
    print("Barrier ok2")
    jsum=MPI.COMM_WORLD.allreduce(i, MPI.SUM)
    j= MPI.COMM_WORLD.bcast(i, root=0)
    print("%s (rank %i) j=%i sumj=%i"%(MPI.Get_processor_name(),MPI.COMM_WORLD.rank,j,jsum))


def test_init(MPI_Manager):
    runcommand="""python -c "from mpi4py import MPI; print(MPI.Get_processor_name())" """
    run(runcommand,
        mpiManager=MPI_Manager,
        mpi_disabled_in_serial_call=False)

def test_run_bcast(MPI_Manager):
    runcommand="""python %s bcast """%(__file__)
    run(runcommand,
        mpiManager=MPI_Manager,
        mpi_disabled_in_serial_call=False)

def print_name():
    print(MPI.Get_processor_name())
    
def test_callparallel_bcast():
    ListJobs=[
        ["node081",print_name,(), {}],
        ["node082",print_name,(), {}],
    ]
    mpi_manager.callParallel(ListJobs)
    
    # ListJobs=[["node081",test_bcast,(), {}],
    #           ["node082",test_bcast,(), {}],
    #           ]
    # mpi_manager.callParallel(ListJobs)
    
    ListJobs=[["node081",test_run_bcast,(), {}],
              ["node082",test_run_bcast,(), {}],
              ]
    mpi_manager.callParallel(ListJobs)


    
if __name__=="__main__":
    #test_init()
    # test_bcast()
    if len(sys.argv)==1:
        # test_callparallel_bcast()
        SetMS=mpi_manager.MSSet("mslist.txt")
        FullSetMS=mpi_manager.MSSet("big-mslist.txt")
        MPI_Manager=mpi_manager.mpi_manager({},SetMS, FullSetMS)
        test_init(MPI_Manager)
        test_run_bcast(MPI_Manager)
    elif sys.argv[1]=="bcast":
        print("   process test bcast")
        test_bcast()
