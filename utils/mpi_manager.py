MPIsize = 0
RANK=0
localrank = 0 # local process ID on the node
localsize = 1 # number of mpi process on the node
hostname = "localhost"
USE_MPI=False

try:
    from mpi4py import MPI
    MPIsize = MPI.COMM_WORLD.size
    RANK=MPI.COMM_WORLD.rank
    shared_comm = MPI.COMM_WORLD.Split_type(MPI.COMM_TYPE_SHARED)
    localrank = shared_comm.rank
    localsize = shared_comm.size
    hostname = MPI.Get_processor_name()

    if MPIsize>1:
        USE_MPI=True
    else:
        print(" mpi4py properly initialised, but size=1, not using MPI mode.")

except ModuleNotFoundError:
    print(" Could not initialise mpi4py ")
    pass
except Error as e:
    raise RuntimeError(e)


LIST_SITES_BEING_USED=None

if USE_MPI:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    host = f"{hostname}@{localrank}"
    print(f"{host}")
    hosts_rank=comm.allgather(host)
    print(f"{hosts_rank}")
    LIST_SITES_BEING_USED=hosts_rank
    print(f"LIST_SITES_BEING_USED={LIST_SITES_BEING_USED}")


import itertools
import os
from DDFacet.Other import ModColor
import pprint
from pathlib import Path

WIDTH_PROMPT=90


def CheckDistributedFS():
    cwd = os.getcwd()
    cwd = os.path.realpath(cwd)
    print(cwd)

    best_match = ""
    fs_type = None

    with open("/proc/self/mounts") as f:
        for line in f:
            parts = line.split()
            # print("=================")
            # print(parts)
            if len(parts) < 3:
                continue
            mount_point = parts[1]
            fstype = parts[2]
            # print("  ",mount_point,fstype)
            if cwd.startswith(mount_point) and len(mount_point) > len(best_match):
                best_match = mount_point
                fs_type = fstype
                # print("   --- ",best_match,fs_type)
    # print("!!!!!",fs_type)

    if fs_type is not None:
        return fs_type.startswith("nfs") or fs_type== "beegfs" or fs_type == "lustre" or fs_type == "gpfs" or fs_type == "ceph"

class MSSet():
    def __init__(self,mslist):
        self.file_nodes_mslist=mslist
        # for MPI use
        mslist=[s.strip() for s in open(mslist).readlines()]

        self.DicoMSName2Node = {}
        nodes2ms = {}
        for iMS,sMS in enumerate(mslist):
            if ":" in sMS:
                node,msname=sMS.split(":")
                # if no core specified, assume @0
                if len(node.split("@")) == 1:
                    node = f"{node}@0"
            else:
                node,msname=None,sMS
            l=nodes2ms.get(node,[])
            l.append(msname)
            nodes2ms[node]=l
            self.DicoMSName2Node[msname]=node

        if (None in nodes2ms.keys()) and len(nodes2ms) == 1 and MPIsize>1:
            # get all node names because None have been specified

            self.ListNodesBeingUsed = LIST_SITES_BEING_USED

            mslist=list(zip(itertools.cycle(LIST_SITES_BEING_USED), nodes2ms[None]))
            del nodes2ms[None]

            for node,ms in mslist:
                l=nodes2ms.get(node,[])
                l.append(ms)
                nodes2ms[node] = l
                self.DicoMSName2Node[ms] = node

        self.DicoNodes2ListMS=nodes2ms
        self.ListNodesBeingUsed = list(nodes2ms.keys())

        print(ModColor.Str(" Data set distribution ".center(WIDTH_PROMPT,"="),col="blue"))
        print(ModColor.Str((" %s "%(str(self.ListNodesBeingUsed))).center(WIDTH_PROMPT,"="),col="blue"))
        pprint.pp(self.DicoNodes2ListMS)


def testFunc(*args,**kwargs):
    host = hostname
    print(host,args,kwargs)

import os
def testParallel():
    ListJobs=[["nancep10.obs-nancay.fr@0",testFunc,(5,),{"e":6}],
              ["nancep10.obs-nancay.fr@0",testFunc,(9,),{"f":6}],
              ["nancep11.obs-nancay.fr@0",testFunc,(77,),{"g":88}],
              ]

    ListJobs=[["cw10055@0",os.system,("CleanSHM.py",), {}],
              ["cw10057@0",os.system,("CleanSHM.py",), {}],
              ]

    ListJobs=[["nancep10.obs-nancay.fr@0",os.system,("CleanSHM.py",), {}],
              ["nancep11.obs-nancay.fr@0",os.system,("CleanSHM.py",), {}],
              ]

    import DDFacet.CleanSHM
    ListJobs=[["cw10042@0",DDFacet.CleanSHM.driver,(), {}],
              ["cw10042@1",DDFacet.CleanSHM.driver,(), {}],
              ]

    callParallel(ListJobs)

def filterHost(jobs):
    site = f"{hostname}@{localrank}"
    #print("FilterHost : %s"%site)
    res=[]
    for RunOnHost, func, args, kwargs in jobs:
        if RunOnHost == f"{site}":
            #print("  [exec] [ME=%s][TARGET=%s]: %s(%s,%s)"%(site,RunOnHost,str(func),str(args),str(kwargs)))
            res.append(func(*args,**kwargs))
        else:
            #print("  [skip] [ME=%s][TARGET=%s]: %s(%s,%s)"%(site,RunOnHost,str(func),str(args),str(kwargs)))
            pass
    return res


def callParallel(ListJobs):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    comm.Barrier()

    if rank == 0:
        print()
        print("".center(WIDTH_PROMPT,"="))
        print(ModColor.Str(" CALL PARALLEL ".center(WIDTH_PROMPT,"="),col="blue"))
        print(ModColor.Str((f" MPI size: {size} ").center(WIDTH_PROMPT,"="),col="blue"))
        print(ModColor.Str((f" : {ListJobs} ")))

    local_results = filterHost(ListJobs)

    all_results = comm.gather(local_results, root=0)

    # Synchronization
    comm.Barrier()

    if rank == 0:
        flat = []
        for r in all_results:
            flat.extend(r)
        print(f"callParallel: ok ({flat})")
        print()
        return flat
    else:
        return None

class mpi_manager():
    def __init__(self,options_cfg,MSSet, FullMSSet):
        self.options=options_cfg
        self.MSSet=MSSet
        self.FullMSSet=FullMSSet
        # check subset consistency


        self.ListNodesBeingUsed=FullMSSet.ListNodesBeingUsed if FullMSSet else MSSet.ListNodesBeingUsed
        self.DicoNodes2WorkDir={}
        self.WorkDir=os.getcwd()
        self.MainSite = f"{hostname}@0"
        self.ddf_nproc = int(self.options.get('ddf_nproc', 1))
        self.UseMPI=False
        self.MPIsize=MPIsize

        if MPIsize>1 and (self.ddf_nproc > 1 or self.ListNodesBeingUsed):
            self.UseMPI=True

        # Scatter mslist and big-mslist.txt

        self.DoScatterGather=(CheckDistributedFS()==False)
        self.scpScatter(MSSet.file_nodes_mslist)
        if FullMSSet is not None: self.scpScatter(FullMSSet.file_nodes_mslist)
        self.createRemoteLocal_mslist()
        self.createRemoteLocal_fullmslist()
        for Node in self.DicoNode2mslist.keys():
            LMS=self.FullMSSet.DicoNodes2ListMS[Node]
            for MSName in self.MSSet.DicoNodes2ListMS[Node]:
                if MSName not in LMS: raise Exception(f"MSs lists not consistent: {MSName} not in {LMS}")


    def callParallel(self,*args,**kwargs):
        return callParallel(*args,**kwargs)

    def givePrefixDDF(self):
        if not self.UseMPI: return ""

        if mpirun_singularity:
            return self.givePrefixDDF_mpirun_singularity()
        elif mpirun:
            return "mpirun -n %i "%(self.ddf_nproc)
        elif srun_singularity:
            pass

    def createRemoteLocal_mslist(self):
        self.DicoNode2mslist={}
        print(self.MSSet.DicoNodes2ListMS)
        for Node in self.MSSet.DicoNodes2ListMS.keys():
            Listms=self.MSSet.DicoNodes2ListMS[Node]
            FName="local_%s_mslist.txt"%Node
            f=open(FName,"w")
            for msname in Listms:
                f.write("%s\n"%msname)
            f.close()
            self.DicoNode2mslist[Node]=FName
            self.scpScatter(FName,Node)

    def createRemoteLocal_fullmslist(self):
        self.DicoNode2fullmslist={}
        if self.FullMSSet:
            for Node in self.FullMSSet.DicoNodes2ListMS.keys():
                Listms=self.FullMSSet.DicoNodes2ListMS[Node]
                FName="local_%s_full_mslist.txt"%Node
                f=open(FName,"w")
                for msname in Listms:
                    f.write("%s\n"%msname)
                f.close()
                self.DicoNode2fullmslist[Node]=FName
                self.scpScatter(FName,Node)

    def scpScatter(self,FileName,NodeDest="all"):
        if not self.UseMPI: return
        if not self.DoScatterGather: return
        rank = comm.Get_rank()
        if rank != 0: return
        if local_rank != 0: return
        if NodeDest=="all":
            for site in self.ListNodesBeingUsed:
                if site==self.MainSite:
                    continue
                Node=site.split('@')[0]
                lrank=site.split('@')[1]
                if lrank != 0: return
                ss="scp -r %s %s:%s"%(FileName,Node,self.WorkDir)
                print("[Scatter] %s"%ss)
                os.system("%s > /dev/null 2>&1"%ss)
        else:
            if NodeDest==self.MainSite: return
            ss="scp -r %s %s:%s"%(FileName,NodeDest.split("@")[0],self.WorkDir)
            print("[Scatter] %s"%ss)
            os.system("%s > /dev/null 2>&1"%ss)

    def scpGatherSolutions(self,SolName,DestDir="",NodeSource="all"):
        if not self.UseMPI: return
        if not self.DoScatterGather: return
        if local_rank != 0: return

        SolsDir=self.options["SolsDir"]
        AbsSolsDir=os.path.abspath(SolsDir)

        for site in self.ListNodesBeingUsed:
            Node=site.split('@')[0]
            lrank=site.split('@')[1]
            if lrank != 0: return
            if Node==self.MainHost:
                continue

            LMS=self.FullMSSet.DicoNodes2ListMS[Node]


            for MSName in LMS:
                MSName = Path(MSName).name # if MSName is given with full path
                os.system("mkdir -p %s/%s"%(SolsDir,MSName))
                ss="scp -r %s:%s/%s/\\*.%s.\\* %s/%s/%s"%(Node,AbsSolsDir,MSName,SolName,self.WorkDir,SolsDir,MSName)
                print("[Gather] %s"%ss)
                os.system("%s > /dev/null 2>&1"%ss)
                #os.system("%s"%ss)

    def scpScatterSolutions(self,MSName,SmoothSolName,SolsAliasName):
        if not self.UseMPI: return
        if not self.DoScatterGather: return
        if local_rank != 0: return
        SolsDir=self.options["SolsDir"]
        AbsSolsDir=os.path.abspath(SolsDir)
        Node=self.FullMSSet.DicoMSName2Node[MSName]
        Node =  Node.split("@")[0]

        MSName = Path(MSName).name # if MSName is given with full path
        os.system("mkdir -p %s/%s"%(SolsDir,MSName))
        if Node!=self.MainHost:
            ss="scp -r %s %s:%s"%(SmoothSolName,Node,self.WorkDir)
            print("[Scatter Sols] %s"%ss)
            #os.system("%s &>/dev/null"%ss)
            os.system("%s > /dev/null 2>&1"%ss)

            ss="ssh %s rm -f %s"%(Node,SolsAliasName)
            print(ss)
            os.system("%s > /dev/null 2>&1"%ss)

            ss="ssh %s ln -s %s %s"%(Node,SmoothSolName,SolsAliasName)
            print(ss)
            os.system("%s > /dev/null 2>&1"%ss)
        else:
            if os.path.islink(SolsAliasName):
                os.unlink(SolsAliasName)
            os.symlink(SmoothSolName,SolsAliasName)


if __name__=="__main__":
    # masterNode = MPI.Get_processor_name()
    testParallel()
    # with MPIPoolExecutor() as executor:
    #     f1=executor.submit(filterHost, "nancep10.obs-nancay.fr" ,ftest, (10,),{})
    #     filterHost("nancep11.obs-nancay.fr" ,ftest, (11,),{})
    #     f1.result()
