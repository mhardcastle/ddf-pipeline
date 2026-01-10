try:
    from mpi4py.futures import MPIPoolExecutor
    from mpi4py import MPI
    from mpi4py.futures import MPICommExecutor, MPIPoolExecutor
    MPIsize = MPI.COMM_WORLD.size
except:
    MPIsize = 0


import itertools
import os
from DDFacet.Other import ModColor
import pprint
from pathlib import Path

WIDTH_PROMPT=90

NFS_MODE=True
NFS_MODE=False

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
            else:
                node,msname=None,sMS
            l=nodes2ms.get(node,[])
            l.append(msname)
            nodes2ms[node]=l
            self.DicoMSName2Node[msname]=node
            
        if (None in nodes2ms.keys()) and len(nodes2ms) == 1 and MPIsize>1:
            # get all node names because None have been specified
            fns=[]
            with MPIPoolExecutor() as executor:
                fns=[]
                fns.append(executor.submit(get_node_name))

            # # CT: This was breaking stuff in nancap - not sure what the logic is
            # self.ListNodesBeingUsed=[f"{MPI.Get_processor_name()}@0"]
            self.ListNodesBeingUsed=[MPI.Get_processor_name()]
            for f in fns:
                self.ListNodesBeingUsed.append(f.result())
                
            mslist=list(zip(itertools.cycle(self.ListNodesBeingUsed), nodes2ms[None]))
            del nodes2ms[None]
            
            for node,ms in mslist:
                l=nodes2ms.get(node,[])
                l.append(ms)
                nodes2ms[node] = l

        self.DicoNodes2ListMS=nodes2ms
        self.ListNodesBeingUsed = list(nodes2ms.keys())
        
        
        
        print(ModColor.Str(" Data set distribution ".center(WIDTH_PROMPT,"="),col="blue"))
        print(ModColor.Str((" %s "%(str(self.ListNodesBeingUsed))).center(WIDTH_PROMPT,"="),col="blue"))
        pprint.pp(self.DicoNodes2ListMS)

        
def testFunc(*args,**kwargs):
    host = MPI.Get_processor_name()
    print(host,args,kwargs)
    
import os
def testParallel():
    ListJobs=[["nancep10.obs-nancay.fr",testFunc,(5,),{"e":6}],
              ["nancep10.obs-nancay.fr",testFunc,(9,),{"f":6}],
              ["nancep11.obs-nancay.fr",testFunc,(77,),{"g":88}],
              ]
    
    ListJobs=[["cw10055",os.system,("CleanSHM.py",), {}],
              ["cw10057",os.system,("CleanSHM.py",), {}],
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
    host = MPI.Get_processor_name()
    print("FilterHost : %s"%host)
    #localrank = os.environ.get("SLURM_LOCALID", "0")
    res=[]
    for RunOnHost, func, args, kwargs in jobs:
        #if RunOnHost == f"{host}@{localrank}":
        if RunOnHost == f"{host}":
            print("  [exec] [ME=%s][TARGET=%s]: %s(%s,%s)"%(host,RunOnHost,str(func),str(args),str(kwargs)))
            res.append(func(*args,**kwargs))
        else:
            # print("  [skip] [ME=%s][TARGET=%s]: %s(%s,%s)"%(host,RunOnHost,str(func),str(args),str(kwargs)))
            pass
    return res
    

def callParallel(ListJobs):
    masterNode = MPI.Get_processor_name()
    print()
    print("".center(WIDTH_PROMPT,"="))
    print(ModColor.Str(" CALL PARALLEL ".center(WIDTH_PROMPT,"="),col="blue"))
    print(ModColor.Str((" Master node: %s "%masterNode).center(WIDTH_PROMPT,"="),col="blue"))
    localrank = os.environ.get("SLURM_LOCALID", "0")
    LJobMasterNode=[]
    Lres0=[]
    with MPIPoolExecutor() as executor:
    #with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
        #mpi_comm_executor = MPICommExecutor(MPI.COMM_WORLD, root=0)
        Lres=[]
        for worker in range(MPI.COMM_WORLD.size-1):
            # submit remote jobs
            f1=executor.submit(filterHost, ListJobs)
            Lres.append(f1)
        
        # run local jobs
        for host,func,args,kwargs in ListJobs:
            #if host == f"{masterNode}@{localrank}":
            if host == masterNode:
                print("  [local %s]: %s(%s,%s)"%(host,str(func),str(args),str(kwargs)))
                Lres0.append(func(*args,**kwargs))
                
        for res in Lres:
            try:
                #res.result()
                Lres0.extend(res.result())
            except Exception as e:
                print(f"e: {e}")
        #MPI.COMM_WORLD.Barrier()
    print(f"callParallel: ok ({Lres0})")
    # print(" END callParallel ".center(WIDTH_PROMPT,"="))
    # print("".center(WIDTH_PROMPT,"="))
    print()
    return Lres0
        
# import itertools
# a= [7, 8, 4, 5, 9, 10]
# b= [1, 5, 6]
# 
# # Zipping using cycle
# res= list(zip(a, itertools.cycle(b)))
# print(res)
def get_node_name():
    # localrank = os.environ.get("SLURM_LOCALID", "0")
    # return f"{MPI.Get_processor_name()}@{localrank}"
    return MPI.Get_processor_name()

class mpi_manager():
    def __init__(self,options_cfg,MSSet, FullMSSet):
        self.options=options_cfg
        self.MSSet=MSSet
        self.FullMSSet=FullMSSet
        self.ListNodesBeingUsed=FullMSSet.ListNodesBeingUsed
        self.DicoNodes2WorkDir={}
        self.WorkDir=os.getcwd()
        self.MainHost = MPI.Get_processor_name()
                
        self.ddf_nproc = int(self.options.get('ddf_nproc', 1))
        self.UseMPI=False
        self.MPIsize=MPIsize
        if MPIsize>1 and (self.ddf_nproc > 1 or self.ListNodesBeingUsed):
            self.UseMPI=True

        # Scatter mslist and big-mslist.txt
        self.scpScatter(MSSet.file_nodes_mslist)
        self.scpScatter(FullMSSet.file_nodes_mslist)
        self.createRemoteLocal_mslist()
        self.createRemoteLocal_fullmslist()

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
        if NodeDest=="all":
            for Node in self.ListNodesBeingUsed:
                if Node==self.MainHost:
                    continue
                ss="scp -r %s %s:%s"%(FileName,Node,self.WorkDir)
                print("[Scatter] %s"%ss)
                os.system("%s > /dev/null 2>&1"%ss)
                os.system("%s > /dev/null 2>&1"%ss)
        else:
            if NodeDest==self.MainHost: return
            ss="scp -r %s %s:%s"%(FileName,NodeDest,self.WorkDir)            
            print("[Scatter] %s"%ss)
            os.system("%s > /dev/null 2>&1"%ss)
            #os.system("%s"%ss)

    def scpGatherSolutions(self,SolName,DestDir="",NodeSource="all"):
        if not self.UseMPI: return
        SolsDir=self.options["SolsDir"]
        AbsSolsDir=os.path.abspath(SolsDir)
        for Node in self.ListNodesBeingUsed:
            if Node==self.MainHost: return
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
        SolsDir=self.options["SolsDir"]
        AbsSolsDir=os.path.abspath(SolsDir)
        Node=self.FullMSSet.DicoMSName2Node[MSName]
        
        MSName = Path(MSName).name # if MSName is given with full path
        os.system("mkdir -p %s/%s"%(SolsDir,MSName))
        if Node!=self.MainHost:
            ss="scp -r %s %s:%s"%(SmoothSolName,Node,self.WorkDir)
            print("[Scatter Sols] %s"%ss)
            #os.system("%s &>/dev/null"%ss)
            os.system("%s > /dev/null 2>&1"%ss)
        
        ss="ssh %s rm %s"%(Node,SolsAliasName)
        print(ss)
        os.system("%s > /dev/null 2>&1"%ss)
        
        ss="ssh %s ln -s %s %s"%(Node,SmoothSolName,SolsAliasName)
        print(ss)
        os.system("%s > /dev/null 2>&1"%ss)
                

    # def createLink(self):
    #     os.symlink(os.path.abspath('%s_%.2f_smoothed.npz'%(ddsols,start_time)),symsolname)

                

def ftest(x):
    print(x)
    
if __name__=="__main__":
    # masterNode = MPI.Get_processor_name()
    testParallel()
    # with MPIPoolExecutor() as executor:
    #     f1=executor.submit(filterHost, "nancep10.obs-nancay.fr" ,ftest, (10,),{})
    #     filterHost("nancep11.obs-nancay.fr" ,ftest, (11,),{})
    #     f1.result()
