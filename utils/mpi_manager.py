from mpi4py.futures import MPIPoolExecutor
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor, MPIPoolExecutor

size = MPI.COMM_WORLD.size




class MSSet():
    def __init__(self,mslist):
        self.file_nodes_mslist=mslist
        # for MPI use
        ListMS=[]
        mslist=[s.strip() for s in open(mslist).readlines()]

        for iMS,sMS in enumerate(mslist):
            DicoMS={}
            if ":" in sMS:
                Node,MSName=sMS.split(":")
                DicoMS["Node"]=Node
                DicoMS["MSName"]=MSName
            else:
                DicoMS["Node"]=None
                DicoMS["MSName"]=sMS
            ListMS.append(DicoMS)
    
        DicoNodes={}
        Lmslist=[]
        for MS in ListMS:
            ThisNode=MS["Node"]
            L=DicoNodes.get(ThisNode,[])
            L.append(MS["MSName"])
            DicoNodes[ThisNode]=L
            Lmslist.append(MS["MSName"])
            
        self.ListDicoMS=ListMS
        self.DicoNodes2ListMS=DicoNodes
        self.ListMS=Lmslist
        
        
        
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
    
    ListJobs=[["nancep10.obs-nancay.fr",os.system,("CleanSHM.py",), {}],
              ["nancep11.obs-nancay.fr",os.system,("CleanSHM.py",), {}],
              ]

    import DDFacet.CleanSHM
    ListJobs=[["nancep10.obs-nancay.fr",DDFacet.CleanSHM.driver,(), {}],
              ["nancep11.obs-nancay.fr",DDFacet.CleanSHM.driver,(), {}],
              ]
    
    callParallel(ListJobs)
    
def filterHost(RunOnHost,func, args, kwargs):
    host = MPI.Get_processor_name()
    if RunOnHost != host: return None
    return func(*args,**kwargs)
    

def callParallel(ListJobs):
    masterNode = MPI.Get_processor_name()
    print(masterNode)
    LJobMasterNode=[]
    with MPIPoolExecutor() as executor:
        Lres=[]
        print("Pool",masterNode)
        
        for Job in ListJobs:
            host,func,args,kwargs=Job
            if host == masterNode:
                LJobMasterNode.append([func,args,kwargs])
            else:
                f1=executor.submit(filterHost, host ,func, args, kwargs)
                Lres.append(f1)
                
        Lres0=[]
        for func,args,kwargs in LJobMasterNode:
            Lres0.append(func(*args,**kwargs))
                
        for res in Lres:
            try:
                res.result()
                Lres0.append(res.result())
            except Exception as e:
                print(f"e: {e}")
    return Lres0
        
class mpi_manager():
    def __init__(self,options_cfg,MSSet):
        self.options=options_cfg
        self.MSSet=MSSet
        self.ListDicoMS=MSSet.ListDicoMS
        self.ListNodesBeingUsed=sorted(list(set([MS.get("Node",None) for MS in self.ListDicoMS])))
        if self.ListNodesBeingUsed==[None]:
            self.ListNodesBeingUsed=None
        self.ddf_nproc = int(self.options.get('ddf_nproc', 1))
        self.UseMPI=False
        if self.ddf_nproc > 1 or self.ListNodesBeingUsed:
            self.UseMPI=True
        self.createRemoteLocal_mslist()
        
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
        for Node in self.MSSet.DicoNodes2ListMS.keys():
            Listms=self.MSSet.DicoNodes2ListMS[Node]
            FName="local_%s_mslist.txt"%Node
            f=open(FName,"w")
            for msname in Listms:
                f.write("%s\n"%msname)
            f.close()
            self.DicoNode2mslist[Node]=FName

    
    def givePrefixDDF_mpirun(self,use_singularity=True):

        # MPI case
        cwd = os.getcwd()
        LocDDF_exec_inContainer="/usr/local/src/DDFacet/DDFacet/"
        Loc_Container=self.options['mpi_Singularity_cmd'] # "singularity exec -B/data -B/home /home/cyril.tasse/DDFSingularity/ddf_dev_np1.22.4.mpi.sif"
        
        try:
            nNodes=int(options['mpi_ddfacet_nodes'])
            sNodes="-np %i"%nNodes
        except:
            LNodes=str(options['mpi_ddfacet_nodes']).split(",")
            nNodes=len(LNodes)
            HostName=socket.gethostname()
            LName=[HostName]
            for ThisNodeName in LNodes:
                if HostName in ThisNodeName: continue
                LName.append(ThisNodeName)
            LNodes=["%s:1"%nameNode for nameNode in LName]
            sNodes=",".join(LNodes)
            sNodes="-np %i --host %s"%(nNodes,sNodes)
        PrefixMPI="mpirun %s -wdir %s %s python %s"%(sNodes,cwd,Loc_Container,LocDDF_exec_inContainer)
        return PrefixMPI


    def givePrefixRunCommandFork(self):
        if not self.mpi_enable: return ""
        run_commands = [
            'srun',
            '--nodes=1',
            '--ntasks=1',
            '--cpus-per-task=40',
            '--time=05:00:00',
            '--hint=nomultithread',
            'singularity',
            'run',
            '-B',
            '/lustre/fsn1/projects/rech/doz/udd71uc/small/small/:/linkhome/rech/genrnu01/$USER/',
            '/lustre/fsn1/singularity/images/udd71uc/ddf.sif'
        ]
        commands_list = s.split()
        run_commands.extend(commands_list)
        s = ' '.join(run_commands)
        pass

    def scpScatter(self,FileName):
        pass

    def scpGather(self,FileName):
        pass    

def ftest(x):
    print(x)
    
if __name__=="__main__":
    # masterNode = MPI.Get_processor_name()
    testParallel()
    # with MPIPoolExecutor() as executor:
    #     f1=executor.submit(filterHost, "nancep10.obs-nancay.fr" ,ftest, (10,),{})
    #     filterHost("nancep11.obs-nancay.fr" ,ftest, (11,),{})
    #     f1.result()
