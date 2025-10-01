def unpackMSList(mslist):
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
    return ListMS

class mpi_manager():
    def __init__(self,options_cfg,mslist):
        pass

    def givePrefixDDF(self):
        if mpirun_singularity:
            return self.givePrefixDDF_mpirun_singularity()
        elif mpirun:
            return self.givePrefixDDF_mpirun_singularity()
        elif srun_singularity:
            pass
                
    
    def givePrefixDDF_mpirun(self,use_singularity=True):
        if not options['mpi_ddfacet']: return ""
        # MPI case
        cwd = os.getcwd()
        LocDDF_exec_inContainer="/usr/local/src/DDFacet/DDFacet/"
        Loc_Container=options['mpi_Singularity_cmd'] #"singularity exec -B/data -B/home /home/cyril.tasse/DDFSingularity/ddf_dev_np1.22.4.mpi.sif"
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
