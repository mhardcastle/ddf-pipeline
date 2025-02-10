import os
try:
    import DDFacet.Other.ModColor as MC
    def Str(*args,**kwargs):
        return MC.Str(*args,**kwargs)
except:
    def Str(ss,**kwargs):
        return ss


from subprocess import PIPE
import subprocess

def getVersion(Name):
    try:
        r=getVersionSubprocess(Name)
        if isinstance(r,list):
            r=" ".join(r)
            if "ModuleNotFoundError" in r:
                return "ModuleNotFound"
        return r.strip()
                
    except Exception as e:
        print('%s: failed with: %s'%(Name, str(e)))

def OS(s):
    if isinstance(s,list):
        s=(" ".join(s))
    os.system("%s > /tmp/OS.log 2>&1"%s)
    return open("/tmp/OS.log","r").readlines()
        
def getVersionSubprocess(Name):
    if Name=="losoto":
        command=["losoto", "--version"]
        #out = check_output(command)
        out=OS(command)
        if "not found" in out[0]: return "Not installed"
        return out[0].strip()
    elif Name=="wsclean":
        command=["wsclean","--version"]
        out=OS(command)
        #result = subprocess.run(command, capture_output=True, text=True)
        #stderr,stdout=result.stderr,result.stdout
        out="".join(out)
        if "not found" in out[0]: return "Not installed"
        
        v=(out.replace("\n"," ").split(" WSClean version ")[1].split("This software package")[0]).replace(" ","")
        return v
    elif Name=="DP3" or Name=="aoflagger":
        command=[Name,"--version"]
        #result = subprocess.run(command, capture_output=True, text=True)
        #stderr,stdout=result.stderr,result.stdout
        out=OS(command)
        if "not found" in out[0]: return "Not installed"
        out="".join(out)
        v=(out.strip().lower().split(Name.lower())[1].replace(" ",""))
        return v
    elif Name=="DDF":
        command=["%s.py"%Name,"--version"]
        #result = subprocess.run(command, capture_output=True, text=True)
        #stderr,stdout=result.stderr,result.stdout
        r=OS(command)
        if "not found" in r[0]: return "Not installed"
        # print("rr",r)
        # #v= stdout.split("\n")[0].split("DDFacet version is ")[1]
        # print(r[-1].split("DDFacet version is "))
        rr="Not installed"
        for l in r:
            if 'DDFacet version is' in l:
                rr=l.split("DDFacet version is ")[1].strip()
                break
        
        return rr
    elif Name=="kMS":
        command=["%s.py"%Name,"--version"]
        # result = subprocess.run(command, capture_output=True, text=True)
        # stderr,stdout=result.stderr,result.stdout
        # return stdout.split("\n")[-2]
        r=OS(command)
        if "not found" in r[0]: return "Not installed"
        #print("rr",r)
        #v= stdout.split("\n")[0].split("DDFacet version is ")[1]
        return r[-1].strip()
    elif Name=="DynSpecMS":
        command=["ms2dynspec.py","--version"]
        result = subprocess.run(command, capture_output=True, text=True)
        stderr,stdout=result.stderr,result.stdout
        return stdout.split("\n")[-2].split("version ")[-1]
    elif Name=="LOFARBeam":
        command=["python","-c",'"import lofar.stationresponse"']
        # print(command)
        # result = subprocess.run(command, capture_output=True, text=True)
        # stderr,stdout=result.stderr,result.stdout

        s=(" ".join(command))
        r=OS(s)
        if len(r)==0:
            r="Installed, no version available"
        return r
    elif Name=="nenupy":
        command=["python","-c",'"import nenupy,sys; print(str(nenupy.__version__),file=sys.stdout)"']
        s=(" ".join(command))
        r=OS(s)
        return r
    #o=os.popen(s).read().strip()
        
    #    return o
    elif Name=="lsmtool" or Name=="LofarStMan":
        command=["python","-c",'"import %s,sys"'%Name]
        s=(" ".join(command))
        #o=os.popen(s).read().strip()
        r=OS(s)
        if len(r)==0:
            r="Installed, no version available"
        return r
    elif Name=="drawMS":
        command=["drawMS.py","--version"]
        r=OS(command)
        if "not found" in r[0]: return "Not installed"
        # print("rr",r)
        # #v= stdout.split("\n")[0].split("DDFacet version is ")[1]
        # print(r[-1].split("DDFacet version is "))
        rr="Not installed"
        for l in r:
            if 'drawMS.py version' in l:
                rr=l.split("drawMS.py version ")[1].strip()
                break
        return rr
    elif Name=="lotss-query":
        command=["python","-c",'"import surveys_db"']
        s=(" ".join(command))
        r=OS(s)
        if len(r)==0:
            r="Installed, no version available"
        return r
    elif Name=="ddf-pipeline":
        command=["python","-c",'"from pipeline_version import version; print(version())"']
        s=(" ".join(command))
        r=OS(s)
        if len(r)==0:
            r="Installed, no version available"
        return r
    else:
        print(Name)

        
        
DicoExe={"losoto":{"Type":"exe",
                   "Name":"self"},
         "DP3":{"Type":"exe",
                "Name":"self"},
         "wsclean":{"Type":"exe",
                    "Name":"self"},
         "aoflagger":{"Type":"exe",
                      "Name":"self"},
         "DDF":{"Type":"exe",
                "Name":"DDF.py"},
         "kMS":{"Type":"exe",
                "Name":"kMS.py"},
         "DynSpecMS":{"Type":"exe",
                      "Name":"ms2dynspec.py"},
         "LOFARBeam":{"Type":"Package",
                      "Name":"lofar.stationresponse"},
         "nenupy":{"Type":"Package",
                   "Name":"nenupy"},
         "lsmtool":{"Type":"exe",
                   "Name":"self"},
         "LofarStMan":{"Type":"Package",
                       "Name":"LofarStMan"},
         "drawMS":{"Type":"exe",
                   "Name":"self"},
         "lotss-query":{"Type":"Package",
                        "Name":"surveys_db"},
         "ddf-pipeline":{"Type":"Package",
                         "Name":"run_full_field_reprocessing_pipeline"}
}

SHIFT=30

def Print_v():
    for Name in DicoExe.keys():

        # D=DicoExe[Name]
        
        version=getVersion(Name)#.rjust(30," ")
        NameS=Name#.rjust(20," ")
        #print(NameS,version)
        pLine([NameS,version],SHIFT)
    
    # if D["Type"]=="exe":
    #     exeName=D["Name"]
    #     if exeName=="self":
    #         exeName=Name
    #     Path=shutil.which(exeName)
    #     if Path is None:
    #         Status=("Not installed")
    #     else:
    #         Status=("Installed    ")
    # else:
    #     try:
    #         exec("import %s"%D["Name"])
    #         Status=("Installed    ")
    #     except:
    #         Status=("Not installed")
    # ss="%20s : %s"%(Name,Status)


W=90
l="%"+"%is"%(W-5)
ll="  |%s|"%(l)
def pLine(s,justify="center",length=None):
    if justify=="center":
        #print(ll%(str(s).center(W-5," ")))
        ls=len(s)
        if length is not None:
            ls=length
        
        S=" "*((W-ls-5)//2)+s+" "*((W-ls-5)//2)
        S0=" "*((W-ls-5)//2)+" "*ls+" "*((W-ls-5)//2)
        if len(S0)%2==0: S+=" "
        print("  |%s|"%S)
        
        #print(ll%(str(s).center(W-5," ")))
    elif isinstance(s,list):
        s0,s1=s
        #print(s0,s1)
        #print(s0.rjust(justify," "))
        ss=(s0.rjust(justify," ")+" : "+s1)
        #print(ss,len(ss))
        F=" "*(W-len(ss)-5)
        print("  |"+ss+F+"|")
        
Sep="="*(W-5)
print()
pLine(Sep)
s="Radio soft singularity image"
pLine(Str(s.upper(),Bold=1,col="green"),length=len(s))
#pLine("Radio soft singularity image".upper())
pLine("")
pLine(["To source an external dev dir"," source setDev.sh <DirName>"],SHIFT)
pLine(["for example"," source setDev.sh /data/$USER/DEV"],SHIFT)
pLine(Sep)
#pLine(["<Software>","<Versions>"],SHIFT)
#Print_v()
#pLine(Sep)
print()
print()
