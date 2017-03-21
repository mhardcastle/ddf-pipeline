import subprocess
import os

def version():
    prevdir = os.getcwd()
    path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(path)
    try:
        result=subprocess.check_output('git describe --tags', shell=True).rstrip()
    except:
        result='unknown'
    os.chdir(prevdir)
    return result

if __name__=='__main__':
    print version()
