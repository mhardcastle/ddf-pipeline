import subprocess
import os

def version():
    path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(path)
    try:
        result=subprocess.check_output('git describe --tags', shell=True).rstrip()
    except:
        result='unknown'
    return result

if __name__=='__main__':
    print version()
