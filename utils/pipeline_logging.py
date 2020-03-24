from __future__ import print_function
import subprocess
import sys
import select
import datetime

def run_log(cmd,logfile,quiet=False):
    logfile = open(logfile, 'w')
    logfile.write('Running process with command: '+cmd+'\n')
    proc=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,universal_newlines=True)
    while True:
        try:
            select.select([proc.stdout],[],[proc.stdout])
        except select.error:
            pass
        line=proc.stdout.readline()
        if line=='':
            break
        if not quiet:
            sys.stdout.write(line)
        ts='{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
        logfile.write(ts+': '+line)
        logfile.flush()
    retval=proc.wait()
    logfile.write('Process terminated with return value %i\n' % retval)
    return retval

if __name__=='__main__':
    v=run_log(' '.join(sys.argv[1:]),'test-log.txt',quiet=False)
    print('Return value was',v)
