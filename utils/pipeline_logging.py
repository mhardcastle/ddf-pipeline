#!/usr/bin/python

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
            ready = select.select([proc.stdout], [], [], 1.0)[0]
        except OSError:
            ready = []
        if ready:
            line = proc.stdout.readline()
            if line == '':
                break
            if not quiet:
                sys.stdout.write(line)
            ts = '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
            logfile.write(ts+': '+line)
            logfile.flush()
        elif proc.poll() is not None:
            break
    retval=proc.wait()
    logfile.write('Process terminated with return value %i\n' % retval)
    return retval

if __name__=='__main__':
    v=run_log(' '.join(sys.argv[2:]),sys.argv[1],quiet=False)
    print('Return value was',v)
