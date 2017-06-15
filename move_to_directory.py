import sys
import os

listname=sys.argv[1]
dirname=listname.replace('-list.txt','_dir')
os.mkdir(dirname)

s=[l.split()[0] for l in open(listname).readlines()]
for name in s:
    os.system('mv '+name+'*.png '+dirname)
    os.system('mv '+name+'*manifest*txt '+dirname)
