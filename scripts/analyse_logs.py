#!/usr/bin/env python

from __future__ import print_function
from builtins import range
from datetime import datetime
import glob
import os
import numpy as np
import matplotlib.pyplot as plt

def get_time_range(file):
    lines=open(file).readlines()
    start=datetime.strptime(lines[1][:19], '%Y-%m-%d %H:%M:%S')
    end=datetime.strptime(lines[-2][:19], '%Y-%m-%d %H:%M:%S')
    return (end-start).total_seconds()

os.chdir('logs')
g=glob.glob('*')
times=[]
for f in g:
    dt=get_time_range(f)
    times.append(dt)

labels=['Miscellaneous','Dynspec','Clipcal','Shift','Wide-field KillMS','KillMS phase 60sb','KillMS amp/phase 60sb', 'KillMS amp/phase full', 'KillMS amp/phase full 2', 'Wide-field DDF dirin', 'Predict wide-field', 'DDF dirin', 'DDF phase 60sb', 'DDF amp-phase 60sb', 'DDF band images','DDF full', 'DDF full 2', 'DDF full low QU','DDF full vlow QU','DDF full low V','DDF full low', 'DDF bootstrap', 'DDF bootstrap single-band','KillMS DIS0','KillMS DDS0','KillMS DIS1','KillMS DDS1','KillMS DIS2','KillMS DDS2','KillMS DDS3','DDF full DI','DDF predict']
fragments=['***','dynspec','ClipCal','shift','wide_killms_p1','killms_p1','killms_ap1','killms_f_ap1','killms_f_ap2','wide_image_dirin', 'wide_image_phase1_predict', 'image_dirin','image_phase1','image_ampphase1','NS_Band','image_full_ampphase1','image_full_ampphase2', 'image_full_low_QU', 'image_full_vlow_QU', 'image_full_high_stokesV','image_full_low','image_bootstrap','image_low','_DIS0','_DDS0','_DIS1','_DDS1','_DIS2','_DDS2','_DDS3','image_full_ampphase_di','DDF-Predict']
sums=np.zeros(len(labels))
print(len(labels),len(fragments))

# classify each file
for j,f in enumerate(g):
    label=0
    for i in range(len(labels)):
        if fragments[i] in f:
            label=i
            break
    
    print(f,times[j],labels[label],fragments[label])
    sums[label]+=times[j]

print('------\n')

pl=[]
ps=[]
print('---------------------------------\nOperation                 Time(s)\n---------------------------------')
for i in range(len(labels)):
    if sums[i]:
        print("%-25s %7.0f" % (labels[i],sums[i]))
        pl.append(labels[i])
        ps.append(sums[i])

tt=np.sum(sums)
print('---------------------------------\n\nTotal time %.0f seconds -- %.2f days' % (tt,tt/86400.0))

cmap = plt.cm.spectral
colors = cmap(np.linspace(0.1, 1., len(pl)))

fig1, ax1 = plt.subplots()
ax1.pie(ps,labels=pl,autopct='%1.1f%%',colors=colors)
ax1.axis('equal')
plt.show()
