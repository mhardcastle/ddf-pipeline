#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import os

mslist=sys.argv[1]
filenames=[l.rstrip() for l in open(mslist).readlines()]
for f in filenames:
    gufile=f+'/gu.npy'
    if os.path.isfile(gufile):
        gu=np.load(gufile)
        plt.plot(gu[0,:,0],label=f)

plt.legend(loc=0)
plt.show()
