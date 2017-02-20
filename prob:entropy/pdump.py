#!/usr/bin/python
#
# CFG pdump.py
#
# import
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
#
# get arguments
# argument 1 is dump number
parser = argparse.ArgumentParser()
parser.add_argument('dumpno', type=int, nargs=1)
parsed = parser.parse_args()   #parses sys.argv by default

# construct filename
fname = "dumps/pdump%03d" % parsed.dumpno[0]

# open data file and get data 
data = np.loadtxt(fname)
#
# set font size
plt.rc('text', usetex=True) 
plt.rc('font', size = 16)
plt.rc('font', family='sans-serif')
#
plt.plot(data[:,0],data[:,1],'k.')
#
# labels
plt.xlabel('$x_1$', weight='bold', fontsize=20)
plt.ylabel('$x_2$', weight='bold', fontsize=20)
#
# and, as in German, the verb comes last
plt.show()
#
