#!/usr/bin/python
#
# first argument is dump number.
# you must edit the file to change which variable is plotted!
#
# CFG dump.py, first go.
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
fname = "dumps/dump%03d" % parsed.dumpno[0]

# open data file and get data 
data = np.loadtxt(fname)
#
# first get line 0 header info
t       = data[0,0]
n1      = data[0,1]
n2      = data[0,2]
startx1 = data[0,3]
startx2 = data[0,4]
dx1     = data[0,5]
dx2     = data[0,6]
tf 	    = data[0,7]
nstep 	= data[0,8]
gam     = data[0,9]
cour 	= data[0,10]
DTd 	= data[0,11]  # dump freq
DTl 	= data[0,12]  # log freq
DTi 	= data[0,13]  # imag freq
DTr 	= data[0,14]  # restart freq
dump_cnt = data[0,15]  
image_cnt = data[0,16]  
rdump_cnt = data[0,17]  
dt  	= data[0,18]
#
# now get the rest of the variables 
data = np.loadtxt(fname, skiprows=1)
#
x1 	= data[:,0]
x2 	= data[:,1]
rho	= data[:,2]
u 	= data[:,3]
v1 	= data[:,4]
v2 	= data[:,5]
v3 	= data[:,6]
B1 	= data[:,7]
B2 	= data[:,8]
B3 	= data[:,9]
divb= data[:,10]
uu0	= data[:,11]  # second u means "up" index
uu1	= data[:,12]
uu2	= data[:,13]
uu3	= data[:,14]
ud0	= data[:,15]  # second d means "down" index
ud1	= data[:,16]
ud2	= data[:,17]
ud3	= data[:,18]
bu0	= data[:,19]  # u means "up" index
bu1	= data[:,20]
bu2	= data[:,21]
bu3	= data[:,22]
bd0	= data[:,23]  # d means "down" index
bd1	= data[:,24]
bd2	= data[:,25]
bd3	= data[:,26]
v1m	= data[:,27]
v1p	= data[:,28]
v2m	= data[:,29]
v2p	= data[:,30]
gdet= data[:,31]
ju0	= data[:,32]
ju1	= data[:,33]
ju2	= data[:,34]
ju3	= data[:,35]
jd0	= data[:,36]
jd1	= data[:,37]
jd2	= data[:,38]
jd3	= data[:,39]
#
# done!
#
# find derived quantities
jsq = ju0*jd0 + ju1*jd1 + ju2*jd2 + ju3*jd3
jdotu = ju0*ud0 + ju1*ud1 + ju2*ud2 + ju3*ud3
Jsq = jsq + jdotu*jdotu
gJsq = gdet*Jsq
#
p = (gam - 1.)*u
lp = np.log10(p)
K = p*rho**(-gam)
lK = np.log10(K)
EF = rho + gam*u
bsq = bu0*bd0 + bu1*bd1 + bu2*bd2 + bu3*bd3
EE = bsq + EF
va2 = bsq/EE
cs2 = gam*(gam - 1.)*u/EF
cms2 = cs2 + va2 - cs2*va2
T = (p/rho)*918.059
lT = np.log10(T)
ptot = p + 0.5*bsq
lbsq = np.log10(bsq + 1.e-20)
lrho = np.log10(rho)
ldivb = np.log10(abs(divb) + 1.e-20)
ibeta = 0.5*bsq/p
libeta = np.log10(ibeta + 1.e-20)
sigma = 0.5*bsq/rho
lsigma = np.log10(sigma + 1.e-20)
#
#
T00 = (rho+p+u+bsq)*uu0*ud0 + (p + 0.5*bsq) - bu0*bd0
T01 = (rho+p+u+bsq)*uu0*ud1 - bu0*bd1
T02 = (rho+p+u+bsq)*uu0*ud2 - bu0*bd2
T03 = (rho+p+u+bsq)*uu0*ud3 - bu0*bd3
#
print('done with preliminaries')
#
# plotting details
#
# set font size
plt.rc('text', usetex=True) 
plt.rc('font', size = 16)
plt.rc('font', family='sans-serif')
#
# plot the variable
# reshape the array for input to contour
plotted_variable = v2
pvar = np.reshape(plotted_variable, (n1, n2)) 
X = np.reshape(x1, (n1, n2)) 
Y = np.reshape(x2, (n1, n2)) 
# prevent dashed negative contours
#plt.rcParams['contour.negative_linestyle'] = 'solid'
#print len(x1var), len(x2var)
#plt.contour(X, Y, pvar, 20, colors = 'k')
#plt.pcolormesh(X, Y, pvar, cmap='afmhot', vmin=-0.01, vmax=0.01)
plt.pcolormesh(X, Y, pvar, cmap='afmhot')
# labels
plt.xlabel('$x_1$', weight='bold', fontsize=20)
plt.ylabel('$x_2$', weight='bold', fontsize=20)
# and, as in German, the verb comes last
plt.show()
#
