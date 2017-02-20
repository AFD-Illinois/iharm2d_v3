#!/usr/bin/python
#
# import
import numpy as np
import matplotlib.pyplot as plt
#

# open data file and get data 
data = np.loadtxt("dumps/dump010", skiprows=1)
#
x1  	= data[:,0]
x2  	= data[:,1]
rho 	= data[:,2]
u 	= data[:,3]
v1 	= data[:,4]
v2 	= data[:,5]
v3 	= data[:,6]
B1 	= data[:,7]
B2 	= data[:,8]
B3 	= data[:,9]
divb	= data[:,10]
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
#
gam = 2.0
scale = 100.
p = (gam - 1.)*u*scale*scale
v = uu1*scale
By = B3*scale
#
# read in athena solution
adata = np.loadtxt("brio.dat", skiprows=7)
ax    = adata[:,1]
arho  = adata[:,2]
aMx   = adata[:,3]
aMy   = adata[:,4]
avx   = aMx/arho
ap    = adata[:,6]
aBx   = adata[:,8]
aBy   = adata[:,9]
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
plt.subplot(221)
plt.xlim(0.0,1.0)
plt.ylim(0.0,1.1)
plt.plot(x1,rho)
#plt.plot(x1,rho,'s')
plt.plot(ax,arho)
plt.xlabel('$x_1$', weight='bold', fontsize=20)
plt.ylabel('$\\rho$', weight='bold', fontsize=20)

plt.subplot(222)
plt.xlim(0.0,1.0)
plt.ylim(-0.5,1.0)
plt.plot(x1,v)
#plt.plot(x1,v,'s')
plt.plot(ax,avx)
plt.xlabel('$x_1$', weight='bold', fontsize=20)
plt.ylabel('$v$', weight='bold', fontsize=20)

plt.subplot(223)
plt.xlim(0.0,1.0)
plt.ylim(0.0,1.1)
plt.plot(x1,p)
#plt.plot(x1,p,'s')
plt.plot(ax,ap)
plt.xlabel('$x_1$', weight='bold', fontsize=20)
plt.ylabel('$p$', weight='bold', fontsize=20)

plt.subplot(224)
plt.xlim(0.0,1.0)
plt.ylim(-1.1,1.1)
plt.plot(x1,By)
#plt.plot(x1,By,'s')
plt.plot(ax,aBy)
plt.xlabel('$x_1$', weight='bold', fontsize=20)
plt.ylabel('$B_y$', weight='bold', fontsize=20)

# and, as in German, the verb comes last
plt.show()
#
