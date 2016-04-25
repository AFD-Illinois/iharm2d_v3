#!/usr/bin/python
#
# import
import numpy as np
import matplotlib.pyplot as plt
#
# define some functions
#
# plot single variable 
def plot_one_var(num, var, label):
	plt.subplot(3,1,num)
	size = var.size
	me = np.mean(var[size/2:size])
	si = np.std(var[size/2:size])
	plt.ylim(me - 3.*si, me + 3.*si)
	plt.plot(t,var)
	plt.xlabel('$t/M$', weight='bold')
	plt.ylabel(label, weight='bold')
#
# set font size
plt.rc('text', usetex=True) 
plt.rc('font', size = 16)
plt.rc('font', family='sans-serif')
#
data = np.loadtxt("ener.out")
#
t = data[:,0]
m = data[:,1]
e = data[:,2]
l = data[:,3]
dm = data[:,4]
de = data[:,5]
dl = data[:,6]
#
sde = de/(dm + 1.e-10)  # specific energy accreted
sdl = dl/(dm + 1.e-10)  # specific angular momentum accreted
#
#
plt.subplots_adjust(left=0.2,hspace=0.45)
#
plot_one_var(3, sdl, '$\dot{L}/\dot{M}$')
plot_one_var(2, sde, '$\dot{E}/\dot{M}$')
plot_one_var(1, dm, '$\dot{M}$')
#
plt.show()
#
#
