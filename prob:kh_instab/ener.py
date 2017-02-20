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
u1 = data[:,3]
u2 = data[:,4]
u3 = data[:,5]
#
plt.subplots_adjust(left=0.2,hspace=0.45)
#
plot_one_var(3, m, '$M$')
plot_one_var(2, e, '$E$')
plot_one_var(1, u3,'$U_3$')
#
plt.show()
#
#
