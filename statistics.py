import numpy as np
import cmath
import random 
from scipy import fft, ifft
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#statistics
value1 = np.loadtxt("ekinetic.dat",usecols=(0),unpack=True)
value2 = np.loadtxt("epotential.dat",usecols=(0),unpack=True)
value3 = np.loadtxt("etotal.dat",usecols=(0),unpack=True)
value4 = np.loadtxt("temperature.dat",usecols=(0),unpack=True)
value5 = np.loadtxt("pressure.dat",usecols=(0),unpack=True)
value6,rad = np.loadtxt("gofr.dat",usecols=(0,2),unpack=True)
value7,tdel7 = np.loadtxt("vacf.dat",usecols=(0,2),unpack=True)
value8,tdel8,qvec = np.loadtxt("ddcf.dat",usecols=(0,2,3),unpack=True)

#n is the number of experiments, n1=nblk,n6=nbins,n7=t_delay,n8=nq
n1 = len(value1)
n6 = int(len(value6)/n1)
n7 = int(len(value7)/n1)
n8 = int(len(value8)/n1/n7)

avg1 = np.mean(value1)
var1 = np.var(value1,ddof=1)
err1 = np.sqrt(var1/n1)

avg2 = np.mean(value2)
var2 = np.var(value2,ddof=1)
err2 = np.sqrt(var2/n1)

avg3 = np.mean(value3)
var3 = np.var(value3,ddof=1)
err3 = np.sqrt(var3/n1)

avg4 = np.mean(value4)
var4 = np.var(value4,ddof=1)
err4 = np.sqrt(var4/n1)

avg5 = np.mean(value5)
var5 = np.var(value5,ddof=1)
err5 = np.sqrt(var5/n1)

#nbins data points long
avg6=np.empty(n6)
var6=np.empty(n6)
err6=np.empty(n6)

#t_delay data points long
avg7=np.empty(n7)
var7=np.empty(n7)
err7=np.empty(n7)

#t_delay by nq data points
#leave it to the plot script to slice it up and make a heatmap
avg8=np.empty(n7*n8)
var8=np.empty(n7*n8)
err8=np.empty(n7*n8)

stat_stats=open('stat_statistics.dat','w')
stat_stats.write("Kinetic energy:                     K = {:f} +/- {:f} \n".format(avg1,err1))
stat_stats.write("Potential energy:                   V = {:f} +/- {:f} \n".format(avg2,err2))
stat_stats.write("Total energy:                       E = {:f} +/- {:f} \n".format(avg3,err3))
stat_stats.write("Temperature:                        T = {:f} +/- {:f} \n".format(avg4,err4))
stat_stats.write("Pressure:                           P = {:f} +/- {:f} \n".format(avg5,err5))
stat_stats.close()

gofr_stats=open('gofr_statistics.dat','w')
vacf_stats=open('vacf_statistics.dat','w')
ddcf_stats=open('ddcf_statistics.dat','w')

for i in range(n6):
  avg6[i]=np.mean(value6[i::n6])
  var6[i]=np.var(value6[i::n6],ddof=1)
  err6[i]=np.sqrt(var6[i]/n1)
  gofr_stats.write("{:f}\t\t\t{:f}\t\t\t{:f}\n".format(rad[i],avg6[i],err6[i]))

for i in range(n7):
  avg7[i]=np.mean(value7[i::n7])
  var7[i]=np.var(value7[i::n7],ddof=1)
  err7[i]=np.sqrt(var7[i]/n1)
  vacf_stats.write("{:f}\t\t\t{:f}\t\t\t{:f}\n".format(tdel7[i],avg7[i],err7[i]))

#need to make this a double loop arranging 
#actually should make it using reshape
for i in range(n7*n8):
  avg8[i]=np.mean(value8[i::n7*n8])
  var8[i]=np.var(value8[i::n7*n8],ddof=1)
  err8[i]=np.sqrt(var8[i]/n1)
  ddcf_stats.write("{:f}\t\t\t{:f}\t\t\t{:f}\t\t\t{:f}\n".format(tdel8[i],qvec[i],avg8[i],err8[i]))

gofr_stats.close()
vacf_stats.close()
ddcf_stats.close()



#calculate S(q,w), numerically integrate to calculate Fourier transform of F(q,t)
#number of omega values, arbitrary can be selected by user
#for GIFT run typically 600
#n9=nq
n9=10
dw=0.1
#time step, when complete this will be accessible from header file
dt=0.001

t,q,Fqt,err = np.loadtxt("ddcf_statistics.dat",usecols=(0,1,2,3),unpack=True)


#F(q,t) is the intermediate scattering function
#need to rearrange contents of Fqt to make it t_delay by nq
Fqt=np.reshape(Fqt,(n7,n8))


#from 0 to approx positive infinity, need enough data for approximation to be justified
sqw=np.empty((n8,n9),dtype=np.complex_)

print("t_delay = n7 = ",n7)
print("nq = n8 = ",n8)
print("nw = n9 = ",n9)






#fourier transform of g(r)
#Sq=fft(avg6)
#plt.plot(Sq)
#plt.show()





