import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import h5py
from matplotlib.colors import LogNorm

plt.rcParams.update({'font.size':14})

filename=sys.argv[1]
f = h5py.File(filename, 'r')
kx = f['/Timestep_0000/kx'].value
ky = f['/Timestep_0000/ky'].value

tsub = np.arange(0,len(f.keys()),20)

#Initialise the arrays containing the centroids
encentvec = np.zeros([2,len(tsub)])
enstcentvec = np.zeros([2,len(tsub)])
zoncentvec = np.zeros([2,len(tsub)])
count = 0
for i in tsub:
    t = str(i)
    tn = t.rjust(4,'0')
    gr = "Timestep_" + tn + "/W_hat"
    whattemp=f[gr].value.view(np.complex)
    toten=0.0
    totenst=0.0
    totzon=0.0
    totrand=0.0
    cenen=0.0
    cenenst=0.0
    cenzon=0.0
    #Calculating the centroids for energy, enstrophy and zonostrophy. Taking care to not double the value for k=0
    for j in range(len(kx)):
        kxd=float(kx[j])
        for k in range(1,len(ky)):
            kyd=float(ky[k])
            enst = 2*(kxd*kxd+kyd*kyd)*(kxd*kxd+kyd*kyd)*np.abs(whattemp[j,k])**2
            en = 2*(kxd*kxd+kyd*kyd)*np.abs(whattemp[j,k])**2
            zon = 2*((kxd*kxd)/(kxd*kxd+kyd*kyd+10.**(-16))**3)*(kxd*kxd+5*kyd*kyd)*np.abs(whattemp[j,k])**2
            rand = 2*np.abs(whattemp[j,k])**2

            toten += en
            totenst += enst
            totzon += zon
            totrand += rand

            cenen += np.array([np.abs(kx[j]),ky[k]])*en
            cenenst += np.array([np.abs(kx[j]),ky[k]])*enst
            cenzon += np.array([np.abs(kx[j]),ky[k]])*zon
        
        k = 0
        kyd=float(ky[k])
        enst = (kxd*kxd+kyd*kyd)*(kxd*kxd+kyd*kyd)*np.abs(whattemp[j,k])**2
        en = (kxd*kxd+kyd*kyd)*np.abs(whattemp[j,k])**2
        zon = ((kxd*kxd)/(kxd*kxd+kyd*kyd+10.**(-16))**3)*(kxd*kxd+5*kyd*kyd)*np.abs(whattemp[j,k])**2
        rand = np.abs(whattemp[j,k])**2

        toten += en
        totenst += enst
        totzon += zon
        totrand += rand

        cenen += np.array([np.abs(kx[j]),ky[k]])*en
        cenenst += np.array([np.abs(kx[j]),ky[k]])*enst
        cenzon += np.array([np.abs(kx[j]),ky[k]])*zon


    print(i,toten,totenst,totzon,totrand)
    encentvec[:,count]=cenen/toten
    enstcentvec[:,count]=cenenst/totenst
    zoncentvec[:,count]=cenzon/totzon
    count+=1


#Adding extra lines on the plot for visualisation purposes.
k0=1
th=np.linspace(np.pi/4,np.pi/2,50)
x = np.sqrt(2)*np.cos(th)
y = np.sqrt(2)*np.sin(th)
x2 = np.linspace(0,1,500)
y2 = np.sqrt(((np.sqrt(k0**2 + k0**2)**3/k0)**(2./3.))*(x2)**(2./3.) - (x2)**(2))
x3 = np.linspace(1,1.57,5200)
y3 = np.sqrt(((np.sqrt(k0**2 + k0**2)**4/k0)**(2./4.))*(x3)**(2./4.) - (x3)**(2))


#Plotting the evolution of the centroids
enplot = plt.plot(encentvec[0,:]/20,encentvec[1,:]/20,color='r',marker='o', markersize=7, linewidth=2,markevery=2)
enstplot = plt.plot(enstcentvec[0,:]/20,enstcentvec[1,:]/20,color='b',marker='s', markersize=7, linewidth=2,markevery=2)
zonplot = plt.plot(zoncentvec[0,:]/20,zoncentvec[1,:]/20,color='g',marker='v', markersize=7, linewidth=2,markevery=2)
axes=plt.gca()
axes.legend(('energy','enstrophy','zonostrophy'),loc=4)
plt.plot(x,y,color='k')
plt.plot(x2,y2,color='k')
plt.plot(x3,y3,color='k')




# axes.legend((enplot,enstplot,zonplot),('energy','enstrophy','zonostrophy'),loc=4)
axes.set_xlim([0,2])
axes.set_ylim([0,2])
plt.xlabel(r'$k_x/k_{0x}$')
plt.ylabel(r'$k_y/k_{0y}$')
plt.show()






