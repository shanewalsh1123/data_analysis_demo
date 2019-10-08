#This code calculates the amount of enstrophy transferred to zonal scales.
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
# import scipy
import h5py
from matplotlib.colors import LogNorm

plt.rcParams.update({'font.size':14})

markerlist=["o", "s", "v", "^", "*"]
for ii in range(len(sys.argv)-1):
    filename=sys.argv[ii + 1]
    f = h5py.File(filename, 'r')
    kx = f['/Timestep_0000/kx'].value
    ky = f['/Timestep_0000/ky'].value

    tsub = np.arange(0,len(f.keys()),10)

    encentvec = np.zeros([2,len(tsub)])
    enstcentvec = np.zeros([2,len(tsub)])
    zoncentvec = np.zeros([2,len(tsub)])
    distenst = np.zeros(len(tsub))
    distengy = np.zeros(len(tsub))
    totensthighkvec = np.zeros(len(tsub))
    enzontransvec = np.zeros(len(tsub))
    count = 0
    for i in tsub:
        t = str(i)
        tn = t.rjust(4,'0')
        gr = "Timestep_" + tn + "/W_hat"
        whattemp=f[gr].value.view(np.complex)
        toten=0.0
        totenst=0.0
        totensthighk=0.0
        enzontrans=0.0
        totzon=0.0
        totrand=0.0
        cenen=0.0
        cenenst=0.0
        cenzon=0.0
        for j in range(len(kx)):
            for k in range(1,len(ky)):
                en = 2*(kx[j]*kx[j]+ky[k]*ky[k])*np.abs(whattemp[j,k])**2

                toten += en

                if((np.abs(kx[j])<2) and ( np.abs(ky[k])<(20*np.sqrt(2)))):
                    enzontrans += en
                
            
            k = 0
            en = (kx[j]*kx[j]+ky[k]*ky[k])*np.abs(whattemp[j,k])**2

            toten += en

            if((np.abs(kx[j])<2) and ( np.abs(ky[k])<(20*np.sqrt(2)))):
                enzontrans += en
        print(ii,i)
        enzontransvec[count] = enzontrans/toten
        count+=1

    tmax =9.12 
    tvec = np.linspace(0,9.12,len(enzontransvec))
    plt.plot(tvec, enzontransvec, marker=markerlist[ii],markevery=7,markersize=7)
axes=plt.gca()
axes.set_xlim([0,9.12])

plt.ylabel(r'$E_{zon}/E_{tot}$')
plt.xlabel(r'$t\sqrt{E}$')
plt.legend([r'$E_1$',r'$E_2$',r'$E_3$',r'$E_4$',r'$E_5$'], loc = 2)
plt.show()






