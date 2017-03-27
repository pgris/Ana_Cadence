from Throughputs import Throughputs
import matplotlib.pyplot as plt

transmission=Throughputs(aerosol=False)

#transmission.Plot_Throughputs()

transmission_new=Throughputs(through_dir='NEW_THROUGH',atmos_dir='NEW_THROUGH',aerosol=False)

#transmission_new.Plot_Throughputs()

filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}

fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
for i,band in enumerate(['u','g','r','i','z','y']):
            #plt.plot(self.lsst_std[band].wavelen,self.lsst_std[band].sb,linestyle='-',color=self.filtercolors[band], label='%s - std' %(band))
            ax.plot(transmission.lsst_system[band].wavelen,transmission.lsst_system[band].sb,linestyle='-',color=filtercolors[band], label='%s - LSST' %(band))
            ax.plot(transmission_new.lsst_system[band].wavelen,transmission_new.lsst_system[band].sb,linestyle='--',color=filtercolors[band], label='%s - LSST_new' %(band))

ax.legend(loc=('upper right'), fontsize='smaller', fancybox=True, numpoints=1)
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Sb (0-1)')
fig.suptitle('System throughput (detector+lenses+mirror)')


toload=['detector', 'lens1', 'lens2', 'lens3', 'm1', 'm2', 'm3']
for filtre in ['u','g','r','i','z','y']:
    toload.append('filter_'+filtre)

figb, axb = plt.subplots(ncols=4, nrows=4, figsize=(10,9))

for j,system in enumerate(toload):
    if j < 4:
        k=0
    if j >= 4 and j < 8:
        k=1
    if j >= 8 and j < 12:
        k=2
    if j >= 12 and j < 16:
        k=3   

    axb[k][j%4].plot(transmission.lsst_telescope[system].wavelen,transmission.lsst_telescope[system].sb,linestyle='-', color='k',label=system)
    axb[k][j%4].plot(transmission_new.lsst_telescope[system].wavelen,transmission_new.lsst_telescope[system].sb,linestyle='-', color='r',label=system)
    axb[k][j%4].legend(loc='best',prop={'size':8})

figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))

axc.plot(transmission.atmos.wavelen, transmission.atmos.sb, 'k:', label='X =%.1f atmos' %(transmission.airmass),linestyle='--')
axc.plot(transmission_new.atmos.wavelen,transmission_new.atmos.sb, 'sr', label='X =%.1f atmos' %(transmission_new.airmass))
axc.set_xlabel('Wavelength (nm)')
axc.set_ylabel('Sb (0-1)')
figc.suptitle('Atmospheric transmission (no aerosol)')
plt.show()
