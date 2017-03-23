from lsst.sims.photUtils import SignalToNoise
from lsst.sims.photUtils import PhotometricParameters
from Throughputs import Throughputs
from lsst.sims.photUtils import Bandpass,Sed
import numpy as np
from astropy.table import vstack,Table
import matplotlib.pyplot as plt

filtre='r'
airmass=1.
FWHMeff=[0.8,1.0,1.2]
filtSkyBrightness=[19., 20., 21.]
#mag_SN=18
transmission=Throughputs()
transmission.Load_Atmosphere(airmass)

step=0.1

m5_tab=Table(names=('seeing','msky','m5'), dtype=('f8', 'f8','f8'))
snr_tab=Table(names=('seeing','msky','m5','mag','snr'), dtype=('f8', 'f8','f8','f8','f8'))

for seeing in FWHMeff:
    for sky in filtSkyBrightness:

        wavelen_min, wavelen_max, wavelen_step=transmission.lsst_system[filtre].getWavelenLimits(None,None,None)
        flatSed = Sed()
        flatSed.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
        flux0=np.power(10.,-0.4*sky)
        flatSed.multiplyFluxNorm(flux0)
        
        photParams = PhotometricParameters()
 
        m5_calc=SignalToNoise.calcM5(flatSed,transmission.lsst_atmos[filtre],transmission.lsst_system[filtre],photParams=photParams,FWHMeff=seeing)
        m5_tab.add_row((seeing,sky,m5_calc))
        for mag_SN in np.arange(18,m5_calc,step):

            snr_m5_through,gamma_through=SignalToNoise.calcSNR_m5(mag_SN,transmission.lsst_atmos[filtre],m5_calc,photParams)

            snr_tab.add_row((seeing,sky,m5_calc,mag_SN,snr_m5_through))
            print 'Res',seeing, sky, mag_SN, m5_calc,1./snr_m5_through

col=['b','r','g']
ls=['solid','dashed','dotted']
figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(8,5))
for i,seeing in enumerate(FWHMeff):
    sel=m5_tab[np.where(m5_tab['seeing']==seeing)]
    ll='seeing : '+str(seeing)
    axa.plot(sel['msky'],sel['m5'],linestyle=ls[i],color=col[i],label=ll)

axa.set_xlabel('msky')
axa.set_ylabel('m5')
axa.legend(loc='upper left',prop={'size': 10})
figa.suptitle('Filter : '+filtre)

figb, axb = plt.subplots(ncols=1, nrows=1, figsize=(8,5))

for j,seeing in enumerate(FWHMeff):
    for i,sky in enumerate(filtSkyBrightness):
        sel=snr_tab[np.where(np.logical_and(snr_tab['seeing']==seeing,snr_tab['msky']==sky))]
        ll='(seeing,msky)=('+str(seeing)+','+str(sky)+')'
        axb.plot(sel['mag'],1./sel['snr'],linestyle=ls[j],color=col[i],label=ll)
axb.set_xlabel('mag')
axb.set_ylabel('1./SNR')
axb.legend(loc='upper left',prop={'size':10})
figb.suptitle('Filter : '+filtre)

plt.show()
