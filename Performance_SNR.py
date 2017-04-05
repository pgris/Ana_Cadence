from lsst.sims.photUtils import SignalToNoise
from lsst.sims.photUtils import PhotometricParameters
from Throughputs import Throughputs
from lsst.sims.photUtils import Bandpass,Sed
import numpy as np
from astropy.table import vstack,Table
import matplotlib.pyplot as plt
from lsst.sims.photUtils import Bandpass,Sed
from Throughputs import Throughputs
from lsst.sims.photUtils import PhotometricParameters

class Get_References:
    
    def __init__(self,airmass=1.,aerosol=False):

        
        #self.transmission.Load_Atmosphere(airmass)

        self.filters=['u','g','r','i','z','y']

        self.paper={}
        self.paper['mbsky']={'u': 22.92 ,'g': 22.27 ,'r':21.20,'i':20.47,'z':19.59,'y':18.63}
        self.paper['Tb']={'u': 0.0379 ,'g': 0.1493 ,'r':0.1386,'i':0.1198,'z':0.0838,'y':0.0413}
        self.paper['Sigmab']={'u': 0.0574 ,'g': 0.1735 ,'r':0.1502,'i':0.1272,'z':0.0872,'y':0.0469}
        self.paper['Seeing']={'u': 0.77 ,'g': 0.73 ,'r':0.70,'i':0.67,'z':0.65,'y':0.63}
        self.paper['Cb']={'u': 421.4,'g': 691.4 ,'r':952.9,'i':1152,'z':1373,'y':1511}
        self.paper['Skyb']={'u': 85.07,'g':467.9 ,'r':1085.2,'i':1800.3,'z':2775.7,'y':3613.4}
        self.paper['mb_Z']={'u': 27.09,'g':28.58 ,'r':28.50,'i':28.34,'z':27.95,'y':27.18}
        self.paper['counts_mb_Z']={'u': 1.0,'g':1.0 ,'r':1.0,'i':1.0,'z':1.0,'y':1.0}
        self.paper['mb_stnd']={'u': 24.22,'g': 25.17 ,'r':24.74,'i':24.38,'z':23.80,'y':22.93}


        transmission_git=Throughputs(aerosol=aerosol)
        transmission_lse40=Throughputs(through_dir='NEW_THROUGH',atmos_dir='NEW_THROUGH',aerosol=aerosol)
        airmass=airmass
        transmission_git.Load_Atmosphere(airmass)
        transmission_lse40.Load_Atmosphere(airmass)

        self.throughput_git={}
        self.throughput_lse40={}

        for key in self.paper.keys():
            self.throughput_git[key]={}
            self.throughput_lse40[key]={}
            for filtre in self.filters:
                self.throughput_git[key][filtre]={}
                self.throughput_lse40[key][filtre]={}
               

        self.Calc_Inputs(transmission_git,self.throughput_git)
        self.Calc_Inputs(transmission_lse40,self.throughput_lse40)
        self.Calc_Sky(self.paper,self.throughput_lse40,transmission_lse40)
        self.Calc_Sky(self.paper,self.throughput_git,transmission_git)

    def Calc_Inputs(self,transmission, tofill):

        for filtre in self.filters:
            
            myup=transmission.darksky.calcInteg(transmission.lsst_system[filtre])
               
            tofill['Tb'][filtre]=self.Calc_Integ(transmission.lsst_atmos[filtre])
            Sigmab=self.Calc_Integ(transmission.lsst_system[filtre])
            tofill['Sigmab'][filtre]=Sigmab
                
            tofill['mbsky'][filtre]=-2.5*np.log10(myup/(3631.*Sigmab))
        
    def Calc_Sky(self,paper,infos,transmission):

        Diameter=6.5 #m
        Deltat=30 #s
        platescale=0.2 #arsec
        gain=2.3

        for filtre in self.filters:

            filtre_trans=transmission.lsst_system[filtre]
            wavelen_min, wavelen_max, wavelen_step=filtre_trans.getWavelenLimits(None,None,None)
            photParams = PhotometricParameters()
            #photParams._exptime=30.
            
            bandpass=Bandpass(wavelen=filtre_trans.wavelen, sb=filtre_trans.sb)

            infos['Skyb'][filtre]=5455*np.power(Diameter/6.5,2.)*np.power(Deltat/30.,2.)*np.power(platescale,2.)*np.power(10.,0.4*(25.-infos['mbsky'][filtre]))*infos['Sigmab'][filtre]
            
            Zb=181.8*np.power(Diameter/6.5,2.)*infos['Tb'][filtre]
            mbZ=25.+2.5*np.log10(Zb)
            flatSed = Sed()
            flatSed.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
            flux0=np.power(10.,-0.4*mbZ)
            flatSed.multiplyFluxNorm(flux0)
            counts = flatSed.calcADU(bandpass, photParams=photParams) #number of counts for exptime
            infos['mb_Z'][filtre]=mbZ
            infos['counts_mb_Z'][filtre]=counts/photParams.exptime
            
            flatSedb = Sed()
            flatSedb.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
            flux0b=np.power(10.,-0.4*infos['mbsky'][filtre])
            flatSedb.multiplyFluxNorm(flux0b)
            #FWHMeff = SignalToNoise.FWHMgeom2FWHMeff(paper['Seeing'][filtre])
            FWHMeff = paper['Seeing'][filtre]
            #m5_calc=SignalToNoise.calcM5(flatSedb,transmission.lsst_atmos[filtre],transmission.lsst_system[filtre],photParams=photParams,FWHMeff=FWHMeff)
            m5_calc=SignalToNoise.calcM5(flatSedb,transmission.lsst_atmos[filtre],transmission.lsst_system[filtre],photParams=photParams,FWHMeff=FWHMeff)
            infos['mb_stnd'][filtre]=m5_calc

    def Calc_Integ(self,bandpass):
        resu=0.
        dlam=0
        for i,wave in enumerate(bandpass.wavelen):
            if i < len(bandpass.wavelen)-1:
                dlam=bandpass.wavelen[i+1]-wave
                resu+=dlam*bandpass.sb[i]/wave
            #resu+=dlam*bandpass.sb[i]

        return resu

resu=Get_References()

round_nums={'mbsky':2,'Tb':4,'Sigmab':4,'Skyb':2,'mb_Z':2,'counts_mb_Z':1,'mb_stnd':2}

for val in ['mbsky','Tb','Sigmab','Skyb','mb_Z','counts_mb_Z','mb_stnd']:
    print val
    for filtre in ['u','g','r','i','z','y']:
        diffa=resu.paper[val][filtre]-resu.throughput_lse40[val][filtre]
        diffb=resu.paper[val][filtre]-resu.throughput_git[val][filtre]
        print np.round(resu.paper[val][filtre],round_nums[val]),np.round(resu.throughput_lse40[val][filtre],round_nums[val]),np.round(diffa,round_nums[val]),np.round(resu.throughput_git[val][filtre],round_nums[val]),np.round(diffb,round_nums[val])

"""
filtre='r'
airmass=1.2
FWHMeff=[0.8,1.0,1.2]
filtSkyBrightness=[19., 20., 21.]
#mag_SN=18
transmission=Throughputs(through_dir='NEW_THROUGH',atmos_dir='NEW_THROUGH',aerosol=False)
transmission.Load_Atmosphere(airmass)

step=0.01

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
            if 1./snr_m5_through < 100.:
                snr_tab.add_row((seeing,sky,m5_calc,mag_SN,snr_m5_through))
            #print 'Res',seeing, sky, mag_SN, m5_calc,1./snr_m5_through

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

mag_ref=28.

for filtre in ['u','g','r','i','z']:
    filtre_trans=transmission.lsst_atmos[filtre]
    
    wavelen_min, wavelen_max, wavelen_step=filtre_trans.getWavelenLimits(None,None,None)
    flatSed = Sed()
    flatSed.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
    # normalize the SED so that it has a magnitude equal to the desired m5
    bandpass=Bandpass(wavelen=filtre_trans.wavelen, sb=filtre_trans.sb)
    fNorm = flatSed.calcFluxNorm(mag_ref, bandpass)
    flatSed.multiplyFluxNorm(fNorm)
    counts = flatSed.calcADU(bandpass, photParams=photParams)
    factor=counts*photParams.gain
    #factor=counts
    #print 'counts',filtre,counts,factor
    flatSedb= Sed()
    flatSedb.setFlatSED(wavelen_min, wavelen_max, wavelen_step)
    mag_refb=mag_ref+2.5*np.log10(factor)
    #print 'ZP',filtre,mag_refb
    fNormb = flatSedb.calcFluxNorm(mag_refb, bandpass)
    flatSedb.multiplyFluxNorm(fNormb)
    countsb = flatSedb.calcADU(bandpass, photParams=photParams)
    print 'ZP',filtre,mag_refb,countsb,countsb*photParams.gain
    
plt.show()
"""
