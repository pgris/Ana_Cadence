import cPickle as pkl
from astropy.table import Table
import sncosmo
import numpy as np
import matplotlib.pyplot as plt
from Throughputs import Throughputs
import astropy.units as u
import pickle as pkl
from astropy import (cosmology, units as u, constants as const)
from optparse import OptionParser
import math

def sigma_rand(gamma,mag,mag5):
    x=np.power(10,0.4*(mag-mag5))
    return (0.04+gamma)*x+gamma*np.power(x,2.)


parser = OptionParser()

parser.add_option("-N", "--nevts", type="int", default=10, help="filter [%default]")
parser.add_option("-m", "--model", type="string", default='salt2-extended', help="filter [%default]")
parser.add_option("-f", "--fieldtype", type="string", default="WFD", help="filter [%default]")
parser.add_option("-n", "--fieldnum", type="int", default=309, help="filter [%default]")

opts, args = parser.parse_args()

idlist=opts.fieldtype+'_'+str(opts.fieldnum)+'_'+str(opts.nevts)+'_'+opts.model
listname='List_'+idlist+'_minion_test.dat'

filenames=[]

opsim_release='minion_1016'
thedir='Sim_'+opsim_release

for line in open(listname, 'r').readlines():
    filenames.append(thedir+'/'+line.strip())

all_obs=[]

for filename in filenames:
    print 'opening',filename
    pkl_file = open(filename,'rb')
    objs = []
    while 1:
        try:
            objs.append(pkl.load(pkl_file))
        except EOFError:
            break
    all_obs.append(objs)

ratio={}
test={}
bands=['u','g','r','i','z','y']

for band in bands:
    ratio[band]=[]

for band in bands:
    test[band]=[]

snr_ratio_a=[]
snr_ratio_b=[]

photgain=2.3

telSeeing = 0.250 # design goal
opticalDesSeeing = 0.08
cameraSeeing = 0.300
scaleToNeff = 1.16
atmNeffFactor = 1.04
FWHM_Sys_Zenith = np.sqrt(telSeeing**2 + opticalDesSeeing**2 + cameraSeeing**2)
filterWave = {'u': 367.0, 'g': 482.5, 'r': 622.2, 'i': 754.5, 'z': 869.1, 'y': 971.0}


for oob in all_obs:
    for i,obj in enumerate(oob):
        if obj['observations'] is not None:
            #print 'there pal',len(obj['observations'])
            #print obj['observations']['err_flux']/obj['observations']['err_flux_opsim']
            print obj['observations'].dtype
            print obj['observations']
            """
            for i in range(len(obj['observations'])):
                val=obj['observations'][i]
                print val
                filtre=val['filter'][0]
                Filter_Wavelength_Correction = np.power(500.0 / filterWave[filtre], 0.3)
                Airmass_Correction = math.pow(val['airmass'],0.6)
                FWHM_Sys = FWHM_Sys_Zenith * Airmass_Correction
                FWHM_Atm = val['rawSeeing']* Filter_Wavelength_Correction * Airmass_Correction
                finSeeing = scaleToNeff * math.sqrt(np.power(FWHM_Sys,2) + atmNeffFactor * np.power(FWHM_Atm,2))
                print 'hello',finSeeing,val['finSeeing'],val['rawSeeing'],val['airmass']
                
                m5_recalc=val['dCm']+val['Cm']+0.5*(val['filtSkyBrightness']-21.)+2.5*np.log10(0.7/val['finSeeing'])-val['katm_opsim']*(val['airmass']-1.)+1.25*np.log10(val['visitExpTime']/30.)
                #print 'hello',val['m5_calc'],val['fiveSigmaDepth'],val['m5_calc']/val['fiveSigmaDepth']
                if val['flux'] >0:
                    test[val['filter'][0]].append((val['m5_calc']/val['fiveSigmaDepth'])[0])
                    #print 'hello',val['err_flux'],val['err_flux_opsim'],val['err_flux']/val['err_flux_opsim']
                    ratio[val['filter'][0]].append(val['err_flux'][0]/val['err_flux_opsim'][0])
                    snr_ratio_a.append(val['snr_SED'][0]/val['snr_m5_through'][0])
                    snr_ratio_b.append(val['snr_m5_through'][0]/val['snr_m5_opsim'][0])

            """
            break
               
print len(ratio)   


posband={}

posband['u']=(0,0)
posband['g']=(0,1)
posband['r']=(1,0)
posband['i']=(1,1)
posband['z']=(2,0)
posband['y']=(2,1)

figa, axa = plt.subplots(ncols=2, nrows=3, figsize=(10,9))

for band in bands: 
    axa[posband[band][0]][posband[band][1]].hist(ratio[band],bins=20,histtype='step')
    
figb, axb = plt.subplots(ncols=2, nrows=3, figsize=(10,9))

for band in bands:
    axb[posband[band][0]][posband[band][1]].hist(test[band],bins=20,histtype='step')

figc, axc = plt.subplots(ncols=2, nrows=1, figsize=(10,9))
axc[0].hist(snr_ratio_a,bins=20,histtype='step')
axc[1].hist(snr_ratio_b,bins=20,histtype='step')

plt.show()
