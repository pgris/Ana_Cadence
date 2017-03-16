import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
import astropy.units as u
import math
import healpy as hp
from astropy.table import Table
from optparse import OptionParser
from Parameters import parameters


def Get_median_finSeeing(observations):
    
    res=[]
    params=parameters()
    for obs in observations:
        filtre=obs['filter'][0]
        seeing=obs['rawSeeing']
        airmass=obs['airmass']
        Filter_Wavelength_Correction = np.power(500.0 / params.filterWave[filtre], 0.3)
        Airmass_Correction = math.pow(obs['airmass'],0.6)
        FWHM_Sys = params.FWHM_Sys_Zenith * Airmass_Correction
        FWHM_Atm = seeing * Filter_Wavelength_Correction * Airmass_Correction
        finSeeing = params.scaleToNeff * math.sqrt(np.power(FWHM_Sys,2) + params.atmNeffFactor * np.power(FWHM_Atm,2))
        res.append(finSeeing)
    return np.median(res)

def Load_Fields(fieldname):

    theres=[] 

    filename='fieldIDs_minion_1016_'+fieldname+'.txt'
    sfile=open(filename, 'r')

    for line in sfile.readlines():
            #print 'hello line',line,line.count('NONIA:')
        if line.count('fieldIds') > 0:
            theres.append(int(line.split(' ')[1].strip()))

    return theres

inputdir='/sps/lsst/data/dev/pgris/Obs_minion_1016'

prefix='Observations'

fieldtypes=['WFD','GalacticPlane','SouthCelestialPole-18','NorthEclipticSpur-18c','DDF']

#fieldtypes=['WFD','GalacticPlane','SouthCelestialPole-18','DDF']

toprocess={}
thedict={}

for typ in fieldtypes:
    toprocess[typ]=Load_Fields(typ)
    thedict[typ]={}

for key,vals in toprocess.items():
    for val in vals:
        keyb=key
        if key == 'DDF':
            keyb='DD'
        name=prefix+'_'+keyb+'_'+str(val)+'.pkl'
        pkl_file = open(inputdir+'/'+name,'rb')
        thedict[key][val]=pkl.load(pkl_file)['dataSlice']

fieldRA={}
fieldDec={}
nobs={}

for typ in fieldtypes:
    fieldRA[typ]=[]
    fieldDec[typ]=[]
    nobs[typ]=0

bands=['u','g','r','i','z','y']
for key,num in thedict.items():
    #print key,num
    for keyb,val in num.items():
        #print val['fieldRA'][0],val['fieldDec'][0],np.rad2deg(val['fieldRA'][0]),np.rad2deg(val['fieldDec'][0])
        fieldRA[key].append(val['fieldRA'][0])
        fieldDec[key].append(val['fieldDec'][0])
        for band in bands:
            sel=val[np.where(val['filter']==band)]
            #print band,np.median(sel['fiveSigmaDepth']),Get_median_finSeeing(sel)
        nobs[key]+=len(val)

ntot_obs=0
for typ in fieldtypes:
    ntot_obs+=nobs[typ]
print 'total obs:',ntot_obs

pourcent={}
for typ in fieldtypes:
    print 'obs',typ,nobs[typ],float(nobs[typ])/float(ntot_obs)
    pourcent[typ]=float(nobs[typ])/float(ntot_obs)

tot_label=[]

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111,projection="aitoff")
#ax = fig.add_subplot(111)

col={}
col['DDF']='r'
col['WFD']='b'
col['GalacticPlane']='k'
col['SouthCelestialPole-18']='g'
col['NorthEclipticSpur-18c']='y'

for typ in fieldtypes:
    # As next step, those coordinates are transformed into an astropy.coordinates
    # astropy.coordinates.SkyCoord object.
    c = SkyCoord(ra=fieldRA[typ]* u.radian, dec=fieldDec[typ]* u.radian, frame='icrs')

    ra_rad = c.ra.wrap_at(180 * u.deg).radian
    dec_rad = c.dec.radian

    """
    ra = coord.Angle(fieldRA[typ],unit=u.radian)
    ra = ra.wrap_at(180*u.degree)
    dec = coord.Angle(fieldDec[typ],unit=u.radian)
    """
    thelabel=typ+' - '+str(round(100.*pourcent[typ],1))+'%'
    tot_label.append(ax.scatter(ra_rad, dec_rad,color=col[typ],label=thelabel))
    #ax.set_xlim(360., 0.)
labs = [l.get_label() for l in tot_label]
ax.legend(tot_label, labs, ncol=2,loc='lower right',prop={'size':5},frameon=False)
ax.legend(bbox_to_anchor=(0.5, -0.1), loc=2, borderaxespad=0.,fontsize=10.)
ax.grid(True)




"""
ra_lsst=coord.Angle(-30.24287,unit=u.degree)
dec_lsst=coord.Angle(-70.74058,unit=u.degree)
ax.scatter(dec_lsst.degree, ra_lsst.degree,color='b')
"""
figb = plt.figure(figsize=(8,6))
axb = figb.add_subplot(111)
for key,num in thedict.items():
    for keyb,val in num.items():
        axb.plot(val['expMJD'],val['airmass'],col[typ]+'.')


plt.show()
