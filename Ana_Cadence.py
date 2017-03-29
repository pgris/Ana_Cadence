import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import ICRS
import math
import healpy as hp
from astropy.table import Table
from optparse import OptionParser
from Parameters import parameters
#from mpl_toolkits.basemap import Basemap


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

def Get_Nearest(orig_table, val):
   
    table=Table(orig_table)
    c = SkyCoord(ra=val['fieldRA']*u.radian, dec=val['fieldDec']*u.radian)  
    catalog = SkyCoord(ra=table['fieldRA']*u.radian, dec=table['fieldDec']*u.radian)  
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)

    #print 'astropy matching',idx,d2d,d3d,len(table),type(idx),table
    theres=[table['fieldID'][int(idx)],table['fieldRA'][int(idx)],table['fieldDec'][int(idx)]]
    table.remove_row(int(idx))
    return table,theres

inputdir='/sps/lsst/data/dev/pgris/Obs_minion_1016'

prefix='Observations'

#fieldtypes=['WFD','GalacticPlane','SouthCelestialPole-18','NorthEclipticSpur-18c','DDF']
fieldtypes=['WFD']

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

Draw_Projection=False
Draw_Ra_Dec=False
Draw_Visits=False
Make_Square=True


if Make_Square:

    pos_tab= Table(names=('fieldID','fieldRA','fieldDec'), dtype=('i4', 'f8','f8'))

    for key,num in thedict.items():
        for keyb,val in num.items():
            #print keyb,val['fieldRA'][0],val['fieldDec'][0]
            pos_tab.add_row((keyb,val['fieldRA'][0],val['fieldDec'][0]))
            
    #print pos_tab
    #print np.min(pos_tab['fieldRA'])

   
    ra_step=2. # degrees
    n_dec_zone=5 # number of region in dec = number of fields to be merged


    ra_min=np.min(pos_tab['fieldRA'])
    
    ra_dec_strip={}
    istrip=-1
    
    while ra_min < 360.-ra_step:
        istrip+=1
        ra_dec_strip[istrip]={}
        ra_max=ra_min+ra_step
        sel=pos_tab[np.where(np.logical_and(np.rad2deg(pos_tab['fieldRA'])>=ra_min,np.rad2deg(pos_tab['fieldRA'])<ra_max))]
        sel.sort('fieldDec')
        num_per_part=len(sel)/n_dec_zone
        ntag=0
        for count in range(n_dec_zone):
            if count == n_dec_zone-1:
                ra_dec_strip[istrip][count]=sel[ntag:]
            else:
               ra_dec_strip[istrip][count]=sel[ntag:ntag+num_per_part] 
            #print count, len(sel),num_per_part,len(ra_dec_strip[istrip][count])
            ntag+=num_per_part

        ra_min+=ra_step
        break
     
    
    
    plt.plot(ra_dec_strip[0][0]['fieldRA'],ra_dec_strip[0][0]['fieldDec'],'bo')
    plt.plot(ra_dec_strip[0][1]['fieldRA'],ra_dec_strip[0][1]['fieldDec'],'ro')
    plt.plot(ra_dec_strip[0][2]['fieldRA'],ra_dec_strip[0][2]['fieldDec'],'go')
    #plt.plot(ra_strip[istrip]['fieldRA'],ra_strip[istrip]['fieldDec'],'bo')


    #now try to match dec pops
    ra_dec_final_combi={}
    all_combo=[]
    for iv in range(len(ra_dec_strip)):
        #print 'o yes',iv
        ra_dec_final_combi[iv]={}
        strip_copy={}
        icombi=-1
        for i in range(1,n_dec_zone):
            strip_copy[i]=ra_dec_strip[iv][i].copy()
        for val in ra_dec_strip[iv][0]:
            icombi+=1
            restable=Table(names=('fieldID','fieldRA','fieldDec'), dtype=('i4', 'f8','f8'))
            restable.add_row((val))
            for  iu in range(1,n_dec_zone):
                strip_copy[iu],resu=Get_Nearest(strip_copy[iu],val)
                restable.add_row((resu))
            ra_dec_final_combi[iv][icombi]=restable

    for key,vval in ra_dec_final_combi.items():
        
        for key,val in vval.items():
            local=[]
            for ik in range(len(val)):
                local.append(val['fieldID'][ik])
            all_combo.append(local)

    print all_combo


    #print 'alors final',len(ra_dec_final_combi),len(ra_dec_final_combi[0])
    figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    colors=['bo','ro','go','k*','b*','r*','g*','ks','bs','rs','gs','ks']
    for iv in range(len(ra_dec_final_combi)):
        for ii in range(len(ra_dec_final_combi[iv])):
            axc.plot(np.rad2deg(ra_dec_final_combi[iv][ii]['fieldRA']),np.rad2deg(ra_dec_final_combi[iv][ii]['fieldDec']),colors[ii])
    
    axc.set_xlim(-1.,360.)
    plt.show()
    


    """
    dec_min=np.min(pos_tab['fieldDec'])
    for i in range(n_dec_zone):
        
        dec_max=dec_min+dec_step
        
        
        dec_min+=step
    """
    


#ax = fig.add_subplot(111)

col={}
col['DDF']='r'
col['WFD']='b'
col['GalacticPlane']='k'
col['SouthCelestialPole-18']='g'
col['NorthEclipticSpur-18c']='y'

if Draw_Projection:
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111,projection="aitoff")
    
    for typ in fieldtypes:
    # As next step, those coordinates are transformed into an astropy.coordinates
    # astropy.coordinates.SkyCoord object.
        c = SkyCoord(ra=np.rad2deg(fieldRA[typ]), dec=np.rad2deg(fieldDec[typ]),unit='degree', frame='icrs')

        ra_rad = c.ra.wrap_at(-1800 * u.deg).radian
    #ra_rad=c.ra.radian
        dec_rad = c.dec.radian


        thelabel=typ+' - '+str(round(100.*pourcent[typ],1))+'%'
        tot_label.append(ax.scatter(ra_rad, dec_rad,color=col[typ],label=thelabel))
    #ax.set_xlim(0., 360.)
    labs = [l.get_label() for l in tot_label]
    ax.legend(tot_label, labs, ncol=2,loc='lower right',prop={'size':5},frameon=False)
    ax.legend(bbox_to_anchor=(0.5, -0.1), loc=2, borderaxespad=0.,fontsize=10.)
    ax.grid(True)

if Draw_Ra_Dec:
    figc = plt.figure(figsize=(8,6))
    axc = figc.add_subplot(111)
    
    for typ in fieldtypes:
        thelabel=typ
        axc.scatter(np.rad2deg(fieldRA[typ]), np.rad2deg(fieldDec[typ]),color=col[typ],label=thelabel)
    #ax.set_xlim(360., 0.)

    axc.legend(loc=2, borderaxespad=0.,fontsize=10.)



"""
ra_lsst=coord.Angle(-30.24287,unit=u.degree)
dec_lsst=coord.Angle(-70.74058,unit=u.degree)
ax.scatter(dec_lsst.degree, ra_lsst.degree,color='b')

figb = plt.figure(figsize=(8,6))
axb = figb.add_subplot(111)

for key,num in thedict.items():
    for keyb,val in num.items():
        axb.plot(val['expMJD'],val['airmass'],col[typ]+'.')
"""

if Draw_Visits:
    nvisits={}
    for typ in fieldtypes:
        nvisits[typ]={}
        for band in bands:
            nvisits[typ][band]=[]
        
    for key,num in thedict.items():
        for keyb,val in num.items():
            for band in bands:
                sel=val[np.where(val['filter']==band)]
                nvisits[key][band].append(len(sel))

    mean_visit={}
    rms_visit={}
    for typ in fieldtypes:
        mean_visit[typ]=[]
        rms_visit[typ]=[]
        for band in bands:
            print typ,band,np.mean(nvisits[typ][band]),np.std(nvisits[typ][band])
        #print nvisits[typ][band]
            mean_visit[typ].append(np.mean(nvisits[typ][band]))
            rms=np.std(nvisits[typ][band])
            if rms < 0.001:
                rms=0.001
            rms_visit[typ].append(rms)
  
    myfmt={}
    myfmt['DDF']='--o'
    myfmt['WFD']='--s'
    myfmt['GalacticPlane']='--*'
    myfmt['SouthCelestialPole-18']='--+'
    myfmt['NorthEclipticSpur-18c']='--p'
    
    filters=[0,1,2,3,4,5]
    fige, axe = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    for typ in fieldtypes:
        axe.errorbar(filters,mean_visit[typ],yerr=rms_visit[typ],fmt=myfmt[typ],color =col[typ],label=typ)

    fontsize=12
    axe.set_ylabel(r'<N$_{visit}$> (per field - 10 years)',{'fontsize': fontsize})
    axe.set_xlabel(r'Filter',{'fontsize': fontsize})
    axe.set_xlim(-0.5, 5.5)
    axe.set_ylim(10.,5000.)
    axe.set_yscale('log')
    axe.legend(loc='center left', fontsize=12.)

plt.show()
