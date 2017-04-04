import numpy as np
from lsst.sims.maf.metrics import BaseMetric

from lsst.sims.photUtils import SignalToNoise
from lsst.sims.photUtils import PhotometricParameters
import sncosmo
from astropy.table import vstack,Table
from astropy import (cosmology, units as u, constants as const)
import astropy.units as u
import matplotlib.pyplot as plt
from lsst.sims.photUtils.EBV import EBVbase
from lsst.sims.photUtils import Bandpass
import math
import scipy
from lsst.sims.photUtils import Sed

from SN_Object import SN_Object
from Throughputs import Throughputs
import cPickle as pkl
from cosmology import GeneralCosmo

bands=['u','g','r','i','z','y']

 #FWHM_500 = seeing at 500 nm
 # FWHM_Sys_Zenith = sqrt(telSeeing**2 + opticalDesSeeing**2 + cameraSeeing**2)
 # Filter_Wavelength_Correction = (500 nm / Filter_Effective_Wavelength)**0.3
 # Airmass_Correction = airmass**0.6
 # FWHM_Sys = FWHM_Sys_Zenith * Airmass_Correction
 # FWHM_Atm = FWHM_500 * Filter_Wavelength_Correction * Airmass_Correction
 # FWHM_Eff = scaleToNeff * sqrt(FWHM_Sys**2 + atmNeffFactor * FWHM_Atm**2)
 # FWHM_Eff is the value in ObsHistory.finSeeing for the observations filter
 #
 # Units = unitless, Format = float, no default
 #

telSeeing = 0.250 # design goal
opticalDesSeeing = 0.08
cameraSeeing = 0.300
# Scaling factors for above seeing calculation
scaleToNeff = 1.16
atmNeffFactor = 1.04
FWHM_Sys_Zenith = np.sqrt(telSeeing**2 + opticalDesSeeing**2 + cameraSeeing**2)

filterWave = {'u': 367.0, 'g': 482.5, 'r': 622.2, 'i': 754.5, 'z': 869.1, 'y': 971.0}


# Read throughputs
transmission=Throughputs()

# Register LSST band pass (system) in sncosmo

for band in bands:
    band=sncosmo.Bandpass(transmission.lsst_system[band].wavelen, transmission.lsst_system[band].sb, name='LSST::'+band,wave_unit=u.nm)
    sncosmo.registry.register(band)


Delay_Obs={'u':0,'g':1,'r':2,'i':3,'z':4,'y':5}

T0=60000
Tmin=T0-30
Tmax=T0+80

#define a set of observations between T0-30 and T0+50

List=['Observations_DD_290.pkl']
thedict={}
outdir='../Obs_minion_1016'

for i,name in enumerate(List):
    pkl_file = open(outdir+'/'+name,'rb')
    thedict[i]=pkl.load(pkl_file)


print thedict[0]['dataSlice'].dtype

obs_filt={}

for band in bands:
    obs_filt[band]=thedict[0]['dataSlice'][np.where(thedict[0]['dataSlice']['filter']==band)]
    print band,len(obs_filt[band])

myobservations=np.zeros((50,1),dtype=thedict[0]['dataSlice'].dtype)

iobs=-1

for Time_obs in np.arange(Tmin, Tmax, 3.):
    for band in bands:
        nchoice= int(np.random.uniform(0,len(obs_filt[band])))
        print nchoice,obs_filt[band][nchoice]
        iobs+=1 
        if len( myobservations) <= iobs:
            myobservations=np.resize(myobservations,(len(myobservations)+50,1))
        for name in myobservations.dtype.names:
            myobservations[name][iobs]=obs_filt[band][name][nchoice]
        #correct for expMJD
        myobservations['expMJD'][iobs]=Time_obs+Delay_Obs[band]/24.
    #break

myobservations=np.resize(myobservations,(iobs+1,1))
print myobservations

ra=myobservations['fieldRA'][0][0]
dec=myobservations['fieldDec'][0][0]

print 'hello',ra,dec
redshift=0.01
redshiftb=0.01
c=0.0
x1=0.0
model='salt2-extended'
version='1.0'

sntype='Ia'

sn_dict={}

sn_dict[0]={}

sn_dict[0]['z']=0.01
sn_dict[0]['c']=0.0
sn_dict[0]['x1']=0.0

sn_dict[1]={}

sn_dict[1]['z']=0.01
sn_dict[1]['c']=-0.1
sn_dict[1]['x1']=0.0

for i,val in sn_dict.items():
    sn_dict['SN']=SN_Object(ra=np.rad2deg(ra),dec=np.rad2deg(dec),z=val['z'],t0=T0,c=val['c'],x1=val['x1'],model=model,version=version,sn_type=sntype)

SN=SN_Object(ra=np.rad2deg(ra),dec=np.rad2deg(dec),z=redshift,t0=T0,c=c,x1=x1,model=model,version=version,sn_type=sntype)
SNb=SN_Object(ra=np.rad2deg(ra),dec=np.rad2deg(dec),z=redshiftb,t0=T0,c=-0.1,x1=x1,model=model,version=version,sn_type=sntype)

table_for_fit={}
table_for_fit_b={}

table_for_fit['error_calc'] = Table(names=('time','flux','fluxerr','band','zp','zpsys'), dtype=('f8', 'f8','f8','S7','f4','S4'))
table_for_fit_b['error_calc'] = Table(names=('time','flux','fluxerr','band','zp','zpsys'), dtype=('f8', 'f8','f8','S7','f4','S4'))


table_mag={}
for i in range(len(sn_dict)):
    table_for_fit[i] = Table(names=('time','flux','fluxerr','band','zp','zpsys'), dtype=('f8', 'f8','f8','S7','f4','S4'))
    table_mag[i]=Table(names=('time','mag','band'), dtype=('f8', 'f8','S7'))

table_mag= Table(names=('time','mag','band'), dtype=('f8', 'f8','S7'))
table_magb= Table(names=('time','mag','band'), dtype=('f8', 'f8','S7'))

mycosmo=GeneralCosmo(0.2865,1.-0.2865,-1.,0.)
ratio_cosmo=np.power(mycosmo.Dl(redshiftb)/mycosmo.Dl(redshift),2.)*(1.+redshiftb)/(1.+redshift)
print 'lumidist',SN.lumidist,SNb.lumidist

for band in bands:
    observations_filtre=myobservations[np.where(myobservations['filter']==band)]
    
    for obs in observations_filtre:
        time=obs['expMJD']
        seeing=obs['rawSeeing']
        m5_opsim=obs['fiveSigmaDepth']

        #time=T0
        sed_SN=SN.get_SED(time)
        #timeb=(time-T0)*(1.+redshiftb)/(1.+redshift)+T0
        sed_SNb=SNb.get_SED(time)
        
        for i in range(len(dict_sn)):
            dict_SN[i]['SED']=dict_sn[i]['SN'].get_SED(time)


        transmission.Load_Atmosphere(obs['airmass'])
        
        flux_SN=sed_SN.calcFlux(bandpass=transmission.lsst_atmos_aerosol[band])
        flux_SNb=sed_SNb.calcFlux(bandpass=transmission.lsst_atmos_aerosol[band])

        for i in range(len(dict_sn)):
            dict_SN[i]['Flux']=dict_SN[i]['SED'].calcFlux(bandpass=transmission.lsst_atmos_aerosol[band])
        


        Filter_Wavelength_Correction = np.power(500.0 / filterWave[band], 0.3)
        Airmass_Correction = math.pow(obs['airmass'],0.6)
        FWHM_Sys = FWHM_Sys_Zenith * Airmass_Correction
        FWHM_Atm = seeing * Filter_Wavelength_Correction * Airmass_Correction
        finSeeing = scaleToNeff * math.sqrt(np.power(FWHM_Sys,2) + atmNeffFactor * np.power(FWHM_Atm,2))
        

        for i in range(len(sn_dict)):
            flux_SN=dict_SN[i]['Flux']
            if flux_SN > 0:
                sed_SN=dict_SN[i]['SED']

                mag_SN=-2.5 * np.log10(flux_SN / 3631.0)
            
            
                FWHMeff = SignalToNoise.FWHMgeom2FWHMeff(finSeeing)
                photParams = PhotometricParameters()
                snr_SN= SignalToNoise.calcSNR_sed(sed_SN,transmission.lsst_atmos_aerosol[band], transmission.darksky, transmission.lsst_system[band], 
                                                      photParams, FWHMeff=FWHMeff, verbose=False)
                m5_calc=SignalToNoise.calcM5(transmission.darksky,transmission.lsst_atmos_aerosol[band],transmission.lsst_system[band],photParams=photParams,FWHMeff=FWHMeff)
                snr_m5_through,gamma_through=SignalToNoise.calcSNR_m5(mag_SN,transmission.lsst_atmos_aerosol[band],m5_calc,photParams)
                snr_m5_opsim,gamma_opsim=SignalToNoise.calcSNR_m5(mag_SN,transmission.lsst_atmos_aerosol[band],m5_opsim,photParams)

                err_flux_SN=flux_SN/snr_SN
           
                table_for_fit[i].add_row((time,flux_SN,err_flux_SN,'LSST::'+band,25.,'ab'))
                table_mag[i].add_row((time,mag_SN,'LSST::'+band))







        #print 'flux_SN',band,flux_SN
        if flux_SN >0:
                            
            mag_SN=-2.5 * np.log10(flux_SN / 3631.0)
            #mag_SN=-2.5 * np.log10(flux_SN )
            
            
            FWHMeff = SignalToNoise.FWHMgeom2FWHMeff(finSeeing)
            photParams = PhotometricParameters()
            snr_SN= SignalToNoise.calcSNR_sed(sed_SN,transmission.lsst_atmos_aerosol[band], transmission.darksky, transmission.lsst_system[band], 
                                                      photParams, FWHMeff=FWHMeff, verbose=False)
            m5_calc=SignalToNoise.calcM5(transmission.darksky,transmission.lsst_atmos_aerosol[band],transmission.lsst_system[band],photParams=photParams,FWHMeff=FWHMeff)
            snr_m5_through,gamma_through=SignalToNoise.calcSNR_m5(mag_SN,transmission.lsst_atmos_aerosol[band],m5_calc,photParams)
            snr_m5_opsim,gamma_opsim=SignalToNoise.calcSNR_m5(mag_SN,transmission.lsst_atmos_aerosol[band],m5_opsim,photParams)

            err_flux_SN=flux_SN/snr_SN
           
            table_for_fit['error_calc'].add_row((time,flux_SN,err_flux_SN,'LSST::'+band,25.,'ab'))
            table_mag.add_row((time,mag_SN,'LSST::'+band))

        if flux_SNb >0:
                            
            mag_SNb=-2.5 * np.log10(flux_SNb / 3631.0)
            #mag_SNb=-2.5 * np.log10(flux_SNb)
            
            FWHMeff = SignalToNoise.FWHMgeom2FWHMeff(finSeeing)
            photParams = PhotometricParameters()
            snr_SNb= SignalToNoise.calcSNR_sed(sed_SNb,transmission.lsst_atmos_aerosol[band], transmission.darksky, transmission.lsst_system[band], 
                                                      photParams, FWHMeff=FWHMeff, verbose=False)
            m5_calc=SignalToNoise.calcM5(transmission.darksky,transmission.lsst_atmos_aerosol[band],transmission.lsst_system[band],photParams=photParams,FWHMeff=FWHMeff)
            snr_m5_through,gamma_through=SignalToNoise.calcSNR_m5(mag_SNb,transmission.lsst_atmos_aerosol[band],m5_calc,photParams)
            snr_m5_opsim,gamma_opsim=SignalToNoise.calcSNR_m5(mag_SNb,transmission.lsst_atmos_aerosol[band],m5_opsim,photParams)

            err_flux_SNb=flux_SNb/snr_SNb
            
            table_for_fit_b['error_calc'].add_row((time,flux_SNb,err_flux_SNb,'LSST::'+band,25.,'ab')) 
            table_magb.add_row((time,mag_SNb,'LSST::'+band))

print 'resultat',table_for_fit['error_calc'],table_for_fit_b['error_calc']

"""
#now do the fit...

z_sim=SN.z
res, fitted_model = sncosmo.fit_lc(table_for_fit['error_calc'], SN.SN_fit_model,['z', 't0', 'x0', 'x1', 'c'],bounds={'z':(z_sim-0.01, z_sim+0.01)})

mbfit=fitted_model._source.peakmag('bessellb','vega')

sncosmo.plot_lc(table_for_fit['error_calc'], model=fitted_model,color='k',pulls=False)

print res,fitted_model

z_simb=SNb.z
resb, fitted_modelb = sncosmo.fit_lc(table_for_fit_b['error_calc'], SN.SN_fit_model,['z', 't0', 'x0', 'x1', 'c'],bounds={'z':(z_simb-0.01, z_simb+0.01)})

mbfitb=fitted_modelb._source.peakmag('bessellb','vega')

sncosmo.plot_lc(table_for_fit_b['error_calc'], model=fitted_modelb,color='k',pulls=False)

print resb,fitted_modelb

print mbfit,mbfitb

"""

figa, axa = plt.subplots(ncols=2, nrows=3, figsize=(10,9))

for j,band in enumerate(['u','g','r','i','z','y']):
    if j<2:
        k=0
    if j>= 2 and j < 4:
        k=1
    if j>=4:
        k=2
           
      
    sel=table_for_fit['error_calc'][np.where(table_for_fit['error_calc']['band']=='LSST::'+band)]
    selb=table_for_fit_b['error_calc'][np.where(table_for_fit_b['error_calc']['band']=='LSST::'+band)] 
   
    print sel
    axa[k][j%2].errorbar(sel['time'],sel['flux'],xerr=0.0000001, yerr=sel['fluxerr'],fmt='-',color = 'k')
    axa[k][j%2].errorbar(selb['time'],selb['flux'],xerr=0.0000001, yerr=selb['fluxerr'],fmt='-',color = 'r')
   

figb, axb = plt.subplots(ncols=2, nrows=3, figsize=(10,9))

for j,band in enumerate(['u','g','r','i','z','y']):
    if j<2:
        k=0
    if j>= 2 and j < 4:
        k=1
    if j>=4:
        k=2
           
    """    
    sel=table_for_fit['error_calc'][np.where(table_for_fit['error_calc']['band']=='LSST::'+band)]
    selb=table_for_fit_b['error_calc'][np.where(table_for_fit_b['error_calc']['band']=='LSST::'+band)] 
    """
    sel=table_mag[np.where(table_mag['band']=='LSST::'+band)]
    selb=table_magb[np.where(table_magb['band']=='LSST::'+band)] 

    print sel
    #axb[k][j%2].errorbar(sel['time'],sel['flux'],xerr=0.0000001, yerr=sel['fluxerr'],fmt='-',color = 'k')
    #axb[k][j%2].errorbar(selb['time'],selb['flux'],xerr=0.0000001, yerr=selb['fluxerr'],fmt='-',color = 'r')
    axb[k][j%2].plot(sel['time'],sel['mag'],'ko')
    axb[k][j%2].plot(selb['time'],selb['mag'],'ro')
    axb[k][j%2].set_ylim(axb[k][j%2].get_ylim()[::-1])

plt.show()
