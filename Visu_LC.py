import cPickle as pkl
import sncosmo
import matplotlib.pyplot as plt
from astropy import (cosmology, units as u, constants as const)
from Throughputs import Throughputs
import numpy as np
from astropy.table import vstack,Table
from lsst.sims.photUtils.EBV import EBVbase
from optparse import OptionParser
from astropy.table import vstack,Table

def Get_coadd(filtc):
    #names_ref=('time','flux','fluxerr','band','zp','zpsys')
    dtype_ref=('f8', 'f8','f8','S7','f4','S4')
    names_ref=filtc.colnames
    #dtype_ref=filtc.dtype

    #print 'before filtering',filtc
    out_table=Table(names=names_ref,dtype=dtype_ref)
    dict_for_coadd={}
    if len(filtc) > 0:
        inum=0
        dict_for_coadd[inum]=Table(names=names_ref,dtype=dtype_ref)
        #print 'timediff',24.*60.*60.*(filtc['time']-filtc['time'][0])
                                
        iloop=0
        #print 'blablabla',dict_for_coadd[inum]
        dict_for_coadd[inum].add_row(filtc[iloop])
                                
        if len(filtc) > 1:
            while iloop < len(filtc)-1:   
                diff_time_sec=24.*60.*60.*(filtc['time'][iloop+1]-filtc['time'][iloop])
                #print 'alors ???',diff_time_sec,inum
                if diff_time_sec > 40.:
                    inum+=1
                    dict_for_coadd[inum]=Table(names=names_ref,dtype=dtype_ref)
                
                dict_for_coadd[inum].add_row(filtc[iloop+1])
                    
                iloop+=1
        #print 'thedict',dict_for_coadd

    for key,vals in dict_for_coadd.items():
        mean_pond=np.sum(vals['flux']*(np.power(vals['fluxerr'],-2))/np.sum(np.power(vals['fluxerr'],-2)))
        sigsum=1./np.sqrt(np.sum(np.power(vals['fluxerr'],-2)))
        out_table.add_row((np.mean(vals['time']),mean_pond,sigsum,vals['band'][0],vals['zp'][0],vals['zpsys'][0]))

    #print 'after filtering',out_table

    return out_table

parser = OptionParser()

parser.add_option("-N", "--nevts", type="int", default=1000, help="filter [%default]")
parser.add_option("-f", "--fieldname", type="string", default="DD", help="filter [%default]")
parser.add_option("-n", "--fieldid", type="int", default=290, help="filter [%default]")
parser.add_option("-s", "--season", type="int", default=7, help="filter [%default]")
parser.add_option("-t", "--sntype", type="string", default="Ia", help="filter [%default]")
parser.add_option("-r", "--rolling", type="int", default=0, help="filter [%default]")
parser.add_option("-z", "--zmin", type="float", default=0.7, help="filter [%default]")
parser.add_option("-Z", "--zmax", type="float", default=0.8, help="filter [%default]")

opts, args = parser.parse_args()

file='SuperNova_'+opts.sntype+'_'+opts.fieldname+'_'+str(opts.fieldid)+'_'+str(opts.zmin)+'_'+str(opts.zmax)+'_'+str(opts.nevts)+'_season_'+str(opts.season)+'_0.pkl'

#file='SuperNova_'+sntype+'_WFD_309_'+str(zmin)+'_'+str(zmax)+'_salt2-extended_10_11.pkl'
opsim_release='minion_1016'
thedir='Sim_'+opsim_release
filenames=[]

filenames.append(thedir+'/'+file)

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


val='error_coadd_calc'

dust = sncosmo.OD94Dust()
fitted_model=sncosmo.Model(source='salt2-extended', effects=[dust, dust],
                    effect_names=['host', 'mw'],
                    effect_frames=['rest', 'obs'])


fitted_model_coadd=sncosmo.Model(source='salt2-extended', effects=[dust, dust],
                    effect_names=['host', 'mw'],
                    effect_frames=['rest', 'obs'])
 # Read throughputs

transmission=Throughputs()

# Register LSST band pass (system) in sncosmo
bands= ['u','g','r','i','z','y']

for filtre in bands:
    band=sncosmo.Bandpass(transmission.lsst_system[filtre].wavelen, transmission.lsst_system[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm)
    sncosmo.registry.register(band)

Visu_LC=False
Compare_Errors=True

if Visu_LC:
    for oob in all_obs:
        for i,obj in enumerate(oob):
            print i,obj['status']
            print i,obj['observations']['flux']
            
            dict_fit=obj['fit']
            if dict_fit is not None:
                dict_tag=dict_fit[val]
                print dict_tag
                if dict_tag['fit_status'] == 'ok':
                    fitted_model.set(z=dict_tag['sncosmo_fitted']['z'])
                    fitted_model.set(t0=dict_tag['sncosmo_fitted']['t0'])
                    fitted_model.set(x0=dict_tag['sncosmo_fitted']['x0'])
                    fitted_model.set(x1=dict_tag['sncosmo_fitted']['x1'])
                    fitted_model.set(c=dict_tag['sncosmo_fitted']['c'])
                    fitted_model.set(hostebv=dict_tag['sncosmo_fitted']['hostebv'])
                    fitted_model.set(hostr_v=dict_tag['sncosmo_fitted']['hostr_v'])
                    fitted_model.set(mwebv=dict_tag['sncosmo_fitted']['mwebv'])
                    fitted_model.set(mwr_v=dict_tag['sncosmo_fitted']['mwr_v'])

               
                
                    if dict_tag['sncosmo_res']['ndof']>0 and dict_tag['sncosmo_res']['chisq']/dict_tag['sncosmo_res']['ndof']>5.:
                        
                        filt=dict_tag['table_for_fit']
                        filtb=filt[np.where(np.logical_and(filt['flux']/filt['fluxerr']>5.,filt['flux']>0.))]

                        T0=dict_tag['sncosmo_fitted']['t0']
                        print 'Nobs',len(filt),len(filtb)
                        n_meas=2
                        nmeas_before={}
                        nmeas_after={}
                    
                
                        for band in ['g','r','i']:
                            filtc=filtb
                            filtc=filtc[np.where(filtc['band']=='LSST::'+band)]
                            filtcf=filtc
                            filtcc=filtcf[np.where(filtcf['time']-T0<=0.)]
                            filtcd=filtcf[np.where(filtcf['time']-T0>0.)]
                        #print band,len(filtc),len(filtcc),len(filtcd)
                            nmeas_before[band]=len(filtcc)
                            nmeas_after[band]=len(filtcd)
                        #nmeas_before[band]=len(filtc)
                        #nmeas_after[band]=len(filtc)
                            """
                            if band == 'g':
                           print Get_coadd(filtc) 
                        """   

                        pass_event=True
                    

                        for band in ['g','r','i']:
                            if nmeas_before[band] < n_meas or nmeas_after[band] < n_meas:
                                pass_event = False
                                #break

                        if pass_event:
                            zdiff=obj['z']-dict_tag['sncosmo_fitted']['z']
                            print 'filtrons ',pass_event,obj['z'],dict_tag['sncosmo_fitted']['z'],zdiff
                            for band in ['g','r','i']:
                                print band,nmeas_before[band],nmeas_after[band]

                            if zdiff < 100000:
                                print 'plotting',zdiff
                                select=dict_tag['table_for_fit'][np.where(dict_tag['table_for_fit']['flux']/dict_tag['table_for_fit']['fluxerr']>5.)]
                                sncosmo.plot_lc(select, model=fitted_model,color='k',pulls=True,errors=dict_tag['sncosmo_res'].errors)
                                
                                """
                                dust = sncosmo.OD94Dust()
                                SN=sncosmo.Model(source='salt2-extended',effects=[dust, dust],
                                                 effect_names=['host', 'mw'],
                                                 effect_frames=['rest', 'obs'])

                                SN.set(z=obj['z'])
                                SN.set(t0=obj['t0'])
                                SN.set(c=obj['c'])
                                SN.set(x1=obj['x1'])

                                #SN.set_source_peakabsmag(-19.3, 'bessellB', 'vega',cosmo=cosmology.WMAP9)
                                SN.set_source_peakabsmag(25, 'bessellB', 'ab',cosmo=cosmology.WMAP9)
                                
                                lsstmwebv = EBVbase()
                                ebvofMW = lsstmwebv.calculateEbv(equatorialCoordinates=np.array([[obj['ra']], [obj['dec']]]))[0]

                                 
                                print 'hello',obj['ra'],obj['dec'],ebvofMW
                                SN.set(mwebv=ebvofMW)
                               
                                print 'simulated',obj['z'],obj['t0'],obj['c'],obj['x1']

                                resb, fitted_modelb = sncosmo.fit_lc(coadd,SN,['z', 't0', 'x0', 'x1', 'c'],bounds={'z':(0.9*zmin, 1.1*zmax)})
                
                                sncosmo.plot_lc(coadd, model=fitted_modelb,color='k',pulls=True,errors=resb.errors)

                                print resb
                                """
                        figbb, axbb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
                        myobs=obj['observations'][np.where(obj['observations']['filter']=='g')]
                        axbb.plot(myobs['expMJD'],myobs['airmass'],'b.')

                        figbc, axbc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
                        axbc.plot(myobs['airmass'],myobs['err_flux'],'b.')
                    #axbc.plot(myobs['airmass'],myobs['snr_m5_through'],'b.')

                        figbd, axbd = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
                        axbd.plot(myobs['filtSkyBrightness'],myobs['err_flux'],'b.')
                        
                        plt.show()
            #break

if Compare_Errors:
     for oob in all_obs:
        for i,obj in enumerate(oob):
            print i,obj['status']
            dict_fit=obj['fit']
            if obj['status'] != 'Killed':
                print i,obj['observations']
                #break
                figa, axa = plt.subplots(ncols=2, nrows=3, figsize=(10,9))

                #for j,band in enumerate(bands):
                for j,band in enumerate(['r']):
                    if j<2:
                        k=0
                    if j>= 2 and j < 4:
                        k=1
                    if j>=4:
                        k=2
                        
                    selobs=obj['observations'][np.where(obj['observations']['filter']==band)]
                    print 'hello',len(selobs['expMJD']),len(selobs['flux']),len(selobs['err_flux']),len(selobs['err_flux_opsim'])
                    if len(selobs) > 0:
                        """
                        axa[k][j%2].errorbar(selobs['expMJD']-selobs['expMJD'].min(),selobs['flux'],yerr=selobs['err_flux'],fmt='.',ecolor='r')
                        axa[k][j%2].errorbar(selobs['expMJD']-selobs['expMJD'].min(),selobs['flux'],yerr=selobs['err_flux_opsim'],fmt='.',ecolor='k')
                        """
                        axa[k][j%2].plot(selobs['expMJD']-selobs['expMJD'].min(),selobs['flux']/selobs['err_flux'],'ko')
                        print 'before',selobs
                        table_for_fit= Table(names=('time','flux','fluxerr','band','zp','zpsys'), dtype=('f8', 'f8','f8','S7','f4','S4'))
                        for val in selobs:
                            table_for_fit.add_row((val['expMJD'],val['flux'],val['err_flux'],'LSST::'+band,25,'ab'))
                        coadd=Get_coadd(table_for_fit)
                        print 'coadd',coadd
                        axa[k][j%2].plot(coadd['time']-coadd['time'].min(),coadd['flux']/coadd['fluxerr'],'ro')
                        #axa[k][j%2].plot(selobs['expMJD']-selobs['expMJD'].min(),selobs['flux']/selobs['err_flux_opsim'],'bo') 

                    if dict_fit is not None:
                    
                        #for pp,val in enumerate(['error_coadd_calc','error_coadd_opsim']):
                        for pp,val in enumerate(['error_coadd_calc']):
                            dict_tag=dict_fit[val]
                            selcoadd=dict_tag['table_for_fit'][np.where(dict_tag['table_for_fit']['band']=='LSST::'+band)]
                            axa[k][j%2].plot(selcoadd['time']-selcoadd['time'].min(),selcoadd['flux']/selcoadd['fluxerr'],'g.')
                            print 'fit status',selcoadd
                plt.show()
