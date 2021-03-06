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

def Get_Seasons(filtc):

    """
    names_ref=('time','flux','fluxerr','band','zp','zpsys')
    dtype_ref=('f8', 'f8','f8','S7','f4','S4')

    dtype_ref=[type[1] for type in filtc.dtype]
    print 'hello',type(filtc),filtc.dtype.names,dtype_ref
    """

    dict_for_seasons={}
    if len(filtc) > 0:
        inum=0
        print filtc.dtype,filtc.shape,type(filtc)
        sorted_data=[]
        for data in filtc:
            if np.isscalar(data["expMJD"]):
                sorted_data.append(data["expMJD"])
            else:
                sorted_data.append(data["expMJD"][0])
        #sorted_data=sorted_data.sort(axis=1)
        ind =np.argsort(sorted_data)
        
        #print ind
        filtc=filtc[ind]
        """
        plt.plot(filtc['expMJD'],filtc['airmass'],'b.')
        plt.show()
        """
        #dict_for_seasons[inum]=Table(names=filtc.dtype.names,dtype=filtc.dtype)
        dict_for_seasons[inum]=np.zeros((60,1),dtype=filtc.dtype)
        #print 'timediff',24.*60.*60.*(filtc['time']-filtc['time'][0])
                                
        iloop=0
        iinside=0
        dict_for_seasons[inum][iinside]=filtc[iloop]
                                
        if len(filtc) > 1:
            while iloop < len(filtc)-1: 
                iinside+=1
                diff_time_days=filtc['expMJD'][iloop+1]-filtc['expMJD'][iloop]
                #print 'alors ???',diff_time_sec,inum
                if diff_time_days > 100.:
                    dict_for_seasons[inum]=np.resize(dict_for_seasons[inum],iinside)
                    inum+=1
                    #dict_for_seasons[inum]=Table(names=filtc.dtype.names, dtype=filtc.dtype)
                    dict_for_seasons[inum]=np.zeros((60,1),dtype=filtc.dtype)
                    iinside=0
                #dict_for_seasons[inum].add_row(filtc[iloop+1])
                if len(dict_for_seasons[inum]) <= iinside:
                    dict_for_seasons[inum]=np.resize(dict_for_seasons[inum],(len(dict_for_seasons[inum])+50,1))

                dict_for_seasons[inum][iinside]=filtc[iloop+1]
   
                iloop+=1
        #print 'thedict',dict_for_seasons
                
                
    return dict_for_seasons

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

def Plot_bands(obj,time_name='expMJD',flux_name='flux',errflux_name='err_flux',filter_name='filter',addit='',opsim=None,T0=-1):

    figa, axa = plt.subplots(ncols=2, nrows=3, figsize=(10,9))
    
    for j,band in enumerate(['u','g','r','i','z','y']):
        if j<2:
            k=0
        if j>= 2 and j < 4:
            k=1
        if j>=4:
            k=2
        bandsel=addit+band
        selobs=obj[np.where(obj[filter_name]==bandsel)]
        if bandsel=='g':
            print selobs[flux_name],selobs[errflux_name]
        #selobs=selobs[np.where(selobs[flux_name]/selobs[errflux_name]>5)]
        #axa[k][j%2].errorbar(selobs[time_name]-selobs[time_name].min(),selobs[flux_name],yerr=selobs[errflux_name],fmt='.',ecolor='r')
        axa[k][j%2].errorbar(selobs[time_name],selobs[flux_name],yerr=selobs[errflux_name],fmt='.',ecolor='r',color='r')
        axca = axa[k][j%2].twinx()
        opsim_sel=opsim[np.where(opsim['filter']==band)]
        axca.plot(opsim_sel['expMJD'],opsim_sel['fieldRA'],'b.') 
        axa[k][j%2].set_xlim(T0-30,T0+50)
       

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

file='SuperNova_'+opts.sntype+'_'+opts.fieldname+'_'+str(opts.fieldid)+'_'+str(opts.zmin)+'_'+str(opts.zmax)+'_'+str(opts.nevts)+'_season_'+str(opts.season)+'_3.pkl'

#file='SuperNova_'+sntype+'_WFD_309_'+str(zmin)+'_'+str(zmax)+'_salt2-extended_10_11.pkl'
opsim_release='minion_1016'
thedir='../Make_Cadence/Sim_'+opsim_release
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


val='error_coadd_through'

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
Compare_Errors=False
Get_m5=False
Visu_Observations=True

if Visu_Observations:

    #load Opsim obs so as to compare with the LC
    opsim_dir='/sps/lsst/data/dev/pgris/Make_Cadence/Obs_minion_1016'
    
    prefix='Observations'
    
    List=[prefix+'_'+opts.fieldname+'_'+str(opts.fieldid)+'.pkl']
    
    thedict={}

    for i,name in enumerate(List):
        pkl_file = open(opsim_dir+'/'+name,'rb')
        thedict[i]=pkl.load(pkl_file)

    seasons=Get_Seasons(thedict[0]['dataSlice'])
    thetab=Table(names=('T0','x1','c','z','nbefore_g','nafter_g'), dtype=('f8', 'f8','f8','f8','i4','i4'))

    thefiltre='z'
    dict_obs={}
    for oob in all_obs:
        for i,obj in enumerate(oob):
            #print i,obj['status']
            if obj['observations'] is not None:
                #print i,obj['observations']['flux']
               
                #Plot_bands(obj['observations'])

                dict_fit=obj['fit']
                if dict_fit is not None:
                    dict_tag=dict_fit[val]
                    #print dict_tag['fit_status']

                    filt=dict_tag['table_for_fit']
                    filtb=filt[np.where(np.logical_and(filt['flux']/filt['fluxerr']>5.,filt['flux']>0.))]
                    #Plot_bands(filtb,time_name='time',errflux_name='fluxerr',filter_name='band',addit='LSST::')
                    sel_for_fit=filtb[np.where(filtb['band']=='LSST::'+thefiltre)]
                                     
                    sel_before=sel_for_fit[np.where(sel_for_fit['time']-obj['t0']<=0.)]
                    sel_after=sel_for_fit[np.where(sel_for_fit['time']-obj['t0']>0.)]
                    n_before=len(sel_before)
                    n_after=len(sel_after)

                    thetab.add_row((obj['t0'],obj['x1'],obj['c'],obj['z'],n_before,n_after))
                    dataSlice=seasons[opts.season]
                    #print 'redshift',obj['z']
                    timelow=obj['t0']-30.
                    timehigh=obj['t0']+50.

                    if obj['t0'] > 0.9999999*61974.7537142 and obj['t0'] < 1.0000001*61974.7537142:
                        #print obj
                        
                        theobservations=obj['observations'][np.where(np.logical_and(obj['observations']['expMJD']>timelow,obj['observations']['expMJD']<timehigh))]
                        thefilt=theobservations[np.where(theobservations['filter']==thefiltre)]
                        print 'hello',len(thefilt),obj['t0'],obj['x1'],obj['c']
                        print thefilt['expMJD'],thefilt['flux'],thefilt['flux']/thefilt['err_flux']
                        
                    if n_before == 1 and obj['t0']-np.min(dataSlice['expMJD'])>=20:
                        
                        observations=dataSlice[np.where(np.logical_and(dataSlice['expMJD']>timelow,dataSlice['expMJD']<timehigh))]
                        #print 'T0',obj['t0']
                        Plot_This=False
                        if Plot_This:
                            Plot_bands(obj['observations'],opsim=seasons[opts.season],T0=obj['t0'])
                            Plot_bands(filtb,time_name='time',errflux_name='fluxerr',filter_name='band',addit='LSST::',opsim=seasons[opts.season],T0=obj['t0'])
                        #Plot_bands(filtb,time_name='time',errflux_name='fluxerr',filter_name='band',addit='LSST::',opsim=observations,T0=obj['t0'])
                        #print obj['observations']
                            plt.show()
              

    observations=seasons[opts.season]
    sel_obs=observations[np.where(observations['filter']=='g')]
    tmin=np.min(sel_obs['expMJD'])
    fontsize=15.
    figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    sel=thetab[np.where(thetab['nbefore_g']>=1)]
    axa.plot(sel['T0']-tmin,sel['nbefore_g'],'bo')
    axa.set_xlabel(r'$T_{0}-T_{obs}$$^{min}$/$T-T_{obs}$$^{min}$',{'fontsize': fontsize})
    axa.set_ylabel(r'N$_{measurements}$ before T$_0$',{'fontsize': fontsize})
    axa.set_ylim(0.5, 5.5)
    axca=axa.twinx()
    #what_obs='fiveSigmaDepth'
    what_obs='fieldRA'
    #what_obs='c'
    #axca.plot(sel['T0']-tmin,sel[what_obs],'ro')
    axca.plot(sel_obs['expMJD']-tmin,sel_obs[what_obs],'ro')
    print 'hh',sel_obs['expMJD']-tmin
    axca.set_ylabel(r''+what_obs,{'fontsize': fontsize})
    #axca.set_ylim(5.4, 6.5)
    figa.suptitle('Field '+opts.fieldname+' - '+str(opts.fieldid)+' Filter '+'g'+' - '+str(opts.zmin)+'< z <'+str(opts.zmax)+' - Season '+str(opts.season))

    figb, axb = plt.subplots(ncols=3, nrows=2, figsize=(10,9))
    selb=thetab[np.where(thetab['nbefore_g']==1)]
    selb=thetab

    #seld=selb[np.where(np.logical_and(selb['T0']-tmin >17,selb['T0']-tmin<26))]
    seld=selb[np.where(np.logical_and(selb['T0']-tmin >44.5,selb['T0']-tmin<45.))]
    print seld['T0'],seld['x1'],seld['c']


    axb[0][0].plot(seld['x1'],seld['nbefore_g'],'bo')
    axb[0][1].plot(seld['c'],seld['nbefore_g'],'bo')
    axb[0][2].plot(seld['z'],seld['nbefore_g'],'bo')

    selc=seld[np.where(seld['nbefore_g']==1)]
    sele=seld[np.where(seld['nbefore_g']==2)]
    
    #print 'ee',selc['T0'],selc['x1'],selc['c'],selc['z']

    #print selc.dtype

    


    axb[1][0].hist(selc['x1'],bins=15,histtype='step',fill=False,color='r')
    axb[1][0].hist(sele['x1'],bins=15,histtype='step',fill=False,color='k')
   
    axb[1][1].hist(selc['c'],bins=15,histtype='step',fill=False,color='r')
    axb[1][1].hist(sele['c'],bins=15,histtype='step',fill=False,color='k')
    
    axb[1][2].plot(selc['x1'],selc['c'],'ro')
    axb[1][2].plot(sele['x1'],sele['c'],'ko')

    """
    figb, axb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    sel=thetab[np.where(sel['nbefore_g']==1)]
    axb.plot(sel['T0']-tmin,sel['x1'],'bo')
    """
    plt.show()

if Visu_LC:
    for oob in all_obs:
        for i,obj in enumerate(oob):
            print i,obj['status']
            if obj['observations'] is not None:
                print i,obj['observations']['flux']
            
            dict_fit=obj['fit']
            if dict_fit is not None:
                dict_tag=dict_fit[val]
                print dict_tag['fit_status']
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

               
                
                    if dict_tag['sncosmo_res']['ndof']>0:
                        
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
                            nmeas_before[band]=len(filtcc)
                            nmeas_after[band]=len(filtcd)

                        pass_event=False
                    
                        test_one = nmeas_before['g'] >= n_meas and nmeas_after['g'] >= n_meas
                        test_two = nmeas_before['r'] >= n_meas and nmeas_after['r'] >= n_meas
                        test_three = nmeas_before['i'] >= n_meas and nmeas_after['i'] >= n_meas
                        if (test_one and test_two) or (test_two and test_three):
                            pass_event=True

                        if pass_event:

                            dust = sncosmo.OD94Dust()                                                                
                            SN=sncosmo.Model(source='salt2-extended',effects=[dust, dust],effect_names=['host', 'mw'],effect_frames=['rest', 'obs']) 
                            SN.set(z=obj['z'])
                            SN.set(t0=obj['t0'])
                            SN.set(c=obj['c']) 
                            SN.set(x1=obj['x1']) 
                            SN.set_source_peakabsmag(-19.3, 'bessellB', 'vega',cosmo=cosmology.WMAP9)
                            lsstmwebv = EBVbase()
                            ebvofMW = lsstmwebv.calculateEbv(equatorialCoordinates=np.array([[obj['ra']], [obj['dec']]]))[0]
                            print 'hello',obj['ra'],obj['dec'],ebvofMW
                            SN.set(mwebv=ebvofMW)
                            print 'simulated',obj['z'],obj['t0'],obj['c'],obj['x1']
                            resb, fitted_modelb = sncosmo.fit_lc(dict_tag['table_for_fit'],SN,['z', 't0', 'x0', 'x1', 'c'],bounds={'z':(obj['z']-0.01,obj['z']+0.01)})
                            sncosmo.plot_lc(dict_tag['table_for_fit'], model=fitted_modelb,color='k',pulls=True,errors=resb.errors) 

                            plt.show()

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

                    colors=['g','b']
                    if dict_fit is not None:
                    
                        #for pp,val in enumerate(['error_coadd_calc','error_coadd_opsim']):
                        for pp,val in enumerate(['error_coadd_calc','error_coadd_through']):
                            dict_tag=dict_fit[val]
                            selcoadd=dict_tag['table_for_fit'][np.where(dict_tag['table_for_fit']['band']=='LSST::'+band)]
                            axa[k][j%2].plot(selcoadd['time']-selcoadd['time'].min(),selcoadd['flux']/selcoadd['fluxerr'],colors[pp]+'.')
                            print 'fit status',val,selcoadd
                plt.show()
if Get_m5:
    dtype=[]
    for oob in all_obs:
        for i,obj in enumerate(oob):
            print i,obj['status']
            dict_fit=obj['fit']
            if obj['status'] != 'Killed':
                dtype=obj['observations'].dtype
                break
                
    

    resu=np.zeros((0,1),dtype=dtype)
    print 'hello',resu.dtype
    iop=-1
    for oob in all_obs:
        iop+=1
        for i,obj in enumerate(oob):
            print i,obj['status']
            dict_fit=obj['fit']
            if obj['status'] != 'Killed'and obj['observations']!= None:
                if iop == 0:
                    print 'bbrbrbrbr',obj['observations']
                    resu=obj['observations']
                else:
                    resu=np.concatenate((resu,obj['observations']),axis=1)

    print 'alors',resu

    figa, axa = plt.subplots(ncols=2, nrows=3, figsize=(10,9))
    print resu.dtype
    for j,band in enumerate(['u','g','r','i','z','y']):
        if j<2:
            k=0
        if j>= 2 and j < 4:
            k=1
        if j>=4:
            k=2
                    
        selobs=resu[np.where(resu['filter']==band)]
        print 'hello',band,len(selobs)
        axa[k][j%2].plot(selobs['snr_m5_opsim'],selobs['snr_m5_through'],'k.')
    plt.show()

