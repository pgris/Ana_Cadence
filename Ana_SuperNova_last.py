import cPickle as pkl
from astropy.table import Table
#import sncosmo
import numpy as np
import matplotlib.pyplot as plt
#from Throughputs import Throughputs
import astropy.units as u
from astropy import (cosmology, units as u, constants as const)
from optparse import OptionParser

"""
filename='SuperNova_4.189756_-1.082474_1000.pkl'
pkl_file = open(filename,'rb')
"""

def Get_coadd(filtc):
    names_ref=('time','flux','fluxerr','band','zp','zpsys')
    dtype_ref=('f8', 'f8','f8','S7','f4','S4')

    out_table=Table(names=names_ref, dtype=dtype_ref)
    dict_for_coadd={}
    if len(filtc) > 0:
        inum=0
        dict_for_coadd[inum]=Table(names=names_ref, dtype=dtype_ref)
        #print 'timediff',24.*60.*60.*(filtc['time']-filtc['time'][0])
                                
        iloop=0
        dict_for_coadd[inum].add_row(filtc[iloop])
                                
        if len(filtc) > 1:
            while iloop < len(filtc)-1:   
                diff_time_sec=24.*60.*60.*(filtc['time'][iloop+1]-filtc['time'][iloop])
                #print 'alors ???',diff_time_sec,inum
                if diff_time_sec > 40.:
                    inum+=1
                    dict_for_coadd[inum]=Table(names=names_ref, dtype=dtype_ref)
                
                dict_for_coadd[inum].add_row(filtc[iloop+1])
                    
                iloop+=1
        #print 'thedict',dict_for_coadd

    for key,vals in dict_for_coadd.items():
        out_table.add_row((np.mean(vals['time']),np.mean(vals['flux']),np.sqrt(np.sum(vals['fluxerr']*vals['fluxerr']))/np.sqrt(float(len(vals))),vals['band'][0],vals['zp'][0],vals['zpsys'][0]))

    return out_table

#filenames=['SuperNova_4.189756_-1.082474_1000_0.01_0.0999999.pkl','SuperNova_4.189756_-1.082474_1000_0.1_0.2.pkl']
#filenames=['SuperNova_4.189756_-1.082474_100_0.4_0.5.pkl']
#filenames=['SuperNova_4.189756_-1.082474_50000_0.4_0.5.pkl']
#filenames=['SuperNova_4.189756_-1.082474_500_0.3_0.4.pkl','SuperNova_4.189756_-1.082474_500_0.2_0.3.pkl','SuperNova_4.189756_-1.082474_500_0.1_0.2.pkl','SuperNova_4.189756_-1.082474_500_0.4_0.5.pkl']
#filenames=['SuperNova_4.189756_-1.082474_500_0.3_0.4.pkl','SuperNova_4.189756_-1.082474_500_0.4_0.5.pkl']

parser = OptionParser()
#parser.add_option("-r", "--ra", type="float", default=0.01, help="filter [%default]")
#parser.add_option("-d", "--dec", type="float", default=0.5, help="filter [%default]")
parser.add_option("-N", "--nevts", type="int", default=10, help="filter [%default]")
parser.add_option("-m", "--model", type="string", default='salt2-extended', help="filter [%default]")
#parser.add_option("-l", "--nlcpoints", type="int", default=5, help="filter [%default]")
parser.add_option("-f", "--fieldtype", type="string", default="WFD", help="filter [%default]")
parser.add_option("-n", "--fieldnum", type="int", default=309, help="filter [%default]")
parser.add_option("-s", "--season", type="int", default=-1, help="filter [%default]")

opts, args = parser.parse_args()

#idlist=opts.fieldtype+'_'+str(opts.fieldnum)+'_'+str(opts.ra)+'_'+str(opts.dec)+'_'+str(opts.nevts)+'_'+opts.model
idlist=opts.fieldtype+'_'+str(opts.fieldnum)+'_'+str(opts.nevts)+'_'+opts.model

if opts.season > -1:
    idlist+='_season_'+str(opts.season)

listname='List_'+idlist+'_minion_test.dat'

filenames=[]

opsim_release='minion_1016'
thedir='Sim_'+opsim_release
outdir='Control_'+opsim_release

for line in open(listname, 'r').readlines():
    filenames.append(thedir+'/'+line.strip())

pkl_file_res = open(outdir+'/output_'+idlist+'.pkl','wb')

#pkl_file_res = open('output_'+idlist+'.pkl','wb')

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


#print len(all_obs),len(all_obs[0])

"""
test=pkl.load(open('SuperNova_4.189756_-1.082474_10000_1.pkl','rb'))
for i,val in enumerate(test):
    print 'ayo',i
"""

vals=['z','t0','c','x1','status','fit_status']
bands= ['u','g','r','i','z','y']

"""
dust = sncosmo.OD94Dust()
fitted_model=sncosmo.Model(source='salt2-extended', effects=[dust, dust],
                    effect_names=['host', 'mw'],
                    effect_frames=['rest', 'obs'])
"""

 # Read throughputs
"""
transmission=Throughputs()

# Register LSST band pass (system) in sncosmo

for filtre in bands:
    band=sncosmo.Bandpass(transmission.lsst_system[filtre].wavelen, transmission.lsst_system[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm)
    sncosmo.registry.register(band)
"""

zdist=[]
zdist_phase=[]
zdist_not_fitted=[]
zdist_fitted_noqual=[]
zdist_fitted=[]
mbvals=[]
#mbvals_sim=[]
zdist_noobs=[]
zcrash=[]
chisq=[]
delta={}
nmeas_band={}
param_sim={}
param_fit={}



for band in bands:
    nmeas_band[band]=[]

for val in ['t0','x1','c','z']:
    delta[val]=[]
    param_sim[val]=[]
    param_fit[val]=[]

thetype=[]
for val in ['t0','x1','c','z']:
        thetype.append((val+'_sim',np.float))

thetype.append(('status',np.dtype('a15')))

for vv in ['','_coadd']:
    thetype.append(('fit_status'+vv,np.dtype('a15')))
    thetype.append(('fit_status'+vv+'_opsim',np.dtype('a15')))
    thetype.append(('chisq'+vv,np.float))
    thetype.append(('ndof'+vv,np.float))
    thetype.append(('mbfit'+vv,np.float))
    thetype.append(('chisq'+vv+'_opsim',np.float))
    thetype.append(('ndof'+vv+'_opsim',np.float))
    thetype.append(('mbfit'+vv+'_opsim',np.float))
    for val in ['t0','x1','c','z']:
        thetype.append((val+'_fit'+vv,np.float))
        thetype.append((val+'_fit'+vv+'_opsim',np.float))
        thetype.append((val+'_fit_error'+vv,np.float))
        thetype.append((val+'_fit_error'+vv+'_opsim',np.float))

for band in bands:
    thetype.append(('nmeas_'+band,np.int))
    thetype.append(('nmeas_before_T0_'+band,np.int))
    thetype.append(('nmeas_after_T0_'+band,np.int))
    thetype.append(('nmeas_coadd_'+band,np.int))
    thetype.append(('nmeas_coadd_before_T0_'+band,np.int))
    thetype.append(('nmeas_coadd_after_T0_'+band,np.int))
    thetype.append(('nmeas_'+band+'_opsim',np.int))
    thetype.append(('nmeas_before_T0_'+band+'_opsim',np.int))
    thetype.append(('nmeas_after_T0_'+band+'_opsim',np.int))
    thetype.append(('nmeas_coadd_'+band+'_opsim',np.int))
    thetype.append(('nmeas_coadd_before_T0_'+band+'_opsim',np.int))
    thetype.append(('nmeas_coadd_after_T0_'+band+'_opsim',np.int))


tab_resu=np.zeros((60,1),dtype=[type for type in thetype])
nfitted=-1

nobs_fit={}
for band in bands:
    nobs_fit[band]=[]


sn_num=0
thedict={}
for oob in all_obs:
    for i,obj in enumerate(oob):
        res_str=''
        """
        if obj['Obs_LC'] != None:
            print len(obj['Obs_LC'])
        """
        #print i,obj['sncosmo_fitted']
        nfitted+=1

        if len(tab_resu) <= nfitted:
            tab_resu=np.resize(tab_resu,(len(tab_resu)+100,1))

        #tab_resu['t0_sim'][nfitted]=obj['t0']
        for vval in ['t0','c','x1','z']:
            tab_resu[vval+'_sim'][nfitted]=obj[vval]

        
        #print 'mes observations',obj['observations'],obj['observations'].dtype
        tab_resu['status'][nfitted]=obj['status']

        """
        if obj['status']=='killed':
            print 'hello',obj['status'],obj['fit']
        """
        #print 'global status',obj['status']
        dict_fit=obj['fit']
        #print dict_fit
        if dict_fit is not None:
            #print 'hello',dict_fit.keys()
            
            for val in ['error_calc','error_opsim','error_coadd_calc','error_coadd_opsim']:
            #for val in ['error_coadd_calc']:
                addit=''
                addita=''
                if val.count('opsim') > 0:
                    addit='_opsim'
                if val.count('coadd') > 0:
                    addita='_coadd'

                dict_tag=dict_fit[val]        
                #print 'fit status',dict_tag['fit_status']
                tab_resu['fit_status'+addita+addit][nfitted]=dict_tag['fit_status']

                table_for_fit=dict_tag['table_for_fit']
                #print 'hello table for fit',table_for_fit

                for band in bands:
                    sel=table_for_fit[np.where(table_for_fit['band']=='LSST::'+band)]
                    tab_resu['nmeas'+addita+'_'+band+addit][nfitted]=len(sel[np.where(np.logical_and(sel['flux']/sel['fluxerr']>5.,sel['flux']>0.))])

                if dict_tag['fit_status'] == 'ok':
                    #print 'mbfit',dict_tag['mbfit']
                    tab_resu['mbfit'+addita+addit][nfitted]=dict_tag['mbfit']
                    resfit=dict_tag['sncosmo_res']
                    corr={}
                    for i,pal in enumerate(dict_tag['sncosmo_res']['vparam_names']):
                        #print i,val
                        corr[pal]=i

                    #print 'alors resfit',resfit
                    tab_resu['chisq'+addita+addit][nfitted]=resfit.chisq
                    tab_resu['ndof'+addita+addit][nfitted]=resfit.ndof
                    for vval in ['t0','c','x1','z']:
                        tab_resu[vval+'_fit'+addita+addit][nfitted]=dict_tag['sncosmo_fitted'][vval] 
                        if resfit['covariance'] is not None:
                            tab_resu[vval+'_fit_error'+addita+addit][nfitted]=np.sqrt(resfit['covariance'][corr[vval]][corr[vval]])
                        else:
                            tab_resu[vval+'_fit_error'+addita+addit][nfitted]=-999.
                    for band in bands:
                        sel=table_for_fit[np.where(table_for_fit['band']=='LSST::'+band)]
                        sel_for_fit=sel[np.where(np.logical_and(sel['flux']/sel['fluxerr']>5.,sel['flux']>0.))]
                        sel_before=sel_for_fit[np.where(sel_for_fit['time']-dict_tag['sncosmo_fitted']['t0']<=0.)]
                        sel_after=sel_for_fit[np.where(sel_for_fit['time']-dict_tag['sncosmo_fitted']['t0']>0.)]
                        tab_resu['nmeas'+addita+'_before_T0_'+band+addit][nfitted]=len(sel_before)
                        tab_resu['nmeas'+addita+'_after_T0_'+band+addit][nfitted]=len(sel_after)
                        #print 'mesures','nmeas'+addita+'_before_T0_'+band+addit,'nmeas'+addita+'_after_T0_'+band+addit,len(sel_before),len(sel_after)
                        """
                        sel_for_coadd=sel
                        sel_before_coadd=sel_for_coadd[np.where(sel_for_coadd['time']-dict_tag['sncosmo_fitted']['t0']<=0.)]
                        sel_after_coadd=sel_for_coadd[np.where(sel_for_coadd['time']-dict_tag['sncosmo_fitted']['t0']>0.)]
                        tab_resu['nmeas_coadd_'+band+addit][nfitted]=len(sel_for_coadd)
                        tab_resu['nmeas_coadd_before_T0_'+band+addit][nfitted]=len(sel_before_coadd)
                        tab_resu['nmeas_coadd_after_T0_'+band+addit][nfitted]=len(sel_after_coadd)
                        """

                    #print 'resfit',resfit
                    #print 'sncosmo_fitted',thedict['sncosmo_fitted']
                    #print 'alors pal',val
                    if val == 'error_coadd_calc':

                        #print 'passed SN',tab_resu['ndof'+addita+addit][nfitted]
                        if tab_resu['ndof'+addita+addit][nfitted] > 0.:

                            thechisq=tab_resu['chisq'+addita+addit][nfitted]/tab_resu['ndof'+addita+addit][nfitted]

                            mbfit=tab_resu['mbfit'+addita+addit][nfitted]
                    

                            if resfit['covariance'] is not None:
                                err_mb=2.5*np.sqrt(resfit['covariance'][2][2])/(dict_tag['sncosmo_fitted']['x0']*np.log(10))

                                sn_num+=1
                                id_SN='SN_'+str(sn_num)
                                print 'supernova',id_SN

                                thedict[id_SN]={}
                                
                                #print dict_tag['sncosmo_res']

                                

                                thedict[id_SN]['host.zcmb']=tab_resu['z_fit'+addita+addit][nfitted]
                                if corr.has_key('z'):
                                    thedict[id_SN]['host.zhelio.err']=np.sqrt(dict_tag['sncosmo_res']['covariance'][corr['z']][corr['z']])
                                else:
                                    thedict[id_SN]['host.zhelio.err']=0.

                                thedict[id_SN]['salt2.RestFrameMag_0_B']=mbfit
                                thedict[id_SN]['salt2.RestFrameMag_0_B.err']=err_mb
                                thedict[id_SN]['salt2.X0']=dict_tag['sncosmo_fitted']['x0']
                                thedict[id_SN]['salt2.X1']=dict_tag['sncosmo_fitted']['x1']
                                thedict[id_SN]['salt2.Color']=dict_tag['sncosmo_fitted']['c']
                                thedict[id_SN]['salt2.CovX1X1']=dict_tag['sncosmo_res']['covariance'][corr['x1']][corr['x1']]
                                thedict[id_SN]['salt2.CovColorColor']=dict_tag['sncosmo_res']['covariance'][corr['c']][corr['c']]
                                thedict[id_SN]['salt2.CovX0X1']=dict_tag['sncosmo_res']['covariance'][corr['x0']][corr['x1']]
                                thedict[id_SN]['salt2.CovColorX0']=dict_tag['sncosmo_res']['covariance'][corr['c']][corr['x0']]
                                thedict[id_SN]['salt2.CovColorX1']=dict_tag['sncosmo_res']['covariance'][corr['c']][corr['x1']]
                                thedict[id_SN]['chisq']=thechisq
                                for band in bands:
                                    thedict[id_SN]['nmeas_'+band]=tab_resu['nmeas'+addita+'_'+band+addit][nfitted]
                                    thedict[id_SN]['nmeas_before_T0_'+band]=tab_resu['nmeas'+addita+'_before_T0_'+band+addit][nfitted]
                                    thedict[id_SN]['nmeas_after_T0_'+band]=tab_resu['nmeas'+addita+'_after_T0_'+band+addit][nfitted]

        if obj['status']=='killed':
            print 'hello',obj['status'],obj['fit'],tab_resu['fit_status'][nfitted]

pkl_file_new = open(outdir+'/SuperNova_'+idlist+'.pkl','wb')
pkl.dump(thedict, pkl_file_new)
pkl_file_new.close()


tab_resu=np.resize(tab_resu,(nfitted+1,1))
pkl.dump(tab_resu, pkl_file_res)

pkl_file_res.close()

print tab_resu.dtype
print tab_resu
