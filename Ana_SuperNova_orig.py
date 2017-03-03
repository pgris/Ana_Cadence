import pickle as pkl
from astropy.table import Table
import sncosmo
import numpy as np
import matplotlib.pyplot as plt
from Throughputs import Throughputs
import astropy.units as u
import pickle as pkl
from astropy import (cosmology, units as u, constants as const)
from optparse import OptionParser

"""
filename='SuperNova_4.189756_-1.082474_1000.pkl'
pkl_file = open(filename,'rb')
"""

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
parser.add_option("-l", "--nlcpoints", type="int", default=5, help="filter [%default]")
parser.add_option("-f", "--fieldtype", type="string", default="WFD", help="filter [%default]")
parser.add_option("-n", "--fieldnum", type="int", default=309, help="filter [%default]")

opts, args = parser.parse_args()

#idlist=opts.fieldtype+'_'+str(opts.fieldnum)+'_'+str(opts.ra)+'_'+str(opts.dec)+'_'+str(opts.nevts)+'_'+opts.model
idlist=opts.fieldtype+'_'+str(opts.fieldnum)+'_'+str(opts.nevts)+'_'+opts.model
listname='List_'+idlist+'_minion.dat'

filenames=[]

opsim_release='minion_1016'
thedir='Sim_'+opsim_release
outdir='Control_'+opsim_release

for line in open(listname, 'r').readlines():
    filenames.append(thedir+'/'+line.strip())

pkl_file_res = open(outdir+'/output_'+idlist+'_'+str(opts.nlcpoints)+'points.pkl','wb')

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

transmission=Throughputs()

# Register LSST band pass (system) in sncosmo

for filtre in bands:
    band=sncosmo.Bandpass(transmission.lsst_system[filtre].wavelen, transmission.lsst_system[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm)
    sncosmo.registry.register(band)

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
    thetype.append((val+'_fit',np.float))

thetype.append(('chisq',np.float))
thetype.append(('ndof',np.float))

thetype.append(('status',np.dtype('a15')))
thetype.append(('fit_status',np.dtype('a15')))

for band in bands:
    thetype.append(('nmeas_'+band,np.int))

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
        #print obj
        for val in vals:
            #print 'alors',val,obj[val]
            res_str+=str(obj[val])+' '
        #print res_str
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

        
        print 'mes observations',obj['observations'],obj['observations'].dtype
        tab_resu['status'][nfitted]=obj['status']
        tab_resu['fit_status'][nfitted]=obj['fit_status']


        
        zdist.append(obj['z'])
    
        if obj['status']=='killed':
            zdist_phase.append(obj['z'])
 
        if obj['status']=='No obs':
            zdist_noobs.append(obj['z']) 

        if obj['status']=='Not fitted' and obj['fit_status']=='unknown':
            zdist_not_fitted.append(obj['z'])
        
        if obj['fit_status']=='crashed':
            zcrash.append(obj['z'])

        if obj['status']=='fitted':
            
            zdist_fitted_noqual.append(obj['z'])
            chisq.append((obj['sncosmo_res'].chisq,obj['sncosmo_res'].ndof))
           
            for band in bands:
                sel=obj['Obs_LC_for_fit'][np.where(obj['Obs_LC_for_fit']['band']=='LSST::'+band)]
                tab_resu['nmeas_'+band][nfitted]=len(sel[np.where(np.logical_and(sel['flux']/sel['fluxerr']>5.,sel['flux']>0.))])

            tab_resu['chisq'][nfitted]=obj['sncosmo_res'].chisq
            tab_resu['ndof'][nfitted]=obj['sncosmo_res'].ndof     
            #tab_resu['t0_fit'][nfitted]=obj['sncosmo_fitted']['t0']
            
            for vval in ['t0','c','x1','z']:
                tab_resu[vval+'_fit'][nfitted]=obj['sncosmo_fitted'][vval]

            if obj['sncosmo_res'].ndof >0:
                
                nmeas={}
                thechisq=obj['sncosmo_res'].chisq/obj['sncosmo_res'].ndof
                #print 'chisq pal',thechisq
                for band in bands:
                    sel=obj['Obs_LC_for_fit'][np.where(obj['Obs_LC_for_fit']['band']=='LSST::'+band)]
                    nmeas[band]=len(sel[np.where(np.logical_and(sel['flux']/sel['fluxerr']>5.,sel['flux']>0.))])
                    #print 'selections',band,nmeas[band]
                
                delta['t0'].append(obj['t0']-obj['sncosmo_fitted']['t0'])
                param_sim['t0'].append(obj['t0'])
                param_fit['t0'].append(obj['sncosmo_fitted']['t0'])
                
                for vval in ['c','x1','z']:
                    delta[vval].append(obj[vval]-obj['sncosmo_fitted'][vval])
                    param_sim[vval].append(obj[vval])
                    param_fit[vval].append(obj['sncosmo_fitted'][vval])
     
                #nsel_r=len(obj['Obs_LC_for_fit'][np.where(obj['Obs_LC_for_fit']['band']=='LSST::r')])
                #nsel_i=len(obj['Obs_LC_for_fit'][np.where(obj['Obs_LC_for_fit']['band']=='LSST::i')])
            
                #if nsel_g > 4 and nsel_r > 4 and nsel_i > 4 and thechisq < 0.5:
                if thechisq < 1 and nmeas['g'] > opts.nlcpoints and nmeas['r'] > opts.nlcpoints and nmeas['i'] > opts.nlcpoints:
                    
                    zdist_fitted.append(obj['z'])
                    
                    """
                    fitted_model.set(z=obj['sncosmo_fitted']['z'])
                    fitted_model.set(t0=obj['sncosmo_fitted']['t0'])
                    fitted_model.set(x0=obj['sncosmo_fitted']['x0'])
                    fitted_model.set(x1=obj['sncosmo_fitted']['x1'])
                    fitted_model.set(c=obj['sncosmo_fitted']['c'])
                    fitted_model.set(hostebv=obj['sncosmo_fitted']['hostebv'])
                    fitted_model.set(hostr_v=obj['sncosmo_fitted']['hostr_v'])
                    fitted_model.set(mwebv=obj['sncosmo_fitted']['mwebv'])
                    fitted_model.set(mwr_v=obj['sncosmo_fitted']['mwr_v'])
            
                    fitted_model.set_source_peakabsmag(-19.3,'bessellb', 'vega',cosmo=cosmology.WMAP9)
                    """
                    
                    #print 'hello',fitted_model.param_names,fitted_model.parameters

                    #mb=fitted_model.bandmag('bessellb','vega',obj['sncosmo_fitted']['t0'])
                    
                    mb=obj['mbfit']
                    #mbvals.append(obj['mbsim']-mb)
                    #mbvals_sim.append(obj['mbsim'])
                    

                    if obj['sncosmo_res']['covariance'] is not None:
                        err_mb=2.5*np.sqrt(obj['sncosmo_res']['covariance'][2][2])/(obj['sncosmo_fitted']['x0']*np.log(10))
                        #print 'hello',obj['status'],obj['fit_status'],thechisq,obj['z'],obj['sncosmo_res']['covariance'],obj['sncosmo_res']['covariance'][2][2],'ahah',obj['sncosmo_fitted']['x0']
                        
                   
                    #print 'mB',mb,err_mb,obj['sncosmo_fitted']['z']
                    #thedict[id_SN]['host.zcmb']=obj['sncosmo_fitted']['z']

                        sn_num+=1
                        id_SN='SN_'+str(sn_num)
        
                        thedict[id_SN]={}


                        thedict[id_SN]['host.zcmb']=obj['z']
                        thedict[id_SN]['salt2.RestFrameMag_0_B']=obj['mbfit']
                        thedict[id_SN]['salt2.RestFrameMag_0_B.err']=err_mb
                        thedict[id_SN]['host.zhelio.err']=np.sqrt(obj['sncosmo_res']['covariance'][0][0])
                        thedict[id_SN]['salt2.X0']=obj['sncosmo_fitted']['x0']
                        thedict[id_SN]['salt2.X1']=obj['sncosmo_fitted']['x1']
                        thedict[id_SN]['salt2.Color']=obj['sncosmo_fitted']['c']
                        thedict[id_SN]['salt2.CovX1X1']=obj['sncosmo_res']['covariance'][3][3]
                        thedict[id_SN]['salt2.CovColorColor']=obj['sncosmo_res']['covariance'][4][4]
                        thedict[id_SN]['salt2.CovX0X1']=obj['sncosmo_res']['covariance'][2][3]
                        thedict[id_SN]['salt2.CovColorX0']=obj['sncosmo_res']['covariance'][4][2]
                        thedict[id_SN]['salt2.CovColorX1']=obj['sncosmo_res']['covariance'][4][3]
        
                    else:
                        err_mb=-999
                        print 'Problem here - Fit ok but non covariance matrix at the output'
                        print 'hello',obj['status'],obj['fit_status'],thechisq,obj['z'],obj['sncosmo_res']['covariance'],obj['sncosmo_fitted']
                        #break
                        
        #err_mB=
        """
        print 'mB',mB,err_mB,len(obj['Obs_LC']),obj['sncosmo_res'].chisq/obj['sncosmo_res'].ndof
        print obj['Obs_LC'],obj['sncosmo_fitted']['x0'],fitted_model.get('x0')
       
        tmin = fitted_model.mintime()
        tmax = fitted_model.maxtime()
        tgrid = np.linspace(tmin, tmax, int(tmax - tmin) + 1)

        astro_table=obj['Obs_LC']
        print 'astrotable',astro_table
        for band in ['g','r','i','z','y']:
            sele=astro_table[np.where(astro_table['band']=='LSST::'+band)]
            time=np.sort(sele['time'])
            print 'time',band,time 
            mflux = fitted_model.bandflux('LSST::'+band, time, zp=25, zpsys='ab')
            print mflux

        #sncosmo.plot_lc(obj['Obs_LC'], model=fitted_model,color='k',pulls=False)
        sncosmo.plot_lc(obj['Obs_LC_for_fit'], model=fitted_model,color='k',pulls=False,errors=obj['sncosmo_res'].errors)
     
        plt.show()
        """

    """
    if i > 10:
        break
    """

#print thedict
#print zdist
pkl_file_new = open(outdir+'/SuperNova_'+idlist+'_'+str(opts.nlcpoints)+'points.pkl','wb')
pkl.dump(thedict, pkl_file_new)
pkl_file_new.close()
print 'Stats',len(zdist),len(zdist_phase),len(zdist_not_fitted),len(zdist_fitted_noqual),len(zdist_fitted),len(zcrash),len(zdist_noobs),nfitted

"""
dictout={}

dictout['zdist']=zdist
dictout['zdist_phase']=zdist_phase
dictout['zdist_not_fitted']=zdist_not_fitted
dictout['zdist_fitted_noqual']=zdist_fitted_noqual
dictout['zdist_fitted']=zdist_fitted
dictout['zcrash']=zcrash
dictout['zdist_noobs']=zdist_noobs
dictout['chi2']=chisq 

for vval in ['t0','c','x1','z']:
    dictout['delta_'+vval]=delta[vval]
    dictout[vval+'_sim']=param_sim[vval]
    dictout[vval+'_fit']=param_fit[vval]


for band in bands:
    dictout['nmeas_'+band]=nmeas_band[band]
"""
tab_resu=np.resize(tab_resu,(nfitted+1,1))


#dictout['tab_resu']=tab_resu

pkl.dump(tab_resu, pkl_file_res)

pkl_file_res.close()

stat_dict={}

stat_dict['Nsimul']=float(len(zdist))
stat_dict['Phase rejection']=float(len(zdist_phase))
stat_dict['Missing obs SNR >5']=float(len(zdist_not_fitted))
stat_dict['Fitted - no quality']=float(len(zdist_fitted_noqual))
stat_dict['Fitted - with quality']=float(len(zdist_fitted))
stat_dict['crash fit']=float(len(zcrash))

stat_dict['No obs between T0-30 and T0+50']=float(len(zdist_noobs))

nverif=0
stat_file = open('Stat_'+idlist+'_'+str(opts.nlcpoints)+'points.txt', 'w')

for key,val in stat_dict.items():
    if key!='Nsimul':
        print key,val,100.*(val/stat_dict['Nsimul']),'%'
        towrite=key+' '+str(val)+' '+str(100.*(val/stat_dict['Nsimul']))+'% \n'
        stat_file.write(towrite)
        if key!='Fitted - with quality':
            nverif+=val


if nverif != stat_dict['Nsimul']:
    print 'We have a problem here',nverif,stat_dict['Nsimul']
    stat_file.write('We have a problem here '+str(nverif)+' '+str(stat_dict['Nsimul']))
else:
    print 'Total sum equal to Nsimul -> ok!'
    stat_file.write('Total sum equal to Nsimul -> ok!')

stat_file.close()


"""
fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(10,9))
ax[0].set_xlabel(r'z',{'fontsize': 20.})
ax[0].set_ylabel(r'Number of Events',{'fontsize': 20.})

ax[0].hist(zdist, bins= 64,histtype='step',color='k')
#ax[0].hist(zdist_fitted, bins= 64,histtype='step',color='r')

ax[1].set_xlabel(r'$\chi^{2}/NDF$',{'fontsize': 20.})
ax[1].set_ylabel(r'Number of Events',{'fontsize': 20.})

ax[1].hist(chisq[0]/chisq[1], bins= 64,histtype='step',color='k')

"""
"""
figb, axb = plt.subplots(ncols=4, nrows=1, figsize=(10,9))

#axb[0].plot( nobs_fit['g'],delta['t0'],'k.')
#axb[1].plot( nobs_fit['g'],delta['z'],'k.')
axb[0].hist(delta['t0'],bins=40,color='k')
axb[0].set_xlabel(r'$\Delta$t0',{'fontsize': 20.})
axb[0].set_ylabel(r'Number of Entries',{'fontsize': 20.})
axb[1].hist(delta['z'],bins=40,color='k')
axb[1].set_xlabel(r'$\Delta$ z',{'fontsize': 20.})
axb[1].set_ylabel(r'Number of Entries',{'fontsize': 20.})
axb[2].hist(delta['x1'],bins=40,color='k')
axb[2].set_xlabel(r'$\Delta$ x1',{'fontsize': 20.})
axb[2].set_ylabel(r'Number of Entries',{'fontsize': 20.})
axb[3].hist(delta['c'],bins=40,color='k')
axb[3].set_xlabel(r'$\Delta$ c',{'fontsize': 20.})
axb[3].set_ylabel(r'Number of Entries',{'fontsize': 20.})
figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))

#axc[0].plot(zdist,[a/b for a,b in zip(zdist_fitted,zdist)], bins= 64,histtype='step',color='r')
num_bins=24
range=[0.0,0.8]

hista, bin_edgesa = np.histogram(zdist,bins=num_bins,range=range)
histb, bin_edgesb = np.histogram(zdist_fitted,bins=num_bins,range=range)
bin_center = (bin_edgesa[:-1] + bin_edgesa[1:]) / 2
print 'hello',bin_center,hista,np.sum(hista),len(bin_edgesa),len(hista),len(histb),len(bin_edgesb)

ratio=[]
ratio_err=[]
norm=[]
norm_err=[]

for a,b in zip(histb,hista):
    if b==0:
        ratio.append(a)
        ratio_err.append(0)
    else:
        effi=float(a)/float(b)
        ratio.append(effi)
        ratio_err.append(np.sqrt(1.-effi)*effi/np.sqrt(float(b)))
    eff=float(a)/float(np.sum(hista))
    norm.append(eff)
    norm_err.append(np.sqrt(1.-eff)*eff/np.sqrt(float(np.sum(hista))))
    
print 'hello',histb,hista,norm_err
print 'alors',np.sum(hista),np.sum(histb)
"""
"""
axc[0].errorbar(bin_center,norm, yerr=norm_err,fmt='.')
axc[0].set_xlabel(r'z ',{'fontsize': 20.})
axc[0].set_ylabel(r'Efficiency',{'fontsize': 20.})
"""

"""
axc.errorbar(bin_center,ratio, yerr=ratio_err,marker='s', mfc='red', mec='red', ms=10)
axc.set_xlabel(r'z ',{'fontsize': 20.})
axc.set_ylabel(r'Efficiency (per z-bin)',{'fontsize': 20.})

figd, axd = plt.subplots(ncols=2, nrows=1, figsize=(10,9))

axd[0].plot(zdist_fitted,mbvals,'k.')
axd[1].hist(mbvals,bins=40)

print 'Mean',np.mean(mbvals)
"""


plt.show()
