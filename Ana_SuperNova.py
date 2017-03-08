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
"""

class Ana_MyMetrics:

    def __init__(self,idlist):

        self.opsim_release='minion_1016'
        self.thedir='Sim_'+self.opsim_release
        self.outdir='Control_'+self.opsim_release
        self.idlist=idlist
        self.params=['t0','c','x1','z']
        self.bands= ['u','g','r','i','z','y']
        self.obs=self.Get_obs('List_'+idlist+'_minion.dat')

        self.tab_resu=self.Define_Table()
        self.num=-1
        self.sn_num=-1
        self.thedict_SN={}

        self.Fill_Tables()

        self.Save_Results()

    def Get_obs(self,listname):

        filenames=[]

        for line in open(listname, 'r').readlines():
            #filenames.append(self.thedir+'/'+line.strip())
            filenames.append(line.strip())

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

        return all_obs


    def Define_Table(self):

        thetype=[]
        for val in ['t0','x1','c','z']:
            thetype.append((val+'_sim',np.float))

        thetype.append(('status',np.dtype('a15')))
        thetype.append(('sn_type',np.dtype('a15')))
        thetype.append(('sn_model',np.dtype('a15')))
        thetype.append(('sn_version',np.dtype('a15')))

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

        for band in self.bands:
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


        return np.zeros((60,1),dtype=[type for type in thetype])
        

    def Fill_Tables(self):

        for oob in self.obs:
            for i,obj in enumerate(oob):

               self.num+=1
               if len(self.tab_resu) <= self.num:
                   self.tab_resu=np.resize(self.tab_resu,(len(self.tab_resu)+100,1))

               self.tab_resu['status'][self.num]=obj['status']
               
               if obj['status'] != "killed" and obj['status'] != 'No obs in [T0-30;T0+50]':
                   #print 'alors',obj['status'],obj.keys()
                   for vvval in ['sn_type','sn_model','sn_version']:
                       self.tab_resu[vvval][self.num]=obj[vvval]
               

               for vval in self.params:
                   self.tab_resu[vval+'_sim'][self.num]=obj[vval]


               dict_fit=obj['fit']

               if dict_fit is not None:
                   self.Fill_Fit_Results(dict_fit)
                   self.Fill_Supernovae(dict_fit)


               #if self.num > 200.:
                #   break

    def Get_Table(self):
        return np.resize(self.tab_resu,(self.num+1,1))


    def Fill_Fit_Results(self,dict_fit):

        #for val in ['error_calc','error_opsim','error_coadd_calc','error_coadd_opsim']:
        for val in ['error_coadd_calc','error_coadd_opsim']:
            addit=''
            addita=''
            if val.count('opsim') > 0:
                addit='_opsim'
            if val.count('coadd') > 0:
                addita='_coadd'

            dict_tag=dict_fit[val]

            self.tab_resu['fit_status'+addita+addit][self.num]=dict_tag['fit_status']

            table_for_fit=dict_tag['table_for_fit']

            for band in self.bands:
                sel=table_for_fit[np.where(table_for_fit['band']=='LSST::'+band)]
                self.tab_resu['nmeas'+addita+'_'+band+addit][self.num]=len(sel[np.where(np.logical_and(sel['flux']/sel['fluxerr']>5.,sel['flux']>0.))])
               
            #print 'there we are',val,dict_tag['fit_status']
            if dict_tag['fit_status'] == 'ok':
                    #print 'mbfit',dict_tag['mbfit']

                for band in self.bands:
                    sel=table_for_fit[np.where(table_for_fit['band']=='LSST::'+band)]
                    sel_for_fit=sel[np.where(np.logical_and(sel['flux']/sel['fluxerr']>5.,sel['flux']>0.))]
                    sel_before=sel_for_fit[np.where(sel_for_fit['time']-dict_tag['sncosmo_fitted']['t0']<=0.)]
                    sel_after=sel_for_fit[np.where(sel_for_fit['time']-dict_tag['sncosmo_fitted']['t0']>0.)]
                    self.tab_resu['nmeas'+addita+'_before_T0_'+band+addit][self.num]=len(sel_before)
                    self.tab_resu['nmeas'+addita+'_after_T0_'+band+addit][self.num]=len(sel_after)

                self.tab_resu['mbfit'+addita+addit][self.num]=dict_tag['mbfit']
                resfit=dict_tag['sncosmo_res']

                corr={}
                for i,pal in enumerate(dict_tag['sncosmo_res']['vparam_names']):
                    corr[pal]=i

                self.tab_resu['chisq'+addita+addit][self.num]=resfit.chisq
                self.tab_resu['ndof'+addita+addit][self.num]=resfit.ndof
                
                for vval in self.params:
                    self.tab_resu[vval+'_fit'+addita+addit][self.num]=dict_tag['sncosmo_fitted'][vval] 
                    if resfit['covariance'] is not None:
                        self.tab_resu[vval+'_fit_error'+addita+addit][self.num]=np.sqrt(resfit['covariance'][corr[vval]][corr[vval]])
                    else:
                        self.tab_resu[vval+'_fit_error'+addita+addit][self.num]=-999.


    def Fill_Supernovae(self, dict_fit):

        totag='error_coadd_calc'

        addit=''
        addita=''
        if totag.count('opsim') > 0:
            addit='_opsim'
        if totag.count('coadd') > 0:
            addita='_coadd'

        dict_tag=dict_fit[totag] 


        if dict_tag['fit_status'] == 'ok':
            resfit=dict_tag['sncosmo_res']
            
            if resfit.ndof >0:
                thechisq=resfit.chisq/resfit.ndof

                mbfit=dict_tag['mbfit']

                if resfit['covariance'] is not None:

                    self.sn_num+=1
                    id_SN='SN_'+str(self.sn_num)
                  
                    corr={}
                    for i,pal in enumerate(dict_tag['sncosmo_res']['vparam_names']):
                        corr[pal]=i

                    err_mb=2.5*np.sqrt(resfit['covariance'][2][2])/(dict_tag['sncosmo_fitted']['x0']*np.log(10))

                    self.thedict_SN[id_SN]={}
                    
                    self.thedict_SN[id_SN]['host.zcmb']=dict_tag['sncosmo_fitted']['z']
                    if corr.has_key('z'):
                        self.thedict_SN[id_SN]['host.zhelio.err']=np.sqrt(dict_tag['sncosmo_res']['covariance'][corr['z']][corr['z']])
                    else:
                        self.thedict_SN[id_SN]['host.zhelio.err']=0.

                    self.thedict_SN[id_SN]['salt2.RestFrameMag_0_B']=mbfit
                    self.thedict_SN[id_SN]['salt2.RestFrameMag_0_B.err']=err_mb
                    self.thedict_SN[id_SN]['salt2.X0']=dict_tag['sncosmo_fitted']['x0']
                    self.thedict_SN[id_SN]['salt2.X1']=dict_tag['sncosmo_fitted']['x1']
                    self.thedict_SN[id_SN]['salt2.Color']=dict_tag['sncosmo_fitted']['c']
                    self.thedict_SN[id_SN]['salt2.CovX1X1']=dict_tag['sncosmo_res']['covariance'][corr['x1']][corr['x1']]
                    self.thedict_SN[id_SN]['salt2.CovColorColor']=dict_tag['sncosmo_res']['covariance'][corr['c']][corr['c']]
                    self.thedict_SN[id_SN]['salt2.CovX0X1']=dict_tag['sncosmo_res']['covariance'][corr['x0']][corr['x1']]
                    self.thedict_SN[id_SN]['salt2.CovColorX0']=dict_tag['sncosmo_res']['covariance'][corr['c']][corr['x0']]
                    self.thedict_SN[id_SN]['salt2.CovColorX1']=dict_tag['sncosmo_res']['covariance'][corr['c']][corr['x1']]
                    self.thedict_SN[id_SN]['chisq']=thechisq
                    
                    table_for_fit=dict_tag['table_for_fit']

                    for band in self.bands:
                        sel=table_for_fit[np.where(table_for_fit['band']=='LSST::'+band)]
                        self.thedict_SN[id_SN]['nmeas_'+band]=len(sel[np.where(np.logical_and(sel['flux']/sel['fluxerr']>5.,sel['flux']>0.))])
                        sel_for_fit=sel[np.where(np.logical_and(sel['flux']/sel['fluxerr']>5.,sel['flux']>0.))]
                        sel_before=sel_for_fit[np.where(sel_for_fit['time']-dict_tag['sncosmo_fitted']['t0']<=0.)]
                        sel_after=sel_for_fit[np.where(sel_for_fit['time']-dict_tag['sncosmo_fitted']['t0']>0.)]
                        self.thedict_SN[id_SN]['nmeas_before_T0_'+band]=len(sel_before)
                        self.thedict_SN[id_SN]['nmeas_after_T0_'+band]=len(sel_after)


    def Save_Results(self):

        pkl_file_res = open(self.outdir+'/output_'+self.idlist+'.pkl','wb')
        pkl.dump(self.Get_Table(), pkl_file_res)
        pkl_file_res.close()


        pkl_file_new = open(self.outdir+'/SuperNova_'+self.idlist+'.pkl','wb')
        pkl.dump(self.thedict_SN, pkl_file_new)
        pkl_file_new.close()


parser = OptionParser()
#parser.add_option("-r", "--ra", type="float", default=0.01, help="filter [%default]")
#parser.add_option("-d", "--dec", type="float", default=0.5, help="filter [%default]")
parser.add_option("-N", "--nevts", type="int", default=10, help="filter [%default]")
#parser.add_option("-m", "--model", type="string", default='salt2-extended', help="filter [%default]")
#parser.add_option("-l", "--nlcpoints", type="int", default=5, help="filter [%default]")
parser.add_option("-f", "--fieldname", type="string", default="WFD", help="filter [%default]")
parser.add_option("-n", "--fieldid", type="int", default=309, help="filter [%default]")
parser.add_option("-s", "--season", type="int", default=-1, help="filter [%default]")
parser.add_option("-t", "--sntype", type="string", default="Ia", help="filter [%default]")
parser.add_option("-r", "--rolling", type="int", default=0, help="filter [%default]")

opts, args = parser.parse_args()

#idlist=opts.fieldtype+'_'+str(opts.fieldnum)+'_'+str(opts.ra)+'_'+str(opts.dec)+'_'+str(opts.nevts)+'_'+opts.model
addit=''
if opts.rolling == 1:
    addit='Rolling_'

idlist=opts.sntype+'_'+addit+opts.fieldname+'_'+str(opts.fieldid)+'_'+str(opts.nevts)

if opts.season > -1:
    idlist+='_season_'+str(opts.season)

print 'Looking at',idlist
myana=Ana_MyMetrics(idlist)

