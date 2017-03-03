import numpy as np
import pickle as pkl
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
from optparse import OptionParser
from matplotlib import cm
import matplotlib as mpl
from matplotlib.mlab import griddata

def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def Select(thedict,nsel_before,nsel_after,addit='',coadding=False,before_after=True):
    coadd=''
    before=''
    after=''
    if coadding:
        coadd='_coadd'
    if before_after:
        before='_before_T0'
        after='_after_T0'


    sela=thedict[np.where(np.logical_and(thedict['status']=='try_fit',thedict['fit_status'+coadd]=='ok'))]
    print 'here pal','fit_status'+coadd+addit,len(sela)
    sela=sela[np.where(sela['ndof'+coadd+addit]>0.)]
    #sela=sela[np.where(sela['chisq'+addit]/sela['ndof'+addit]<50.)]
    #sela=sela[np.where(sela['c_fit_error'+coadd]>-10.)]
    for band in ['g','r','i']:
            #sela=sela[np.where(sela['nmeas_'+band+addit]>=nsel)]
        #print 'what','nmeas'+coadd+before+'_'+band+addit,'nmeas'+coadd+after+'_'+band+addit
        sela=sela[np.where(sela['nmeas'+coadd+before+'_'+band+addit]>=nsel_before)]
        sela=sela[np.where(sela['nmeas'+coadd+after+'_'+band+addit]>=nsel_after)]
       

    return sela

def Select_new(thedict,dict_combi,addit='',coadding=False,before_after=True):
    coadd=''
    before=''
    after=''
    if coadding:
        coadd='_coadd'
    if before_after:
        before='_before_T0'
        after='_after_T0'


    sela=thedict[np.where(np.logical_and(thedict['status']=='try_fit',thedict['fit_status'+coadd]=='ok'))]
    print 'here pal','fit_status'+coadd+addit,len(sela)
    sela=sela[np.where(sela['ndof'+coadd+addit]>0.)]
    #sela=sela[np.where(sela['chisq'+addit]/sela['ndof'+addit]<50.)]
    #sela=sela[np.where(sela['c_fit_error'+coadd]>-10.)]
   
    selb=sela

    band='r' 
    selb=selb[np.where(selb['nmeas'+coadd+before+'_'+band+addit]>=dict_combi[band][0])]
    selb=selb[np.where(selb['nmeas'+coadd+after+'_'+band+addit]>=dict_combi[band][1])]

    band_one='g'
    band_two='i'

    selb=selb[np.where(np.logical_or(np.logical_and(selb['nmeas'+coadd+before+'_'+band_one+addit]>=dict_combi[band_one][0],selb['nmeas'+coadd+after+'_'+band_one+addit]>=dict_combi[band_one][1]),np.logical_and(selb['nmeas'+coadd+before+'_'+band_two+addit]>=dict_combi[band_two][0],selb['nmeas'+coadd+after+'_'+band_two+addit]>=dict_combi[band_two][1])))]

    return selb

def Select_multiple(thedict,nsel_before,nsel_after,addit='',coadding=True,before_after=True):
    coadd=''
    before=''
    after=''
    if coadding:
        coadd='_coadd'
    if before_after:
        before='_before_T0'
        after='_after_T0'

    sela=thedict
    for band in ['g','r','i']:
        sela=sela[np.where(sela['nmeas'+coadd+before+'_'+band+addit]>=nsel_before[band])]
        sela=sela[np.where(sela['nmeas'+coadd+after+'_'+band+addit]>=nsel_after[band])]

    return sela

def Select_band(thedict,band,nsel_before,nsel_after,addit='',coadding=True,before_after=True):
    
    coadd=''
    before=''
    after=''
    if coadding:
        coadd='_coadd'
    if before_after:
        before='_before_T0'
        after='_after_T0'

    sela=thedict
    print 'there man',len(thedict),coadd,before,after,addit,nsel_before,nsel_after
    sela=sela[np.where(sela['nmeas'+coadd+before+'_'+band+addit]>=nsel_before)]
    sela=sela[np.where(sela['nmeas'+coadd+after+'_'+band+addit]>=nsel_after)]

    return sela


def Histo_ratio(thedict, nsel_before,nsel_after,what,addit='',coadd=False,before_after=True):

    num_bins=24
    range=[0.0,0.8]
    
    hista, bin_edgesa = np.histogram(thedict['z_sim'],bins=num_bins,range=range)

    print 'nums',len(thedict['z_sim'])
    if what == 'Detection efficiency':
        #test=thedict[np.where(np.logical_and(thedict['z_sim']>=0.1,thedict['z_sim']<0.2))]
        sela=Select(thedict,nsel_before,nsel_after,addit,coadd,before_after)
        #print 'hello pal ',len(sela)

    if what == 'Detection efficiency IR':
        sela=thedict[np.where(thedict['fit_status'+addit]=='ok')]
        sela=sela[np.where(sela['ndof'+addit]>0.)]
        sela=sela[np.where(sela['chisq'+addit]/sela['ndof'+addit]<1.)]
        for band in ['i','z','y']:
            sela=sela[np.where(sela['nmeas_'+band]>=nsel)]

    if what == 'Missing obs SNR>5':
        totag='fit_status'
        if coadd:
            totag+='_coadd'
        totag+=addit
        sela=thedict[np.where(np.logical_and(thedict['status']=='try_fit',thedict[totag]=='Noobs'))]

    if what == 'crash fit':
        totag='fit_status'
        if coadd:
            totag+='_coadd'
        totag+=addit
        sela=thedict[np.where(np.logical_and(thedict['status']=='try_fit',thedict[totag]=='crash'))]
        #print 'crash fit',totag,len(sela)

    if what == 'No obs between T0-30 and T0+50':
        #sela=thedict[np.where(thedict['status']=='No obs in [T0-30;T0+50]')]
        sela=thedict[np.where(thedict['status']=='No obs in [T0-3')]
        #print 'alors T0?',len(sela)
    if what == 'Phase rejection':
        #sela=thedict[np.where(thedict['status']=='No obs in [T0-30;T0+50]')]
        sela=thedict[np.where(thedict['status']=='killed')]

    histb, bin_edgesb = np.histogram(sela['z_sim'],bins=num_bins,range=range)
    bin_center = (bin_edgesa[:-1] + bin_edgesa[1:]) / 2
    print 'hello',bin_center,hista,np.sum(hista),len(bin_edgesa),len(hista),len(histb),len(bin_edgesb)

    return bin_center,hista,histb

def Histo_ratio_new(thedict, dict_combi,what,addit='',coadd=False,before_after=True):

    num_bins=24
    range=[0.0,0.8]
    
    hista, bin_edgesa = np.histogram(thedict['z_sim'],bins=num_bins,range=range)

    print 'nums',len(thedict['z_sim'])
    if what == 'Detection efficiency':
        #test=thedict[np.where(np.logical_and(thedict['z_sim']>=0.1,thedict['z_sim']<0.2))]
        sela=Select_new(thedict,dict_combi,addit,coadd,before_after)
        #print 'hello pal ',len(sela)

    if what == 'Detection efficiency IR':
        sela=thedict[np.where(thedict['fit_status'+addit]=='ok')]
        sela=sela[np.where(sela['ndof'+addit]>0.)]
        sela=sela[np.where(sela['chisq'+addit]/sela['ndof'+addit]<1.)]
        for band in ['i','z','y']:
            sela=sela[np.where(sela['nmeas_'+band]>=nsel)]

    if what == 'Missing obs SNR>5':
        totag='fit_status'
        if coadd:
            totag+='_coadd'
        totag+=addit
        sela=thedict[np.where(np.logical_and(thedict['status']=='try_fit',thedict[totag]=='Noobs'))]

    if what == 'crash fit':
        totag='fit_status'
        if coadd:
            totag+='_coadd'
        totag+=addit
        sela=thedict[np.where(np.logical_and(thedict['status']=='try_fit',thedict[totag]=='crash'))]
        #print 'crash fit',totag,len(sela)

    if what == 'No obs between T0-30 and T0+50':
        #sela=thedict[np.where(thedict['status']=='No obs in [T0-30;T0+50]')]
        sela=thedict[np.where(thedict['status']=='No obs in [T0-3')]
        #print 'alors T0?',len(sela)
    if what == 'Phase rejection':
        #sela=thedict[np.where(thedict['status']=='No obs in [T0-30;T0+50]')]
        sela=thedict[np.where(thedict['status']=='killed')]

    histb, bin_edgesb = np.histogram(sela['z_sim'],bins=num_bins,range=range)
    bin_center = (bin_edgesa[:-1] + bin_edgesa[1:]) / 2
    print 'hello',bin_center,hista,np.sum(hista),len(bin_edgesa),len(hista),len(histb),len(bin_edgesb)

    return bin_center,hista,histb

def Histo_error(hista,histb):
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
    
    return ratio, ratio_err,norm,norm_err

thedict={}


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

#fichname='output_DD_290_6.097944_-1.10516_3000_salt2-extended_4points.pkl'

opts, args = parser.parse_args()

addit=''
if opts.rolling == 1:
    addit='Rolling_'
ids=opts.sntype+"_"+addit+opts.fieldname+'_'+str(opts.fieldid)+'_'+str(opts.nevts)
if opts.season > -1:
    ids+='_season_'+str(opts.season)

fichname='output_'+ids+'.pkl'

opsim_release='minion_1016'
imin=4
imax=5

for i in range(imin,imax):
    pkl_file = open('Control_'+opsim_release+'/'+fichname,'rb')
    thedict[i]=pkl.load(pkl_file)

Plot_eff_z=False
Plot_eff_z_per_season=False
Plot_eff_z_all=True
Plot_noobs_z=False
Plot_Deltas=False
Plot_sim_vs_fit=False
Plot_sim=True
Plot_nmeas_vs_z=False
Gime_Stats=False
Loop_Combis=False
Loop_Combis_new=False
Draw_Chisquare=False
Read_Fill=False

tab_resu=thedict[i]

#print tab_resu.dtype


fontsize=15.

combis=[(2,1)]
#combis=[(0,0),(2,1),(2,2)]
#combis=[(3,3),(4,4),(5,5)]

combi_dict={}
combi_dict['g']=(2,2)
combi_dict['r']=(2,2)
combi_dict['i']=(2,2)

myfmt=['--','-.',':']
mymarkers=['o','v','s']
coadd=True
before_after=True

#print tab_resu


if Plot_eff_z_per_season:

    figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))

    tot_label=[]
    
    dict_season={}
    
    sfile=open('Seasons/Seasons_'+addit+opts.fieldname+'_'+str(opts.fieldid)+'.txt', 'r')
    
    for line in sfile.readlines():
        sais=int(line.split(' ')[0])
        TSeason_min=float(line.split(' ')[1])
        TSeason_max=float(line.split(' ')[2])
        dict_season[sais]=(TSeason_min,TSeason_max)

    print 'Nombre de season',len(dict_season)
    colors=['black','blue','green','yellow','red','black','blue','green','yellow','red']
    for key, val in dict_season.items():

        tab_resub=tab_resu[np.where(np.logical_and(tab_resu['t0_sim']>= val[0],tab_resu['t0_sim']<= val[1]))]
        print 'season',key,len(tab_resu)
    
        for i,comb in enumerate(combis):
            bin_center, hista, histb=Histo_ratio_new(tab_resub,combi_dict,'Detection efficiency','',coadd,before_after)
            ratio, ratio_err,norm,norm_err=Histo_error(hista,histb)
            if not before_after:
                ll=r'$N_{points}^{LC}$='+str(comb[0])
            else:
                ll=r'$N_{points}^{LC}$=('+str(comb[0])+','+str(comb[1])+')'
            ll='Season '+str(key+1)
            marker='.'
            if key > 4:
                marker='v'
            tot_label.append(axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=marker, mfc=colors[key], mec=colors[key], ms=8, linestyle=myfmt[i],color='k',label=ll))

        """
        thecol='g'
        axca = axc.twinx()
    
        bin_center, hista, histb=Histo_ratio(tab_resu,-1,-1,'Missing obs SNR>5','',coadd,before_after)
        ratio, ratio_err,norm,norm_err=Histo_error(hista,histb)
        ll=r'Missing obs SNR>5'
        tot_label.append(axca.errorbar(bin_center,ratio, yerr=ratio_err,marker='>', mfc='k', mec='k', ms=5, linestyle='-',color='k',label=ll))
    

        bin_centerb, histab, histbb=Histo_ratio(tab_resu,-1,-1,'No obs between T0-30 and T0+50','',coadd,before_after)
        ratiob, ratio_errb,normb,norm_errb=Histo_error(histab,histbb)
        ll=r'No obs in [T0-30,T0+50]'
        tot_label.append(axca.errorbar(bin_centerb,ratiob, yerr=ratio_errb,marker='<', mfc=thecol, mec=thecol, ms=5, linestyle='--',color='g',label=ll))
    
        bin_centerc, histac, histbc=Histo_ratio(tab_resu,-1,-1,'Phase rejection','',coadd,before_after)
        ratioc, ratio_errc,normc,norm_errc=Histo_error(histac,histbc)
        ll=r'Phase rejection'
        tot_label.append(axca.errorbar(bin_centerc,ratioc, yerr=ratio_errc,marker='s', mfc=thecol, mec=thecol, ms=5, linestyle='--',color='g',label=ll))
        
        bin_centerd, histad, histbd=Histo_ratio(tab_resu,-1,-1,'crash fit','',coadd,before_after)
        ratiod, ratio_errd,normd,norm_errd=Histo_error(histad,histbd)
        ll=r'crash fit'
        tot_label.append(axca.errorbar(bin_centerd,ratiod, yerr=ratio_errd,marker='|', mfc='m', mec='m', ms=5, linestyle='-',color='m',label=ll))
        """

    axc.set_xlabel(r'z ',{'fontsize': 20.})
    axc.set_ylabel(r'Efficiency (per z-bin)',{'fontsize': 20.})
    #axc.legend(tot_label,loc='lower right')

    labs = [l.get_label() for l in tot_label]
   
    axc.set_xlim(0, 0.8)
    axc.set_ylim(0, 1.)
    #axc.legend(tot_label, labs, ncol=2,bbox_to_anchor=(1.0, 0.16),loc='lower right',prop={'size':10},frameon=False)
    axc.legend(tot_label, labs, ncol=2,loc='upper left',prop={'size':12},frameon=False)
    #axca.legend()
    axc.set_title('Field type:'+opts.fieldname+' Field ID: '+str(opts.fieldid))

if Plot_eff_z_all:

    figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))

    tot_label=[]
    
    #before_after=False
    #for i in range(2,5):
    for i,comb in enumerate(combis):
        bin_center, hista, histb=Histo_ratio_new(tab_resu,combi_dict,'Detection efficiency','',coadd,before_after)
        ratio, ratio_err,norm,norm_err=Histo_error(hista,histb)
        if not before_after:
            ll=r'$N_{points}^{LC}$='+str(comb[0])
        else:
            ll=r'$N_{points}^{LC}$=('+str(comb[0])+','+str(comb[1])+')'
        ll=r'$N_{points}^{LC}$=(2,2) in (g,r) or (r,i)'
        tot_label.append(axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=mymarkers[i], mfc='red', mec='red', ms=5, linestyle=myfmt[i],color='k',label=ll))

    for i,comb in enumerate(combis):
        bin_center, hista, histb=Histo_ratio_new(tab_resu,combi_dict,'Detection efficiency','_opsim',coadd,before_after)
        ratio, ratio_err,norm,norm_err=Histo_error(hista,histb)
        if not before_after:
            ll=r'$N_{points}^{LC}$ (opsim)='+str(comb[0])
        else:
            ll=r'$N_{points}^{LC}$ (opsim)=('+str(comb[0])+','+str(comb[1])+')'
        ll=r'$N_{points}^{LC}$ (opsim) =(2,2) in (g,r) or (r,i)'
        tot_label.append(axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=mymarkers[i], mfc='b', mec='b', ms=5, linestyle=myfmt[i],color='k',label=ll))
        """
    for i in range(3,6):
        bin_center, hista, histb=Histo_ratio(tab_resu,i,'Detection efficiency IR')
        ratio, ratio_err,norm,norm_err=Histo_error(hista,histb)
        ll=r'$N_{points}^{LC}='+str(i)+'$'
        axc.errorbar(bin_center,ratio, yerr=ratio_err,marker='.', mfc='blue', mec='blue', ms=15, linestyle=myfmt[i-3],color='k', label=ll)
        """
    
    thecol='g'
    axca = axc.twinx()
    
    bin_center, hista, histb=Histo_ratio(tab_resu,-1,-1,'Missing obs SNR>5','',coadd,before_after)
    ratio, ratio_err,norm,norm_err=Histo_error(hista,histb)
    ll=r'Missing obs SNR>5'
    tot_label.append(axca.errorbar(bin_center,ratio, yerr=ratio_err,marker='>', mfc='k', mec='k', ms=5, linestyle='-',color='k',label=ll))
    
    bin_centera, histaa, histba=Histo_ratio(tab_resu,-1,-1,'Missing obs SNR>5','_opsim',coadd,before_after)
    ratioa, ratio_erra,norma,norm_erra=Histo_error(histaa,histba)
    ll=r'Missing obs SNR>5 (Opsim)'
    tot_label.append(axca.errorbar(bin_centera,ratioa, yerr=ratio_erra,marker='o', mfc='k', mec='k', ms=5, linestyle='-',color='k',label=ll))

    bin_centerb, histab, histbb=Histo_ratio(tab_resu,-1,-1,'No obs between T0-30 and T0+50','',coadd,before_after)
    ratiob, ratio_errb,normb,norm_errb=Histo_error(histab,histbb)
    ll=r'No obs in [T0-30,T0+50]'
    tot_label.append(axca.errorbar(bin_centerb,ratiob, yerr=ratio_errb,marker='<', mfc=thecol, mec=thecol, ms=5, linestyle='--',color='g',label=ll))
    
    bin_centerc, histac, histbc=Histo_ratio(tab_resu,-1,-1,'Phase rejection','',coadd,before_after)
    ratioc, ratio_errc,normc,norm_errc=Histo_error(histac,histbc)
    ll=r'Phase rejection'
    tot_label.append(axca.errorbar(bin_centerc,ratioc, yerr=ratio_errc,marker='s', mfc=thecol, mec=thecol, ms=5, linestyle='--',color='g',label=ll))
 
    bin_centerd, histad, histbd=Histo_ratio(tab_resu,-1,-1,'crash fit','',coadd,before_after)
    ratiod, ratio_errd,normd,norm_errd=Histo_error(histad,histbd)
    ll=r'crash fit'
    tot_label.append(axca.errorbar(bin_centerd,ratiod, yerr=ratio_errd,marker='|', mfc='m', mec='m', ms=5, linestyle='-',color='m',label=ll))
    
    bin_centere, histe, histbe=Histo_ratio(tab_resu,-1,-1,'crash fit','_opsim',coadd,before_after)
    ratioe, ratio_erre,norme,norm_erre=Histo_error(histe,histbe)
    ll=r'crash fit (Opsim)'
    tot_label.append(axca.errorbar(bin_centere,ratioe, yerr=ratio_erre,marker='s', mfc='m', mec='m', ms=5, linestyle='-',color='m',label=ll))

    for tl in axca.get_yticklabels():
        tl.set_color(thecol)

    axc.set_xlabel(r'z ',{'fontsize': 20.})
    axc.set_ylabel(r'Efficiency (per z-bin)',{'fontsize': 20.})
    #axc.legend(tot_label,loc='lower right')

    labs = [l.get_label() for l in tot_label]
   
    axc.set_xlim(0, 0.8)
    axc.set_ylim(0, 1.)
    #axc.legend(tot_label, labs, ncol=2,bbox_to_anchor=(1.0, 0.16),loc='lower right',prop={'size':10},frameon=False)
    axc.legend(tot_label, labs, ncol=2,loc='upper left',prop={'size':12},frameon=False)
    #axca.legend()
    axc.set_title('Field type:'+opts.fieldname+' Field ID: '+str(opts.fieldid))

if Plot_eff_z:

    figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))


    myfmt=['--','-.',':']
    mymarkers=['o','v','s']
    for i in range(3,6):
        bin_center, hista, histb=Histo_ratio(tab_resu,i,'Detection efficiency')
        ratio, ratio_err,norm,norm_err=Histo_error(hista,histb)
        ll=r'$N_{points}^{LC}='+str(i)+'$'
        axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=mymarkers[i-3], mfc='red', mec='red', ms=5, linestyle=myfmt[i-3],color='k', label=ll)

    for i in range(3,6):
        bin_center, hista, histb=Histo_ratio(tab_resu,i,'Detection efficiency','_opsim')
        ratio, ratio_err,norm,norm_err=Histo_error(hista,histb)
        ll=r'$N_{points}^{LC} (opsim)='+str(i)+'$'
        axc.errorbar(bin_center,ratio, yerr=ratio_err,marker=mymarkers[i-3], mfc='b', mec='b', ms=5, linestyle=myfmt[i-3],color='k', label=ll)

        """
    for i in range(3,6):
        bin_center, hista, histb=Histo_ratio(tab_resu,i,'Detection efficiency IR')
        ratio, ratio_err,norm,norm_err=Histo_error(hista,histb)
        ll=r'$N_{points}^{LC}='+str(i)+'$'
        axc.errorbar(bin_center,ratio, yerr=ratio_err,marker='.', mfc='blue', mec='blue', ms=15, linestyle=myfmt[i-3],color='k', label=ll)
        """
    axc.set_xlabel(r'z ',{'fontsize': 20.})
    axc.set_ylabel(r'Efficiency (per z-bin)',{'fontsize': 20.})
    axc.legend()

if Plot_noobs_z:
   
    figca, axca = plt.subplots(ncols=1, nrows=1, figsize=(10,9))

    bin_center, hista, histb=Histo_ratio(tab_resu,i,'Missing obs SNR>5')
    ratio, ratio_err,norm,norm_err=Histo_error(hista,histb)
    ll=r'Missing obs SNR>5'
    axca.errorbar(bin_center,ratio, yerr=ratio_err,marker='.', mfc='red', mec='red', ms=15, linestyle='-',color='k', label=ll)

    bin_centerb, histab, histbb=Histo_ratio(tab_resu,i,'No obs between T0-30 and T0+50')
    ratiob, ratio_errb,normb,norm_errb=Histo_error(histab,histbb)
    ll=r'No obs between T0-30 and T0+50'
    axca.errorbar(bin_centerb,ratiob, yerr=ratio_errb,marker='.', mfc='red', mec='red', ms=15, linestyle='--',color='k', label=ll)

    axca.set_xlabel(r'z ',{'fontsize': 20.})
    axca.set_ylabel(r'Efficiency (per z-bin)',{'fontsize': 20.})             
    axca.legend(loc='upper left')

if Plot_sim:
    figbb, axbb = plt.subplots(ncols=2, nrows=2, figsize=(10,9))

    for (j,vals) in enumerate(['t0','z','x1','c']):
        if j==0 or j==2:
            k=0
        else:
            k=1

        axbb[j/2][k].hist(tab_resu[vals+'_sim'],bins=10,histtype='step')
            
        
        axbb[j/2][k].set_xlabel(r''+vals+'_sim',{'fontsize': fontsize})
        axbb[j/2][k].set_ylabel(r'Number of Entries',{'fontsize': fontsize})
        print vals,np.mean(tab_resu[vals+'_sim']),np.std(tab_resu[vals+'_sim'])

if Plot_Deltas:

    figb, axb = plt.subplots(ncols=2, nrows=2, figsize=(10,9))

    cut={}
    cut['t0']=0.5
    cut['z']=0.02
    cut['x1']=0.3
    cut['c']=0.1

    n_standard=2.
    addita='_coadd'
    sela=tab_resu[np.where(np.logical_and(tab_resu['status']=='try_fit',tab_resu['fit_status'+addita]=='ok'))]
    sela=sela[np.where(sela['ndof'+addita]>0.)]
    #sela=sela[np.where(sela['chisq'+addita]/sela['ndof'+addita]<1.)]

    mcol=['k','r','g']

    #before_after=False

    for (j,vals) in enumerate(['t0','z','x1','c']):
        if j==0 or j==2:
            k=0
        else:
            k=1

        for i,comb in enumerate(combis):
            #mysel=[vvals for vvals in thedict[i]['delta_'+vals] if np.abs(vvals)<cut[vals]]
          
            selb=Select(sela,comb[0],comb[1],'',coadd,before_after)

            mean=np.mean(selb[vals+'_sim']-selb[vals+'_fit'+addita])
            rms=np.std(selb[vals+'_sim']-selb[vals+'_fit'+addita])
            #mysel=selb[np.where(np.abs(selb[vals+'_sim']-selb[vals+'_fit'])<n_standard*rms)]
            #mysel=selb[np.where(np.abs(selb[vals+'_sim']-selb[vals+'_fit'])<cut[vals])]
            mysel=selb
            #print 'hello',mysel.dtype.names
            sel=(mysel[vals+'_sim']-mysel[vals+'_fit'+addita])/mysel[vals+'_fit_error'+addita]
            #axb[j/2][k].hist(thedict[i]['delta_'+vals],bins=40,histtype='step')
            
            #ll=r'$N_{points}^{LC}='+str(i+1)+str("%.2f %(mean)")+' '+str(rms)+'$'
            if not before_after:
                ll=r'$N_{points}^{LC}=$'+str(comb[0])+' : '+str(round(mean,4))+'$\pm$'+str(round(rms,4))
            else:
               ll=r'$N_{points}^{LC}=$('+str(comb[0])+','+str(comb[1])+') : '+str(round(mean,4))+'$\pm$'+str(round(rms,4)) 

            print 'hello',i,vals,mean,rms,len(sel)
            if len(sel) > 0:
                axb[j/2][k].hist(sel,bins=40,histtype='step',label=ll,color=mcol[i])
            #print 'hello',i,vals,np.mean(thedict[i]['delta_'+vals]),np.std(thedict[i]['delta_'+vals])
            
            axb[j/2][k].set_xlabel(r'$\Delta$'+vals,{'fontsize': fontsize})
            axb[j/2][k].set_ylabel(r'Number of Entries',{'fontsize': fontsize})

        axb[j/2][k].legend(loc='upper left',prop={'size':8})

    figb.suptitle('Field type:'+opts.fieldname+' Field ID: '+str(opts.fieldid))
    figc, axc = plt.subplots(ncols=2, nrows=2, figsize=(10,9))


    for (j,vals) in enumerate(['t0','z','x1','c']):
        if j==0 or j==2:
            k=0
        else:
            k=1

        #colors=['k.','r.','g.']
        for i,comb in enumerate(combis):
            #mysel=[vvals for vvals in thedict[i]['delta_'+vals] if np.abs(vvals)<cut[vals]]
          
            
            selb=Select(sela,comb[0],comb[1],'',coadd,before_after)
            sel=selb[vals+'_sim']-selb[vals+'_fit'+addita]
            #axc[j/2][k].hist(thedict[i]['delta_'+vals],bins=40,histtype='step')
            
            #ll=r'$N_{points}^{LC}='+str(i+1)+str("%.2f %(mean)")+' '+str(rms)+'$'
            ll=r'$N_{points}^{LC}=$'+str(i)+' : '+str(round(mean,4))+'$\pm$'+str(round(rms,4))
            print 'hello',i,vals,mean,rms
            #axc[j/2][k].hist(sel,bins=40,histtype='step',label=ll)
            #axc[j/2][k].plot(selb['chisq']/selb['ndof'],sel,'k.',label=ll)
            axc[j/2][k].plot(selb[vals+'_sim'],sel,mcol[i]+'.',label=ll)
            #print 'hello',i,vals,np.mean(thedict[i]['delta_'+vals]),np.std(thedict[i]['delta_'+vals])
            
            axc[j/2][k].set_ylabel(r'$\Delta$'+vals,{'fontsize': fontsize})
            #axc[j/2][k].set_xlabel(r'$\chi^2$',{'fontsize': fontsize})
            axc[j/2][k].set_xlabel(r''+vals+'_sim',{'fontsize': fontsize})

        #axc[j/2][k].legend(loc='upper right',prop={'size':8})       
    figc.suptitle('Field type:'+opts.fieldname+' Field ID: '+str(opts.fieldid))

if Plot_nmeas_vs_z:

    sela=tab_resu[np.where(np.logical_and(tab_resu['status']=='try_fit',tab_resu['fit_status']=='ok'))]
    sela=sela[np.where(sela['ndof']>0.)]
    sela=sela[np.where(sela['chisq']/sela['ndof']<1.)]

    
    figa, axa = plt.subplots(ncols=2, nrows=3, figsize=(10,9))

    for j,band in enumerate(['u','g','r','i','z','y']):
        if j<2:
            k=0
        if j>= 2 and j < 4:
            k=1
        if j>=4:
            k=2
            
        axa[k][j%2].plot(sela['z_sim'],sela['nmeas_'+band],'k.')
        axa[k][j%2].plot(sela['z_sim'],sela['nmeas_coadd_'+band],'r.')
        
        axa[k][j%2].set_ylabel(r'$N_{meas} (SNR>5)$',{'fontsize': fontsize})
            #axc[j/2][k].set_xlabel(r'$\chi^2$',{'fontsize': fontsize})
        axa[k][j%2].set_xlabel(r'$z_sim$',{'fontsize': fontsize})
    
    figb, axb = plt.subplots(ncols=2, nrows=3, figsize=(10,9))

    for j,band in enumerate(['u','g','r','i','z','y']):
        if j<2:
            k=0
        if j>= 2 and j < 4:
            k=1
        if j>=4:
            k=2
            
        axb[k][j%2].plot(sela['z_sim'],sela['nmeas_'+band]-sela['nmeas_coadd_'+band],'k.')
        



if Plot_sim_vs_fit:

    sela=tab_resu[np.where(np.logical_and(tab_resu['status']=='try_fit',tab_resu['fit_status']=='ok'))]
    sela=sela[np.where(sela['ndof']>0.)]
    sela=sela[np.where(sela['chisq']/sela['ndof']<1.)]
    figa, axa = plt.subplots(ncols=2, nrows=2, figsize=(10,9))
    colors=['k.','b.','r.']
    
    for (j,vals) in enumerate(['t0','z','x1','c']):
        if j==0 or j==2:
            k=0
        else:
            k=1
        for i in range(3,6):
            for band in ['g','r','i']:
                sela=sela[np.where(sela['nmeas_'+band]>i)]


    #axa.plot(thedict[i]['zdist_fitted'],thedict[i]['delta_z'],'k.')
            axa[j/2][k].plot(sela[vals+'_sim'],sela[vals+'_fit'],colors[i-3])

        fontsize=15.    
        axa[j/2][k].set_ylabel(r''+vals+'_fit',{'fontsize': fontsize})
        axa[j/2][k].set_xlabel(r''+vals+'_sim',{'fontsize': fontsize})  


"""
i=4
chitot=[]
ztot=[]

sel=thedict[i]['tab_fitted'][np.where(thedict[i]['tab_fitted']['ndof']>0.)]

plt.plot(sel['chisq']/sel['ndof'])
"""

resu=thedict[imin]

#print resu['status']
ntot=float(len(resu))
print 'nentries',ntot
"""
for status in [('killed','Phase rejection'),('No obs','No obs between T0-30 and T0+50'),('Not fitted','Missing obs SNR >5'),('fitted','Fitted - no quality'),('crashed','crash fit')]:
    if status[0] != 'crashed':
        sela=resu[np.where(resu['status']==status[0])]
        print status[1],len(sela),round(100.*float(len(sela))/ntot,2),'%'
    else:
        if status[0]!= 'Not fitted':
            sela=resu[np.where(resu['fit_status']==status[0])]
            print status[1],len(sela),round(100.*float(len(sela))/ntot,2),'%'
        else:
            sela=resu[np.where(np.logical_and(resu['status']==status[0],resu['fit_status']=='unknown'))]
            print status[1],len(sela),round(100.*float(len(sela))/ntot,2),'%'

    if status[0] == 'fitted':
        sela=resu[np.where(resu['status']==status[0])]
        sela=sela[np.where(sela['ndof']>0.)]
        sela=sela[np.where(sela['chisq']/sela['ndof']<1.)]
        for num in range(3,6):
            for band in ['g','r','i']:
                sela=sela[np.where(sela['nmeas_'+band]>=num)]
            print 'Fitted with qual - npoints=',num,len(sela),round(100.*float(len(sela))/ntot,2),'%'
"""

if Gime_Stats:

    coadding=''
    before=''
    after=''
    if coadd:
        coadding='_coadd'
    if before_after:
        before='_before_T0'
        after='_after_T0'

    #mysel=resu[np.where(np.logical_and(resu['z_sim']>=0.1,resu['z_sim']<0.2))]
    mysel=resu
    ntot=len(mysel)
    for status in [(('killed',''),'Phase rejection'),(('No obs in [T0-3',''),'No obs in range[T0-30,T0+50]'),(('try_fit','Noobs'),'Missing obs SNR >5'),(('try_fit','ok'),'Fitted - no quality'),(('try_fit','crash'),'fit crashed')]:
        #sela=resu[np.where(resu['status']==status[0][0])]
        sela=mysel[np.where(mysel['status']==status[0][0])]
        if status[0][0] == 'try_fit':
            sela=sela[np.where(sela['fit_status'+coadding]==status[0][1])]
            
            print status[1],len(sela),round(100.*float(len(sela))/ntot,2),'%'
        if status[0][0] == 'killed':
            print status[1],len(sela),round(100.*float(len(sela))/ntot,2),'%'

        if status[0][0] == 'try_fit' and status[0][1]== 'ok':
           
            selb=sela[np.where(sela['fit_status'+coadding]==status[0][1])]
            selb=selb[np.where(selb['ndof'+coadding]>0.)]
            #selb=selb[np.where(selb['chisq']/selb['ndof']<1.)]
            #print 'alors?',len(sela),len(selb)
            for comb in combis:
                selc=Select(selb,comb[0],comb[1],'',coadd,before_after)
                print 'Fitted with qual - npoints=',comb[0],comb[1],len(selc),round(100.*float(len(selc))/ntot,2),'%'

if Draw_Chisquare:

    
    coadding=''
    before=''
    after=''
    if coadd:
        coadding='_coadd'
    if before_after:
        before='_before_T0'
        after='_after_T0'
    selb=resu[np.where(np.logical_and(resu['status']=='try_fit',resu['fit_status'+coadding]=='ok'))]
    selb=selb[np.where(selb['ndof'+coadding]>0.)]
    print 'alooooooo',len(selb)
    figbb, axbb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    tot_label=[]
    for zcut in np.arange(0.,0.6,0.1):
        zmax=zcut+0.1
        selc=selb[np.where(np.logical_and(selb['z_sim']>=zcut,selb['z_sim']<zmax))]
    #print selc['chisq']/selc['ndof']
        print 'bbbbbb',len(selc),np.min(selc['chisq'+coadding]/selc['ndof'+coadding]),np.median(selc['chisq'+coadding]/selc['ndof'+coadding])
        selc=selc[np.where(selc['chisq'+coadding]/selc['ndof'+coadding]<1.e17)]
        
        ll=r''+str(zcut)+'<z<'+str(zmax)
        tot_label.append(axbb.hist(selc['chisq'+coadding]/selc['ndof'+coadding],bins=10000000,normed=1,histtype='step',cumulative=True,label=ll))
    
    axbb.set_xlim(1.e-5, 50.)
    
    axbb.set_xlabel(r'$\chi^2/NDof$',{'fontsize': 20.})
    axbb.set_ylabel(r'Cumulative',{'fontsize': 20.})
    axbb.set_xscale('log')
    axbb.legend()
"""
labs = [l.get_label() for l in tot_label]
axbb.legend(tot_label, labs, ncol=2,loc='upper left',prop={'size':12},frameon=False)
"""
if Loop_Combis_new:

    zmin=0.4
    zmax=0.5

    band_one='r'
    band_two='i'

    selb=resu[np.where(np.logical_and(resu['z_sim']>=zmin,resu['z_sim']<zmax))]

    n_norm=float(len(selb))
    selb=selb[np.where(np.logical_and(selb['status']=='try_fit',selb['fit_status_coadd']=='ok'))]
    selb=selb[np.where(selb['ndof_coadd']>0.)]

    combi_baft=[]
    
    for i in range(2,5):
        for j in range(2,5): 
            combi_baft.append((i,j))

    dict_effi={}

    for combi in combi_baft:
        selc=Select_band(selb,band_one,combi[0],combi[1])
        print 'first sel',len(selc),combi[0],combi[1]
        dict_effi[combi]={}
        for combib in combi_baft:
            dict_effi[combi][combib]={}
            dict_effi[combi][combib]['chisq']=[]
            dict_effi[combi][combib]['effi']=[]
            dict_effi[combi][combib]['error_effi']=[]

            seld=Select_band(selc,band_two,combib[0],combib[1])
            print 'second sel',len(seld)

            print 'hello',np.median(seld['chisq_coadd']/seld['ndof_coadd']),np.mean(seld['chisq_coadd']/seld['ndof_coadd']),np.max(seld['chisq_coadd']/seld['ndof_coadd'])

            max_chisquare=np.max(seld['chisq_coadd']/seld['ndof_coadd'])  

            chisq_step=max_chisquare/50.
            for chisq in np.arange(0.,max_chisquare,chisq_step):
                selcc=seld[np.where(seld['chisq_coadd']/seld['ndof_coadd']<chisq)]
                dict_effi[combi][combib]['chisq'].append(chisq)
                effi=float(len(selcc))/n_norm
                dict_effi[combi][combib]['effi'].append(effi)
                dict_effi[combi][combib]['error_effi'].append(np.sqrt(effi*(1.-effi)/n_norm))


    markers=['.',',','o','v','^','<','>','*','p'] 
    colors=['b','g','r','c','m','y','k','0.75','0.5']

    icombi_one=-1
    tot_label=[]
    figc, axc = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    idd=-1
    for key,val in dict_effi.items():
        print 'combi one',key
        icombi_two=-1
        for keyb,valb in val.items():
            print 'combi2',keyb,valb['effi']
            icombi_two+=1
            #axc.errorbar(valb['chisq'],valb['effi'],yerr=valb['error_effi'],fmt='-',color = 'r')
        #break
            print 'icombi2',icombi_two
            draw=False
            for effi in valb['effi']:
                if effi>0.5:
                    draw=True
                    break
            #tot_label.append(axc.errorbar(valb['chisq'],valb['effi'],yerr=valb['error_effi'],marker=markers[icombi_two], mfc=colors[icombi_two], mec=colors[icombi_two], ms=8, linestyle='--',color='k',label=ll))
            #draw=True
            if draw:
                idd+=1
                ll=str(key)+' '+str(keyb)
                tot_label.append(axc.errorbar(valb['chisq'],valb['effi'],yerr=valb['error_effi'],marker=markers[icombi_two], mfc=colors[idd], mec=colors[idd], ms=8, linestyle='--',color='k',label=ll))
  
    labs = [l.get_label() for l in tot_label]
   
    #axc.set_xlim(0, 0.8)
    axc.set_ylim(0.45, 0.65)
    #axc.legend(tot_label, labs, ncol=2,bbox_to_anchor=(1.0, 0.16),loc='lower right',prop={'size':10},frameon=False)
    axc.legend(tot_label, labs, ncol=2,loc='upper left',prop={'size':12},frameon=False)
                

    axc.set_xlabel(r'$\chi^2/Ndof$',{'fontsize': 20.})
    axc.set_ylabel(r'Efficiency ('+str(zmin)+' < z < '+str(zmax)+')',{'fontsize': 20.})

    axc.set_title('Field type:'+opts.fieldname+' Field ID: '+str(opts.fieldid))
    plt.show()
            

if Loop_Combis:

    selb=resu[np.where(np.logical_and(resu['z_sim']>=0.1,resu['z_sim']<0.2))]

    n_norm=float(len(selb))
    selb=selb[np.where(np.logical_and(selb['status']=='try_fit',selb['fit_status']=='ok'))]
    selb=selb[np.where(selb['ndof']>0.)]
    #selc=selb[np.where(selb['chisq']/selb['ndof']<50.)]
    selc=selb
    print 'effi',n_norm,float(len(selc))/float(len(selb))

    combi_baft=[]
    for i in range(1,5):
       for j in range(1,5): 
           combi_baft.append((i,j))

    print 'hello',np.median(selc['chisq']/selc['ndof']),np.mean(selc['chisq']/selc['ndof']),np.max(selc['chisq']/selc['ndof'])

    max_chisquare=np.max(selc['chisq']/selc['ndof'])

    mytype=[('i_combi',np.int), ('chisq_cut', np.float),('n_before_g', np.int),('n_after_g', np.int),('n_before_r', np.int),('n_after_r', np.int),('n_before_i', np.int),('n_after_i', np.int),('n_norm', np.float),('n_select', np.float),('effi', np.float)]

    tab_resu=np.zeros((60,1),dtype=[type for type in mytype])
    
    ievt=-1
    ichisq=-1
    chisq_step=max_chisquare/100.
    for chisq in np.arange(0.,max_chisquare,chisq_step):
        ichisq+=1
        i_combi=-1
        selcc=selc
        selcc=selcc[np.where(selcc['chisq']/selcc['ndof']<chisq)]
        nsel_before={}
        nsel_after={}
        for combi_g in combi_baft:
            nsel_before['g']=combi_g[0]
            nsel_after['g']=combi_g[1]
            for combi_r in combi_baft:
                nsel_before['r']=combi_r[0]
                nsel_after['r']=combi_r[1]
                for combi_i in combi_baft:
                    nsel_before['i']=combi_i[0]
                    nsel_after['i']=combi_i[1]
                    
                    i_combi+=1
                    seld=Select_multiple(selcc,nsel_before,nsel_after)
                    #print 'hh',chisq,ichisq,i_combi,nsel_before['g'],nsel_after['g'],nsel_before['r'],nsel_after['r'],nsel_before['i'],nsel_after['i'],len(seld),n_norm,round(100.*float(len(seld))/n_norm,2),'%'
                
                    ievt+=1

                    if len(tab_resu) <= ievt:
                        tab_resu=np.resize(tab_resu,(len(tab_resu)+100,1))

                    tab_resu['i_combi'][ievt]=i_combi
                    tab_resu['chisq_cut'][ievt]=chisq
                    tab_resu['n_norm'][ievt]=n_norm
                    tab_resu['n_select'][ievt]=float(len(seld))
                    tab_resu['effi'][ievt]=100.*float(len(seld))/n_norm
                    for band in ['g','r','i']:
                        tab_resu['n_before_'+band][ievt]=nsel_before[band]
                        tab_resu['n_after_'+band][ievt]=nsel_after[band]
                   
    tab_resu=np.resize(tab_resu,(ievt+1,1))
    
    print tab_resu,len(tab_resu),ichisq,np.max(tab_resu['i_combi']),len(combi_baft),combi_baft

    """
    figbb, axbb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    
    axbb.hist(selc['chisq']/selc['ndof'],bins=1000,normed=1,histtype='step',cumulative=True)
    axbb.set_xlim([0.,0.999*max_chisquare])

    """

    """
    fig = plt.figure(figsize=plt.figaspect(0.5))

    
    chisq=[]
    icombi=[]
    effi=[]

    for i in range(len(tab_resu)):
        chisq.append(tab_resu['chisq_cut'][i][0])
        icombi.append(tab_resu['i_combi'][i][0])
        effi.append(tab_resu['effi'][i][0])


    nx=len(chisq)
    ny=len(icombi)
    
    X, Y = np.meshgrid(tab_resu['chisq_cut'], tab_resu['i_combi'])
    Z = np.array(effi).reshape(nx,ny).T

    im = axbb.pcolormesh(X,Y,Z)
    figbb.colorbar(im, ax=axbb)
    """ 
    chisq=[]
    icombi=[]
    effi=[]

    for i in range(len(tab_resu)):
        chisq.append(tab_resu['chisq_cut'][i][0])
        icombi.append(tab_resu['i_combi'][i][0])
        effi.append(tab_resu['effi'][i][0])

    
    figbb, axbb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    
    xedges=[]
    yedges=[]
    xedgesa=[]
    yedgesa=[]
    for i in np.arange(0.,max_chisquare,chisq_step):
        xedges.append(i)

    #for i in range(0,9*9*9):
    for i in range(0,16*16*16):
        yedges.append(i) 

    #plt.colorbar()
    H, x, y, p =axbb.hist2d(chisq,icombi, bins=(xedges, yedges),weights=effi)
    H = H.T
    X, Y = np.meshgrid(xedges,yedges)
    pp=axbb.pcolormesh(X, Y, H,cmap=cm.jet, vmin=0., vmax=100.)

    plt.colorbar(pp)
    imax=np.argmax(tab_resu['effi'], axis=0)
    print 'Done',imax,np.max(tab_resu['effi']),tab_resu[imax]
    """
    nx=len(chisq)
    ny=len(icombi)
    
    
    X, Y = np.meshgrid(tab_resu['chisq_cut'], tab_resu['i_combi'])
    Z=effi
    Z =Z.reshape(nx,ny)
    """





if Read_Fill:

    
  
    #var_to_dump=['t0_sim', 'x1_sim', 'c_sim', 'z_sim', 'status', 'fit_status_coadd', 'chisq_coadd', 'ndof_coadd', 't0_fit_coadd', 't0_fit_error_coadd', 'x1_fit_coadd', 'x1_fit_error_coadd', 'c_fit_coadd', 'c_fit_error_coadd', 'z_fit_coadd', 'z_fit_error_coadd', 'nmeas_coadd_u', 'nmeas_coadd_before_T0_u', 'nmeas_coadd_after_T0_u','nmeas_coadd_g', 'nmeas_coadd_before_T0_g', 'nmeas_coadd_after_T0_g', 'nmeas_coadd_r', 'nmeas_coadd_before_T0_r', 'nmeas_coadd_after_T0_r', 'nmeas_coadd_i', 'nmeas_coadd_before_T0_i', 'nmeas_coadd_after_T0_i', 'nmeas_coadd_z', 'nmeas_coadd_before_T0_z', 'nmeas_coadd_after_T0_z', 'nmeas_coadd_y', 'nmeas_coadd_before_T0_y', 'nmeas_coadd_after_T0_y']

    var_to_dump=['t0_sim', 'x1_sim', 'c_sim', 'z_sim', 'chisq_coadd', 'ndof_coadd', 't0_fit_coadd', 't0_fit_error_coadd', 'x1_fit_coadd', 'x1_fit_error_coadd', 'c_fit_coadd', 'c_fit_error_coadd', 'z_fit_coadd', 'z_fit_error_coadd', 'nmeas_coadd_u', 'nmeas_coadd_before_T0_u', 'nmeas_coadd_after_T0_u','nmeas_coadd_g', 'nmeas_coadd_before_T0_g', 'nmeas_coadd_after_T0_g', 'nmeas_coadd_r', 'nmeas_coadd_before_T0_r', 'nmeas_coadd_after_T0_r', 'nmeas_coadd_i', 'nmeas_coadd_before_T0_i', 'nmeas_coadd_after_T0_i', 'nmeas_coadd_z', 'nmeas_coadd_before_T0_z', 'nmeas_coadd_after_T0_z', 'nmeas_coadd_y', 'nmeas_coadd_before_T0_y', 'nmeas_coadd_after_T0_y','sn_type','sn_model','sn_version']

    mytab=tab_resu[np.where(np.logical_and(tab_resu['status']=='try_fit',tab_resu['fit_status_coadd']=='ok'))]

    names='#'
    for name in var_to_dump:
        names+=name+' '


    outfile= open('Supernovae/Supernovae_'+opts.fieldname+'_'+str(opts.fieldid)+'_'+opts.sntype+'.txt', 'w')
    outfile.write(names+'\n')
    

    for ival in range(len(mytab)):
        toprint=''
        for var in var_to_dump:
            toprint+=str(mytab[var][ival])+' '

        print toprint
        outfile.write(toprint+'\n')
        
    outfile.close()




plt.show()
