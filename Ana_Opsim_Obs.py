import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
import math
import healpy as hp
from astropy.table import Table
from optparse import OptionParser
from Parameters import parameters
 

def Get_fiveSigmaDepth_coadd(nmeas,m5):
    
    g=2.3
    alpha=1./g
    Nosq=(0.04*np.power(10,-0.4*m5)-alpha)*np.power(10,-0.4*m5)
    print 'hello',Nosq
    No=float(nmeas)*np.sqrt(Nosq)
    
    delta=alpha*alpha+4.*No*No*0.04/float(nmeas)
    print 'delta',delta
    sol=(-alpha+np.sqrt(delta))/(2*No*No)
    print 'solution',sol
    return 2.5*np.log10(sol)
    


def Get_mean_finSeeing(observations):
    
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
    return np.mean(res)

def Get_coadd(filt):
   
    
    dict_for_coadd={}
    filtc=filt.copy()
    filtc.reshape((filtc.size,1))
    if len(filtc) > 0:
        inum=0
        dict_for_coadd[inum]=np.zeros((0,1),filtc.dtype)
        #print 'timediff',24.*60.*60.*(filtc['time']-filtc['time'][0])
                                
        iloop=0
        #print 'blablabla',dict_for_coadd[inum]
       
        dict_for_coadd[inum]=np.vstack([dict_for_coadd[inum],filtc[iloop]])
                                
        if len(filtc) > 1:
            while iloop < len(filtc)-1:   
                diff_time_sec=24.*60.*60.*(filtc['expMJD'][iloop+1]-filtc['expMJD'][iloop])
                #print 'alors ???',diff_time_sec,inum
                if diff_time_sec > 40.:
                    inum+=1
                    dict_for_coadd[inum]=np.zeros((0,1),filtc.dtype)
                
                dict_for_coadd[inum]=np.vstack([dict_for_coadd[inum],filtc[iloop]])
                    
                iloop+=1
        #print 'thedict',dict_for_coadd

    return dict_for_coadd

def Get_Obs_per_day(season):
 
    Tmin=np.min(season['expMJD'])
    Tmax=np.max(season['expMJD'])

    days_obs={}
    days_no_obs={}
    bands=['u','g','r','i','z','y']

    for band in bands:
        days_no_obs[band]=[]

    season.sort(order='expMJD')

    for time in np.arange(Tmin,Tmax,1.):
        sel_obs=season[np.where(np.logical_and(season['expMJD']>=time,season['expMJD']<time+1.))]
        
        if len(sel_obs) > 0:
            days_obs[time-Tmin]=sel_obs
            mean_time=np.mean(sel_obs['expMJD'])

            for j,band in enumerate(bands):
                selb=sel_obs[np.where(sel_obs['filter']==band)]
                if len(selb)==0:
                    days_no_obs[band].append((time-Tmin,mean_time))
            
    
    return days_obs, days_no_obs


def Copy_Season(season,shift=0.,frac=1.):
    
    season.sort(order='expMJD')
    season_copy=np.zeros((60,1),dtype=season.dtype)
    season_remain=np.zeros((60,1),dtype=season.dtype)

    excluded=[]

    n_tobe_excluded=int((1.-frac)*len(season))

    while len(excluded)< n_tobe_excluded:
        aleat=np.random.randint(0,len(season))
        if aleat not in excluded:
            excluded.append(aleat)

    #print 'excluded',len(excluded),len(season),(1.-frac)*len(season),float(len(excluded))/float(len(season)),excluded
    #print 'alors man',len(season),season
    inum=-1
    inum_remain=-1
    for i in range(len(season)):
        if i not in excluded:
            inum+=1
            if len(season_copy) <= inum:
                season_copy=np.resize(season_copy,(len(season_copy)+50,1))
            for data_name in season.dtype.names:
                season_copy[data_name][inum]=season[data_name][i]
                season_copy['expMJD'][inum]=season['expMJD'][i]-shift
                
        else:
            inum_remain+=1
            if len(season_remain) <= inum_remain:
                season_remain=np.resize(season_remain,(len(season_remain)+50,1))

            for data_name in season.dtype.names:
                season_remain[data_name][inum_remain]=season[data_name][i]
                season_remain['expMJD'][inum_remain]=season['expMJD'][i] 


    season_copy=np.resize(season_copy,(inum+1,1))
    season_remain=np.resize(season_remain,(inum_remain+1,1))

    #print 'final result',len(season_copy),len(season_remain)
    return season_copy,season_remain

def Get_Seasons(filtb):
 
    #Get the seasons of a given set of obs
    #A season is defined by a set of observations 
    # for which the consecutive difference in time is lower than 100 days
    #input : fitc = set of observations
    #output: dict of seasons; key = seson number (starting at 0), val = corresponding set of observations
    
    dict_for_seasons={}
    filtc=filtb.copy()

    filtc.sort(order='expMJD')
        

    if len(filtc) > 0:
        inum=0
        dict_for_seasons[inum]=np.zeros((0,1),dtype=filtc.dtype)
           
                                
        iloop=0
        iinside=0
        dict_for_seasons[inum]=np.vstack([dict_for_seasons[inum],filtc[iloop]])
            
        if len(filtc) > 1:
            while iloop < len(filtc)-1: 
                iinside+=1
                diff_time_days=filtc['expMJD'][iloop+1]-filtc['expMJD'][iloop]
                if diff_time_days > 100.:
                       
                    inum+=1
                    dict_for_seasons[inum]=np.zeros((0,1),dtype=filtc.dtype)
                       

                dict_for_seasons[inum]=np.vstack([dict_for_seasons[inum],filtc[iloop+1]])
   
                iloop+=1
                
        
        return dict_for_seasons


def Get_Seasons_old(filtc):

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

#List=['Observations_DD_290.0_6.097944_-1.10516.pkl','Observations_WFD_309.0_4.189756_-1.082474.pkl']

telSeeing = 0.250 # design goal
opticalDesSeeing = 0.08
cameraSeeing = 0.300
scaleToNeff = 1.16
atmNeffFactor = 1.04
FWHM_Sys_Zenith = np.sqrt(telSeeing**2 + opticalDesSeeing**2 + cameraSeeing**2)
filterWave = {'u': 367.0, 'g': 482.5, 'r': 622.2, 'i': 754.5, 'z': 869.1, 'y': 971.0}

kAtm = {'u': 0.50,
         'g': 0.21,
         'r': 0.13,
         'i': 0.10,
         'z': 0.07,
         'y': 0.18} 
 

msky = {'u': 22.95,
        'g': 22.24,
        'r': 21.20,
        'i': 20.47,
        'z': 19.60,
        'y': 18.63} 

Cm = {'u':22.94,
      'g':24.46,
      'r':24.48,
      'i':24.34,
      'z':24.18,
      'y':23.73}

dCm_infinity = {'u':0.56,
                'g':0.12,
                'r':0.06,
                'i':0.05,
                'z':0.03,
                'y':0.02}



#List=['Observations_DD_290.pkl','Observations_DD_744.pkl','Observations_DD_1427.pkl','Observations_DD_2412.pkl','Observations_DD_2786.pkl']

parser = OptionParser()
parser.add_option("-f", "--fieldname", type="string", default="DD", help="filter [%default]")
parser.add_option("-F", "--fieldid", type="int", default=290, help="filter [%default]")
parser.add_option("-r", "--rolling", type="int", default=0, help="filter [%default]")
parser.add_option("-n", "--nmergers", type="int", default=3, help="filter [%default]")
parser.add_option("-p", "--merge_factor", type="int", default=80, help="merge factor (%) [%default]")

opts, args = parser.parse_args()
#outdir='/sps/lsst/data/dev/pgris/Obs_minion_1016'
outdir='/sps/lsst/data/dev/pgris/Make_Cadence/Obs_minion_1016'
addit=''

prefix='Observations'
if opts.rolling:
    prefix='Rolling_Cadence'
    addit='_'+str(opts.nmergers)+'_'+str(opts.merge_factor)

List=[prefix+'_'+opts.fieldname+'_'+str(opts.fieldid)+addit+'.pkl']

bands=['u','g','r','i','z','y']

thedict={}

fieldid={}
for i,name in enumerate(List):
    pkl_file = open(outdir+'/'+name,'rb')
    thedict[i]=pkl.load(pkl_file)
    #fieldid[i]=np.int(name.split('_')[2].split('.')[0])

print thedict[0].keys(),thedict[0]['dataSlice'].dtype.names

data={}
data['RA']=[]
data['Dec']=[]
data['dust']=[]

"""
figa, axa = plt.subplots(ncols=2, nrows=3, figsize=(18, 12))
cols=['b.','r.']
dict_pos={}

dict_pos['u']=(0,0)
dict_pos['g']=(0,1)
dict_pos['r']=(1,0)
dict_pos['i']=(1,1)
dict_pos['z']=(2,0)
dict_pos['y']=(2,1)

for key,vals in dict_pos.items():

    for i in range(2):
        
        select=thedict[i]['dataSlice'][np.where(thedict[i]['dataSlice']['filter']==key)]
        axa[vals[0]][vals[1]].plot(select['expMJD'],select['airmass'],cols[i])
"""

"""
ival=1
Tmin_obs=np.min(thedict[ival]['dataSlice']['expMJD'])
Tmax_obs=np.max(thedict[ival]['dataSlice']['expMJD'])

n_noobs=0
nvals=50000

for i in range(nvals):
    T0=np.random.uniform(Tmin_obs,Tmax_obs)
    timelow=T0-30
    timehigh=T0+50
    select=thedict[ival]['dataSlice'][np.where(thedict[ival]['dataSlice']['expMJD']>=timelow)]
    select=select[np.where(select['expMJD']<=timehigh)]

    if len(select) < 1.:
        n_noobs+=1

print n_noobs,float(n_noobs)/float(nvals)
"""

Draw_Molleid=False
Check_m5=False
Make_Rolling=False
Gime_Seasons=False
Draw_Seasons=False
Ana_Cadence=True
Ana_Season=False
Dump_in_File=False


if Draw_Molleid:

    for i in range(1):
        print 'dust',i,thedict[i]['ebvofMW']
        data['RA'].append(np.rad2deg(thedict[i]['dataSlice']['fieldRA'][0]))
        data['Dec'].append(np.rad2deg(thedict[i]['dataSlice']['fieldDec'][0]))
        data['dust'].append(thedict[i]['ebvofMW'])


        """
        minda=np.min(thedict[1]['dataSlice']['expMJD'])
        
        axa.plot(thedict[1]['dataSlice']['expMJD'],thedict[1]['dataSlice']['airmass'],'b.')
        
        minda=np.min(thedict[0]['dataSlice']['expMJD'])
        
        axa.plot(thedict[0]['dataSlice']['expMJD'],thedict[0]['dataSlice']['airmass'],'r.')
        """
        
        
        ra = coord.Angle(data['RA'],unit=u.degree)
        ra = ra.wrap_at(180*u.degree)
        dec = coord.Angle(data['Dec'],unit=u.degree)
        
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111, projection="mollweide")
        ax.scatter(ra.radian, dec.radian)
        ax.grid(True)

        """
        # lon_0 is central longitude of projection.
        # resolution = 'c' means use crude resolution coastlines.
        m = Basemap(projection='moll',lon_0=0,resolution='c')
        m.drawcoastlines()
        m.fillcontinents(color='coral',lake_color='aqua')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-90.,120.,30.))
        m.drawmeridians(np.arange(0.,420.,60.))
        m.drawmapboundary(fill_color='aqua')
        plt.title("Mollweide Projection")
        """




        figc, axc = plt.subplots(ncols=2, nrows=2, figsize=(10,9))
        
        for i in range(1):
            #print 'hello',i,fieldid[i]
            axc[0][0].hist(thedict[i]['dataSlice']['airmass'],bins=20,histtype='step')
            axc[0][1].hist(thedict[i]['dataSlice']['filtSkyBrightness'],bins=20,histtype='step')


if Gime_Seasons:
 
    seasons=Get_Seasons(thedict[0]['dataSlice'])

    outputdir='/sps/lsst/data/dev/pgris/Seasons'
    after='Seasons'
    if prefix.count('Rolling'):
        after+='_Rolling'
    outfile= open(outputdir+'/'+after+'_'+opts.fieldname+'_'+str(opts.fieldid)+addit+'.txt', 'w')

    for iseason,season in seasons.items(): 
        print 'min max',iseason,np.min(season['expMJD']),np.max(season['expMJD'])
        theline=str(iseason)+' '+str(np.min(season['expMJD']))+' '+str(np.max(season['expMJD']))
        outfile.write(theline+'\n')

    outfile.close()
    

if Check_m5:
    ratio=[]
    filters=[]

    filt_dict={'u':0,'g':1,'r':2,'i':3,'z':4,'y':5}
    for obs in thedict[0]['dataSlice']:
        filtre=obs['filter'][0]
        Filter_Wavelength_Correction = np.power(500.0 / filterWave[filtre], 0.3)
        Airmass_Correction = math.pow(obs['airmass'],0.6)
        FWHM_Sys = FWHM_Sys_Zenith * Airmass_Correction
        FWHM_Atm = obs['rawSeeing']* Filter_Wavelength_Correction * Airmass_Correction
        finSeeing = scaleToNeff * math.sqrt(np.power(FWHM_Sys,2) + atmNeffFactor * np.power(FWHM_Atm,2))
        print 'hello',finSeeing,obs['rawSeeing'],obs['airmass']
        Tscale = obs['visitExpTime']/ 30.0 * np.power(10.0, -0.4*(obs['filtSkyBrightness'] - msky[filtre]))
        dCm = dCm_infinity[filtre] - 1.25*np.log10(1 + np.power(10.,0.8*dCm_infinity[filtre]- 1.)/Tscale)

        m5_recalc=dCm+Cm[filtre]+0.5*(obs['filtSkyBrightness']-21.)+2.5*np.log10(0.7/finSeeing)-kAtm[filtre]*(obs['airmass']-1.)+1.25*np.log10(obs['visitExpTime']/30.)
        print m5_recalc,obs['fiveSigmaDepth']
        ratio.append(obs['fiveSigmaDepth']/m5_recalc)
        filters.append(filt_dict[filtre])

    figc, axc = plt.subplots(ncols=2, nrows=2, figsize=(10,9))
    axc[0][0].hist(ratio,bins=40,histtype='stepfilled')
    axc[0][1].plot(filters,ratio,'k.')

    figd, axd = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    axd.plot(thedict[0]['dataSlice']['expMJD'],thedict[0]['dataSlice']['airmass'],'b.')

    seasons=Get_Seasons(thedict[0]['dataSlice'])

    outfile= open('Seasons/Seasons_'+opts.fieldtype+'_'+str(opts.fieldnum)+'.txt', 'w')

    for iseason,season in seasons.items():
        axd.plot(season['expMJD'],season['airmass'],'r.') 
        print 'min max',iseason,np.min(season['expMJD']),np.max(season['expMJD'])
        theline=str(iseason)+' '+str(np.min(season['expMJD']))+' '+str(np.max(season['expMJD']))
        outfile.write(theline+'\n')

    outfile.close()



    TSeason_min=-1
    TSeason_max=-1
    
    tseason=11
    if tseason > -1:
        sfile=open('Seasons/Seasons_'+opts.fieldtype+'_'+str(opts.fieldnum)+'.txt', 'r')
        for line in sfile.readlines():
            if int(line.split(' ')[0]) == tseason:
                TSeason_min=line.split(' ')[1]
                TSeason_max=line.split(' ')[2]
                break
            
        if TSeason_min==-1 and TSeason_max==-1:
            print 'Big problem : season not found - Stop'
            

    print 'hello',TSeason_min,TSeason_max


if Make_Rolling:

    seasons=Get_Seasons(thedict[0]['dataSlice'])
    rolling_cadence={}

    for i in range(len(seasons)):
        rolling_cadence[i]=np.zeros((0,1),dtype=seasons[0].dtype)
        

    season_keep,season_remain=Copy_Season(seasons[0])
    rolling_cadence[0]=np.concatenate((rolling_cadence[0],season_keep))

    season_keep,season_remain=Copy_Season(seasons[1])
    rolling_cadence[1]=np.concatenate((rolling_cadence[1],season_keep))

    season_keep,season_remain=Copy_Season(seasons[2],365.,0.8)
    rolling_cadence[1]=np.concatenate((rolling_cadence[1],season_keep))
    rolling_cadence[2]=np.concatenate((rolling_cadence[2],season_remain))

    season_keep,season_remain=Copy_Season(seasons[3],730.,0.8)
    rolling_cadence[1]=np.concatenate((rolling_cadence[1],season_keep))
    rolling_cadence[3]=np.concatenate((rolling_cadence[3],season_remain))

    season_keep,season_remain=Copy_Season(seasons[4])
    rolling_cadence[4]=np.concatenate((rolling_cadence[4],season_keep))
        
    season_keep,season_remain=Copy_Season(seasons[5],365.,0.8)
    rolling_cadence[4]=np.concatenate((rolling_cadence[4],season_keep))
    rolling_cadence[5]=np.concatenate((rolling_cadence[5],season_remain))

    season_keep,season_remain=Copy_Season(seasons[6],730.,0.8)
    rolling_cadence[4]=np.concatenate((rolling_cadence[4],season_keep))
    rolling_cadence[6]=np.concatenate((rolling_cadence[6],season_remain))
    
    season_keep,season_remain=Copy_Season(seasons[7])
    rolling_cadence[7]=np.concatenate((rolling_cadence[7],season_keep))
        
    season_keep,season_remain=Copy_Season(seasons[8],365.,0.8)
    rolling_cadence[7]=np.concatenate((rolling_cadence[7],season_keep))
    rolling_cadence[8]=np.concatenate((rolling_cadence[8],season_remain))

    season_keep,season_remain=Copy_Season(seasons[9],730.,0.8)
    rolling_cadence[7]=np.concatenate((rolling_cadence[7],season_keep))
    rolling_cadence[9]=np.concatenate((rolling_cadence[9],season_remain))


    #print 'hello',seasons[9]['expMJD'],season_remain['expMJD']
    figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    fontsize=15.

    axa.scatter(thedict[0]['dataSlice']['expMJD'],thedict[0]['dataSlice']['airmass'],facecolors='none', edgecolors='r',marker='*')

    for key,val in seasons.items():
        #print key,np.min(val['expMJD']),np.max(val['expMJD'])
        tot=''
        for band in ['u','g','r','i','z','y']:
            select=val[np.where(val['filter']==band)]
            tot+=str(len(select))+' '

        print key,len(val),tot

    for i,val in seasons.items():
        #print i,len(val),len(rolling_cadence[i])
        tot=''
        for band in ['u','g','r','i','z','y']:
            select=rolling_cadence[i][np.where(rolling_cadence[i]['filter']==band)]
            tot+=str(len(select))+' '

        if i == 0:
            rolling_cat=rolling_cadence[i]
        else:
            rolling_cat=np.concatenate((rolling_cat,rolling_cadence[i]))


        print i,len(rolling_cadence[i]),tot

        axa.scatter(rolling_cadence[i]['expMJD'],rolling_cadence[i]['airmass'],facecolors='none', edgecolors='b')
        axa.set_ylabel(r'airmass',{'fontsize': fontsize})
        axa.set_xlabel(r'MJD',{'fontsize': fontsize})
    
        #print rolling_cat,len(rolling_cat)

 
        pkl_file = open('Rolling_Cadence_'+opts.fieldtype+'_'+str(opts.fieldnum)+'.pkl','w')
        pkl.dump(rolling_cat, pkl_file)
        pkl_file.close()


        newseasons=Get_Seasons(rolling_cat)
        outfile= open('Seasons/Seasons_Rolling_'+opts.fieldtype+'_'+str(opts.fieldnum)+'.txt', 'w')

        for iseason,season in newseasons.items():
            #axd.plot(season['expMJD'],season['airmass'],'r.') 
            print 'min max',iseason,np.min(season['expMJD']),np.max(season['expMJD'])
            theline=str(iseason)+' '+str(np.min(season['expMJD']))+' '+str(np.max(season['expMJD']))
            outfile.write(theline+'\n')

        outfile.close()
        

    
    """
    seasons=Get_Seasons(thedict[0]['dataSlice'])

    print seasons[0].dtype,len(seasons)

    for key,val in seasons.items():
        #print key,np.min(val['expMJD']),np.max(val['expMJD'])
        tot=''
        for band in ['u','g','r','i','z','y']:
            select=val[np.where(val['filter']==band)]
            tot+=str(len(select))+' '

        print key,len(val),tot


    
    #selecb=thedict[0]['dataSlice'][np.where(thedict[0]['dataSlice']['expMJD']-np.min(thedict[0]['dataSlice']['expMJD'])<200.)]
    #print 'alors',len(selecb)


    
    figb, axb = plt.subplots(ncols=2, nrows=3, figsize=(10,9))

    for j,band in enumerate(['u','g','r','i','z','y']):
        if j<2:
            k=0
        if j>= 2 and j < 4:
            k=1
        if j>=4:
            k=2

    
        select=thedict[0]['dataSlice'][np.where(thedict[0]['dataSlice']['filter']==band)]
        axb[k][j%2].scatter(select['expMJD']-np.min(select['expMJD']),select['airmass'],facecolors='none', edgecolors='r')
        """


if Draw_Seasons:
    
    figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    fontsize=15.

    axa.scatter(thedict[0]['dataSlice']['expMJD'],thedict[0]['dataSlice']['airmass'],facecolors='none', edgecolors='r',marker='*')

    sfile=open('MJD_'+opts.fieldtype+'_'+str(opts.fieldnum)+'.txt', 'w')
    for data in thedict[0]['dataSlice']:
        print >> sfile,data['expMJD'],data['filter']
    sfile.close()

    rolling_cat=pkl.load(open('Rolling_Cadence_'+opts.fieldtype+'_'+str(opts.fieldnum)+'.pkl','rb'))

    pfile=open('MJD_'+opts.fieldtype+'_'+str(opts.fieldnum)+'_Rolling.txt', 'w')
    for data in rolling_cat:
        print >> pfile,data['expMJD'][0],data['filter'][0]
    pfile.close()


    axa.scatter(rolling_cat['expMJD'],rolling_cat['airmass'],facecolors='none', edgecolors='b')

    axa.set_xlabel(r'expMJD',{'fontsize': 20.})
    axa.set_ylabel(r'Airmass',{'fontsize': 20.})

if Ana_Season:

    def Get_Table(sel_season):

        Tmin=np.min(sel_season['expMJD'])
        Tmax=np.max(sel_season['expMJD'])

        table_meas = Table(names=('day','n_u','n_g','n_r','n_i','n_z','n_y'), dtype=('i4','i4','i4','i4','i4','i4','i4'))
        ntot_obs={}
        for time in np.arange(Tmin,Tmax,1.):
            val={}
            nmeas=0
            for j,band in enumerate(['u','g','r','i','z','y']):
                selb=sel_season[np.where(sel_season['filter']==band)]
                ntot_obs[band]=len(selb)
                selb.sort(order='expMJD')
                selc=selb[np.where(np.logical_and(selb['expMJD']>=time,selb['expMJD']<time+1.))]
                val[band]=len(selc)
                nmeas+=len(selc)
            if nmeas > 0:
                table_meas.add_row((time-Tmin,val['u'],val['g'],val['r'],val['i'],val['z'],val['y']))

        return table_meas, ntot_obs

    #outdir=['/sps/lsst/data/dev/pgris/Make_Cadence/Obs_minion_1016_orig','/sps/lsst/data/dev/pgris/Make_Cadence/Obs_minion_1016_noresh','/sps/lsst/data/dev/pgris/Make_Cadence/Obs_minion_1016']
    outdir=['/sps/lsst/data/dev/pgris/Make_Cadence/Obs_minion_1016_orig','/sps/lsst/data/dev/pgris/Make_Cadence/Obs_minion_1016']

    prefix='Observations'
    fichname=prefix+'_'+opts.fieldname+'_'+str(opts.fieldid)+'.pkl'
    
    thedict={}
    seasons={}
    table_meas={}
    ntot_obs={}
    bands=['u','g','r','i','z','y']

    num_season=1
    for i,val in enumerate(outdir):
        thedict[i]=pkl.load(open(val+'/'+fichname,'rb'))
        seasons[i]=Get_Seasons(thedict[i]['dataSlice'])

    for key,val in seasons.items():
        sel_season=val[num_season]
        selb=sel_season[np.where(sel_season['filter']=='r')]
        print 'hello',key,len(selb)
        table_meas[key],ntot_obs[key]=Get_Table(sel_season)
            
    for band in bands:
        for key,val in ntot_obs.items():
            for keyb,valb in val.items():
                if keyb == band:
                    print keyb,valb
      
    figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))  
    axa.plot(thedict[0]['dataSlice']['expMJD'],thedict[0]['dataSlice']['airmass'],'bo')
    axa.plot(thedict[1]['dataSlice']['expMJD'],thedict[1]['dataSlice']['airmass'],'r*')
    
    test={}
    for i in range(2):
        test[i]=thedict[i]['dataSlice'][(np.where(np.logical_and(thedict[i]['dataSlice']['expMJD']>=59950,thedict[i]['dataSlice']['expMJD']<60250)))]
    

    print 'alors',len(test[0]),len(test[1]),len(thedict[0]['dataSlice']),len(thedict[1]['dataSlice'])

   


    """
    figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(10,9))
    filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}

    for band in ['u','g','r','i','z','y']:
        axa.plot(table_meas['day'],table_meas['n_'+band],filtercolors[band]+'o',label=band)

    axa.legend(loc='upper right')
    axa.set_xlabel(r'Time',{'fontsize': 20.})
    axa.set_ylabel(r'$N_{Obs}$',{'fontsize': 20.})
    axa.set_ylim(-0.5, 4.5)

    stat={}
    for band in ['u','g','r','i','z','y']:
        stat[band]={}
        for nobs in range(0,5):
            sel=table_meas[np.where(table_meas['n_'+band]==nobs)]
            stat[band][nobs]=len(sel)
            print band,nobs,len(sel),ntot_obs[band]
    """    
   

    sel_season=seasons[0][1]
    figb, axb = plt.subplots(ncols=2, nrows=3, figsize=(10,9))

    for j,band in enumerate(['u','g','r','i','z','y']):
        if j<2:
            k=0
        if j>= 2 and j < 4:
            k=1
        if j>=4:
            k=2

        selb=sel_season[np.where(sel_season['filter']==band)]
        selb.sort(order='expMJD')
        diffs=[jo-io for io, jo in zip(selb['expMJD'][:-1], selb['expMJD'][1:])]
        #axb[k][j%2].plot(selb['expMJD'],selb['airmass'],'bo')
        print 'hello',k,j%2,band,diffs
        axb[k][j%2].hist(diffs,bins=80)
 
    plt.show()

if Ana_Cadence:

    seasons=Get_Seasons(thedict[0]['dataSlice'])
    

    dtypes=[('season_id',np.int),('night_id',np.int)]
    
    for band in bands:
        dtypes.append(('n_'+band,np.int))
        dtypes.append(('mean_diff_time_'+band,np.float))
        dtypes.append(('rms_diff_time_'+band,np.float))
        
    tab_resu=np.zeros((60,1),dtype=[dtype for dtype in dtypes])  
    
    ievt=-1
    nnight_tot={}

    for key,season in seasons.items():
        #print key,season
        min_time=np.min(season['expMJD'])
        max_time=np.max(season['expMJD'])
        tmin=int(min_time)
        nnight_tot[key]=int(max_time-min_time)
        while tmin <=max_time:
            tmax=tmin+1.
            select_time=season[np.where(np.logical_and(season['expMJD']>=tmin,season['expMJD']<tmax))]
            print 'eee',len(select_time),tmin,tmax

            if len(select_time) > 0:
                ievt+=1
            
                if len(tab_resu) <= ievt:
                    tab_resu=np.resize(tab_resu,(len(tab_resu)+100,1))

                tab_resu['season_id'][ievt]=key
                tab_resu['night_id'][ievt]=tmin-min_time+1

                for j,band in enumerate(bands):
                    select=select_time[np.where(select_time['filter']==band)]
                    #print key,band,len(select)
                    diff_time=[]

                    if len(select) > 0.:
                        for i in range(len(select)-1):
                            diff_time.append(select['expMJD'][i+1]-select['expMJD'][i])
                    
               
                        tab_resu['n_'+band][ievt]=len(select)
                        tab_resu['mean_diff_time_'+band][ievt]=24.*3600.*np.mean(diff_time)
                        tab_resu['rms_diff_time_'+band][ievt]=24.*3600.*np.std(diff_time)

            tmin+=1

    


        #break

    tab_resu=np.resize(tab_resu,(ievt+1,1))    

    #print tab_resu

    moy_visit={}
    rms_visit={}
    nobs={}
    frequence_obs={}
    frequence_abs={}
    filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}
    fontsize=15

    for iseason in range(10):
    
        sel=tab_resu[np.where(tab_resu['season_id']==iseason)]
    
        moy_visit[iseason]=[]
        rms_visit[iseason]=[]
        nobs[iseason]=[]
        frequence_obs[iseason]=[]
        frequence_abs[iseason]=[]
        filt_num=[]

        Plot_this=False

        if iseason==6:
            Plot_this=True


        if Plot_this:
            figb, axb = plt.subplots(ncols=1, nrows=3, figsize=(10,9))
            figb.suptitle('Field type:'+opts.fieldname+' Field ID: '+str(opts.fieldid)+' - Season '+str(iseason+1))
            

        if iseason == 0:
            print sel

        for j,band in enumerate(bands):
            
            if Plot_this:
                axb[0].plot(sel['night_id'],sel['n_'+band],filtercolors[band]+'o')
                axb[1].errorbar(sel['night_id'],sel['mean_diff_time_'+band],yerr=sel['rms_diff_time_'+band],fmt='None',color = filtercolors[band])
            sell=sel[np.where(sel['n_'+band]>0.)]
            moy_visit[iseason].append(np.mean(sell['n_'+band]))
            rms_visit[iseason].append(np.std(sell['n_'+band]))
            nobs[iseason].append(len(sell['n_'+band]))
            frequence_obs[iseason].append(float(len(sel)/float(len(sell['n_'+band]))))
            frequence_abs[iseason].append(float(nnight_tot[iseason])/float(len(sell['n_'+band])))
            filt_num.append(j)

        if Plot_this: 
            axb[0].set_ylabel(r'N$_{visit}$',{'fontsize': fontsize})
            axb[0].set_xlabel(r'Night number',{'fontsize': fontsize})
            axb[1].set_ylabel(r'<$\Delta T$ (two visits)> (s)',{'fontsize': fontsize})
            axb[1].set_xlabel(r'Night number',{'fontsize': fontsize})
            axb[1].set_ylim(20., 40.)

            axb[2].errorbar(filt_num,moy_visit[iseason],yerr=rms_visit[iseason],fmt='o',ms=5,color='k')
            axb[2].set_xlim(-0.5, 5.5)
            axb[2].set_ylim(5, 30)
            axb[2].set_ylabel(r'<N$_{visit}$>',{'fontsize': fontsize})
            axb[2].set_xlabel(r'Filter',{'fontsize': fontsize})


    figc, axc = plt.subplots(ncols=1, nrows=2, figsize=(10,9))

    myls=['-', '--', '-.', ':','-', '--', '-.', ':','-', '--', '-.', ':']
    mycols=['b','g','r','c','m','y','k','b','g','r','c','m','y','k']

    tot_label=[]
    tit_label=[]
    for iseason in range(10):
        ll='Season '+str(iseason+1)
        tot_label.append(axc[0].errorbar(filt_num,moy_visit[iseason],yerr=rms_visit[iseason],fmt='--o',ms=5,color=mycols[iseason],label=ll))
        tit_label.append(axc[1].errorbar(filt_num,nobs[iseason],linestyle=myls[iseason],color=mycols[iseason],label=ll))
        
    

    axc[0].set_xlim(-0.5, 5.5) 
    axc[0].set_ylabel(r'<N$_{visit}$> (per night)',{'fontsize': fontsize})
    axc[0].set_xlabel(r'Filter',{'fontsize': fontsize})
    labs = [l.get_label() for l in tot_label]
    axc[0].legend(tot_label, labs, ncol=2,loc='upper left',prop={'size':12},frameon=False)
    axc[1].set_xlim(-0.5, 5.5) 
    axc[1].set_ylim(0,70)
    axc[1].set_ylabel(r'$N_{nights}$ ($N^{obs} > 0$)',{'fontsize': fontsize})
    axc[1].set_xlabel(r'Filter',{'fontsize': fontsize})

    
    labs = [l.get_label() for l in tit_label]
    axc[1].legend(tit_label, labs, ncol=2,loc='upper left',prop={'size':12},frameon=False)
    
    figd, axd = plt.subplots(ncols=1, nrows=2, figsize=(10,9))
    tot_label=[]
    tit_label=[]
    
    for iseason in range(10):
        ll='Season '+str(iseason+1)
        #tot_label.append(axd[0].errorbar(filt_num,frequence_obs[iseason],fmt='--o',ms=5,color=mycols[iseason],label=ll))
        tot_label.append(axd[0].errorbar(filt_num,frequence_obs[iseason],linestyle=myls[iseason],color=mycols[iseason],label=ll))
        tit_label.append(axd[1].errorbar(filt_num,frequence_abs[iseason],linestyle=myls[iseason],color=mycols[iseason],label=ll))
        
    axd[0].set_xlim(-0.5,5.5) 
    axd[0].set_ylabel(r'1/<$N_{nights}^{obs}$> (obs nights only)',{'fontsize': fontsize})
    axd[0].set_xlabel(r'Filter',{'fontsize': fontsize})
    labs = [l.get_label() for l in tot_label]
    axd[0].legend(tot_label, labs, ncol=2,loc='upper right',prop={'size':12},frameon=False)
    axd[1].set_xlim(-0.5,5.5) 
    axd[1].set_ylim(1.,40.)
    axd[1].set_ylabel(r'1/<$N_{nights}^{obs}$>',{'fontsize': fontsize})
    axd[1].set_xlabel(r'Filter',{'fontsize': fontsize})

    
    labs = [l.get_label() for l in tit_label]
    axd[1].legend(tit_label, labs, ncol=2,loc='upper right',prop={'size':12},frameon=False)

    


    """

    figb, axb = plt.subplots(ncols=2, nrows=3, figsize=(10,9))

    for j,band in enumerate(['u','g','r','i','z','y']):
        if j<2:
            k=0
        if j>= 2 and j < 4:
            k=1
        if j>=4:
            k=2

    
        select=seasons[0][np.where(seasons[0]['filter']==band)]
        axb[k][j%2].scatter(select['expMJD']-np.min(seasons[0]['expMJD']),select['airmass'],facecolors='none', edgecolors='r')
    
    """
        


if Dump_in_File:
#band:
#mjd:
#exptime:
#realization:
#seeing:
#moon_frac:
#sky:
#kAtm:
#airmass:
#m5sigmadepth_mean:
#m5sigmadepth_recalc:
#Nexp:
    legend='#band #mjd #exptime #seeing #moon_frac #sky #kAtm #airmass #m5sigmadepth #Nexp'
    
    todisplay=['filter','expMJD','expTime','rawSeeing','moonPhase','filtSkyBrightness','airmass']

    toprocess=thedict[0]['dataSlice']
    params=parameters()
    print legend


    outputfile  = open(prefix+'_'+opts.fieldname+'_'+str(opts.fieldid)+addit+'.txt','wb') 
    outputfile.write(legend+'\n')
    for band in bands:
        coadd=Get_coadd(toprocess[np.where(toprocess['filter']==band)])
        #print 'ee',band,len(coadd)
        for key,val in coadd.items():
            filtre=val['filter'][0][0]
            toprint = 'LSST::'+filtre+' '
            toprint+=str(format(np.mean(val['expMJD']),'.7f'))+' '
            toprint+=str(int(np.sum(val['visitExpTime'])))+' '
            toprint+=str(format(np.mean(Get_mean_finSeeing(val)),'.7f'))+' '
            toprint+=str(format(np.mean(val['moonPhase']),'.7f'))+' '
            toprint+=str(format(np.mean(val['filtSkyBrightness']),'.7f'))+' '
            toprint+=str(params.kAtm[filtre])+' '
            toprint+=str(format(np.mean(val['airmass']),'.7f'))+' '
            toprint+=str(format(np.mean(val['fiveSigmaDepth']),'.7f'))+' '
            toprint+=str(len(val))
            #m5_coadd=Get_fiveSigmaDepth_coadd(len(val),np.mean(val['fiveSigmaDepth']))
            m5_coadd=Get_fiveSigmaDepth_coadd(1.,np.mean(val['fiveSigmaDepth']))
            print 'hello',np.mean(val['fiveSigmaDepth']),m5_coadd
            outputfile.write(toprint+'\n')
            print toprint
        
    outputfile.close()




plt.show()
