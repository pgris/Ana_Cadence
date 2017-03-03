import matplotlib.pyplot as plt
import scipy
import scipy.stats
import pickle as pkl
import numpy as np

pkl_file = open('output.pkl','rb')
thedict=pkl.load(pkl_file)
dist_name='chi2'

fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10,9))

chi2=[]
for vals in thedict['chi2']:
    if vals[1]>10. and vals[1]<12. and vals[0]< 100.:
        chi2.append(vals[0])

ax.hist(chi2, bins= 100,histtype='step',color='k')

"""
dist = getattr(scipy.stats, 'chi2')
param = dist.fit(thedict['chi2'])
print param
size = np.max(thedict['chi2'])
x = np.arange(0.,1.,0.1)
pdf_fitted = dist.pdf(x, *param[:-2], loc=param[-2], scale=param[-1])*120.
plt.plot(pdf_fitted, label=dist_name)
"""
   
xvals = np.arange(0.,10.,0.1)
 
chi2_th={}

for x in xvals:
    for k in range(1,8):
        if not chi2_th.has_key(k):
            chi2_th[k]=[]
        rat=float(k)/2.
        chi2_th[k].append(np.power(0.5,rat)*np.power(x,rat-1.)*np.exp(-x/2)/scipy.special.gamma(rat)/float(k))

figb, axb = plt.subplots(ncols=1, nrows=1, figsize=(10,9))

print len(xvals),len(chi2_th[1])
axb.plot(xvals,chi2_th[1],'k.')
axb.plot(xvals,chi2_th[2],'r.')
axb.plot(xvals,chi2_th[3],'b.')
axb.plot(xvals,chi2_th[7],'g.')

print xvals,chi2_th[1],chi2_th[2],chi2_th[3],np.mean(chi2_th[2]),np.mean(chi2_th[3])

plt.show()
