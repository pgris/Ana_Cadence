from math import pow
import numpy as np
import matplotlib.pyplot as plt

mean=-0.062
sig_m=0.032
sig_p=0.032

xmin=mean-5.*sig_m
xmax=mean+5.*sig_p

pas=1.e-4

nsteps=int((xmax-xmin)/pas)

xvals=[]
weights=[]

for i in range(nsteps):
    x=xmin+float(i)*pas
    if x < mean:
        res=np.exp(-np.power(x-mean,2)/(2*np.power(sig_m,2.)))
    else:
        res=np.exp(-np.power(x-mean,2)/(2*np.power(sig_p,2.)))

    xvals.append(x)
    weights.append(res)


print xvals,xmin,xmax,nsteps
weights=weights/np.sum(weights)



plt.plot(xvals,weights,'b.')

"""
xtest=[]

for i in range(2000):
    xtest.append(np.random.choice(xvals,1,p=weights)[0])

plt.hist(xtest,bins=200,normed=1,histtype='stepfilled')
"""

plt.show()
