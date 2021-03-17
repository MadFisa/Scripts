#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals, print_function
__doc__ = """
Script that does default visualizations (marginal plots, 1-d and 2-d).

Author: Johannes Buchner (C) 2013-2019
"""
import numpy
from numpy import exp, log
import matplotlib.pyplot as plt
import sys, os
import json
import pymultinest
import corner

grb=input('enter grb_name:')

import matplotlib.font_manager
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 12})


if len(sys.argv) != 2:
	sys.stderr.write("""SYNOPSIS: %s <output-root>

	output-root: 	Where the output of a MultiNest run has been written to.
	            	Example: chains/1-
%s""" % (sys.argv[0], __doc__))
	sys.exit(1)

prefix = sys.argv[1]
print('model "%s"' % prefix)
if not os.path.exists(prefix + 'params.json'):
	sys.stderr.write("""Expected the file %sparams.json with the parameter names.
For example, for a three-dimensional problem:

["Redshift $z$", "my parameter 2", "A"]
%s""" % (sys.argv[1], __doc__))
	sys.exit(2)
parameters = json.load(open(prefix + 'params.json'))
n_params = len(parameters)

a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = prefix)
s = a.get_stats()

json.dump(s, open(prefix + 'stats.json', 'w'), indent=4)

print('  marginal likelihood:')
print('    ln Z = %.1f +- %.1f' % (s['global evidence'], s['global evidence error']))
print('  parameters:')

for p, m in zip(parameters, s['marginals']):
	lo, hi = m['1sigma']
	med = m['median']
	sigma = (hi - lo) / 2
	if sigma == 0:
		i = 3
	else:
		i = max(0, int(-numpy.floor(numpy.log10(sigma))) + 1)
	fmt = '%%.%df' % i
	fmts = '\t'.join(['    %-15s' + fmt + " +- " + fmt])
	print(fmts % (p, med, sigma))

print('creating marginal plot ...')
data = a.get_data()[:,2:]
weights = a.get_data()[:,0]



import numpy as np
percent_burn=85 # 50% burn-in
Ndata=len(data) # No of elements in data array
# Bitwise check if the index is greater than the index that corresponds to percent burn-in
indices_mask= (np.arange(Ndata)).astype(int) > ((percent_burn/100)*np.ones(Ndata)*Ndata).astype(int)
mask = weights > 1e-5

# Bitwise check if both weight mask and indices_mask are true.
corner.corner(data[mask & indices_mask,:], weights=weights[mask & indices_mask],
labels=parameters, show_titles=True, color='darkblue',title_fmt='.4f')
plt.savefig('GRB_'+grb + '_corner.pdf')
plt.savefig('GRB_'+grb + '_corner.png')
plt.show()
'''
#mask = weights.cumsum() > 1e-5
mask = weights > 1e-4

corner.corner(data[mask,:], weights=weights[mask],
	labels=parameters, show_titles=True, color='k',title_fmt='.3f')
plt.savefig(prefix + 'corner.pdf')
plt.savefig(prefix + 'corner.png')
plt.close()
'''
    
    
    
    
    
    
    
    
    
#-----------------------------------------------------------------------------------------------------    
#-----------------------------------------------------------------------------------------------------
from model_ism import lgfmtd,lgnumtd,lgnuctd,lgnuatd

import par  
import numpy as np
from scipy.integrate import quad

# Constants in cgs
mp=1.6726231*(10**-24) #g
c=2.99792458*(10**10) #cm/s
e= 4.8032*(10**-10) # esu
me=9.1095*(10**-28) # g
sigma_t=6.6524*(10**-25) # cm**2
c1=299792.458 #speed of light in km/s 
H0=70.0     #Hubble Constant in km/s.Mpc

omega_m=0.25
omega_lam=0.7
 
def Integrand(z): 

	return ((1+z)**2*(1+omega_m*z)-z*(2+z)*omega_lam)**(-0.5)

def DL(z): 
    d_l=(1+z)*(c1/H0)*quad(Integrand, 0, z)[0]  #The luminosity distance in Mpc
    return d_l*3.086e+24               #The luminosity distance in cm



# Parameters of GRB
z= par.z #redshift
DL=DL(z) #cm # Luminosity Distance
eE=par.eE 
phi=par.phi  # fluence
k_bol=par.k_bol   # for bat:5 and for fermi:1
e_iso=(k_bol*4*(np.pi)*np.power(DL,2)*phi)/(1+z)
print(e_iso)




parm=[]  
err = []
for x, m in zip(parameters, s['marginals']):
    lo, hi = m['1sigma']
    sigma = (hi - lo) / 2
    err.append(sigma)
    parm.append(m['median'])
    
    
     
thetaj, n, p, eB, eE, en = parm[0],parm[1],parm[2],parm[3],parm[4],parm[5]
#print(thetaj, n, p, eB)    
tobs=np.logspace(1.,8,100)
tj=(1+z)*np.power( ((3*(10**en)*(thetaj)**8)/((2**5)*np.pi*(10**n)*mp*(c**5))), (1./3.))
lgfmtd_= np.ones(len(tobs))
lgnumtd_=np.ones(len(tobs))
lgnuctd_=np.ones(len(tobs))
lgnuatd_=np.ones(len(tobs))
for i in range(0,len(tobs)):
    lgfmtd_=lgfmtd(10**en,10**eB,10**n,tobs[i],tj)
    lgnumtd_[i]=lgnumtd(10**en,10**eB,p,10**eE,tobs[i],tj)
    lgnuctd_[i]=lgnuctd(10**en,10**eB,10**n,tobs[i],tj)
    lgnuatd_[i]=lgnuatd(tobs[i],thetaj,10**en, 10**n, p, 10**eB, 10**eE,tj,10**lgnumtd_[i],10**lgnuctd_[i])    
f=open('freq_evol_'+grb+'.txt','w')
f.write('#time(s)\t  nu_m(Hz)\t  nu_c(Hz)\t nu_a(Hz)\n')
for i in range(len(tobs)):
       f.write('%1.8e \t%1.8e \t%1.8e \t%1.8e \n' %(tobs[i], 10**lgnumtd_[i] ,10**lgnuctd_[i] , 10**lgnuatd_[i]))
f.close()    





eta = e_iso/(e_iso+10**en)
eta_1 = e_iso/(e_iso+10**((en+err[5])))
eta_2 = e_iso/(e_iso+10**((en-err[5])))
e_err1 = eta - eta_1
e_err2 = eta_2 - eta
err = (e_err1+e_err2)/2
print('eta: ',eta,err)
#------------------------------------------------------------------------------------------------------------    
#------------------------------------------------------------------------------------------------------------









