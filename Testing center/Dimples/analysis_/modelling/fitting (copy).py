# This code fits multiwavelength light curve data of GRBs using pymultinest
# Synchrotron Self Absorption is included. Reverse shock is currently not.
from model_ism import lcflux
import numpy as np
import json
import par
from model_ism import lcflux
from scipy import optimize
import os

# You will need to install pymultinest module for python
# http://johannesbuchner.github.io/PyMultiNest/install.html
from pymultinest.solve import solve
ndim=6 #Number of parameters in the model
datafile='mltnst_dim6' # Name of the 'PyMultiNest' datafile
# This is NOT the obs datafile
# A number of files with this prefix will be created in the directory after
# running pymultinest



#Importing data
grb= input('enter grb_name:')
txt = np.genfromtxt('grbdata_'+grb+'.txt', dtype=float, )
tobs=txt[:,0]
nuobs=txt[:,1]
fdata= txt[:,2]
erdata = txt[:,3]

# grbdata.txt is the observed datafile. It has 4 columns seperated by commas.
# tobs (s), nuobs (Hz), F_nu (mJy), error (mJy)

theta_up = par.theta_up
theta_lo = par.theta_lo
n0_up = par.n0_up
n0_lo = par.n0_lo
p_up = par.p_up
p_lo = par.p_lo
eB_up = par.eB_up
eB_lo = par.eB_lo
en_up = par.en_up
en_lo = par.en_lo
eE_up = par.eE_up
eE_lo = par.eE_lo
print(theta_up,theta_lo,n0_up,n0_lo,p_up,p_lo,eB_up,eB_lo)

#Multinest specific functions
def prior(cube):
    # Make sure the order of these cube paramaters are the same as the parmlist
    # in lcflux function namely thetaj, en, n, p, eB, eE
    cube[0] = cube[0]*(np.radians(theta_up)-np.radians(theta_lo)) +np.radians(theta_lo)   #thetaJ
    cube[1] = cube[1]*(n0_up - n0_lo) + n0_lo  #n0
    cube[2] = cube[2]*(p_up - p_lo) + p_lo  #pFS
    cube[3] = cube[3]*(eB_up - eB_lo) + eB_lo   #epsB
    cube[4] = cube[4]*(eE_up - eE_lo) + eE_lo 
    cube[5] = cube[5]*(en_up - en_lo) + en_lo  #epsE

    return cube

#multinest function
def loglike(cube):
    chi2=0.0
    parmlist=np.ones(ndim)

    for i in range(ndim):
          parmlist[i] = cube[i]
    thetaj, n, p, eB, eE, en = parmlist

    chi2 = (((lcflux(tobs,nuobs,parmlist) - fdata)/erdata)**2).sum()
    mxlik =  -0.5*chi2-0.5*sum(np.log(2.*np.pi*erdata**2))
    return mxlik

# This code both runs PyMultiNest and generates forest plots
# To run PyMultiNest uncomment the below block while commenting the plotting block

#--------------------------running multinest------------------
nlive=par.nlive # Number of live walkers
tol=0.3

parameters=[r'$\theta_{j}$', r'$log n_{0}$','p',r'$log \epsilon_{B}$',r'$log \epsilon_{E}$',r'$log Ek$']
#parameters=['s',r'$p_{FS}$',r'$\lambda$',r'$\iota$',r'$log\/\nu_m$',r'$log\/\nu_c$',r'$log\/f_m$',r'$log\/\tau_m$']

if not os.path.exists(grb):
    os.makedirs(grb)
solve(loglike, prior, n_dims=ndim, outputfiles_basename=grb+'/'+ datafile + '_', resume = False, verbose = False,n_live_points=nlive,sampling_efficiency=0.3)
json.dump(parameters, open(grb+'/'+datafile + '_params.json', 'w')) # save parameter names

#-----------------------multinest OVER--------------------------



