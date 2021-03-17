import numpy as np
import matplotlib.pyplot as plt
from progThinISM import progthinism as prog
import json
from pymultinest.solve import solve

ndim =7
datafile='mltnst_dim7_Xradmmscint'

#-------------Fixing parameters ------------------
#ppFS,lgeps_e,lgeps_B,lgE52,lgdens
def prior(cube):
    cube[0] = cube[0]*(3.-2.01) +2.01  #pFS
    cube[1] = cube[1]*(1.+3.) -3.  #lgeps_e
    cube[2] = cube[2]*(-0.5+6.) -6.  #lgeps_B
    cube[3] = cube[3]*(0.7-0.05) +0.05  #thetaJ
    cube[4] = cube[4]*(2.+5) -5  #lgE52
    cube[5] = cube[5]*(4.+4) -4.  #lgdens
    cube[6] = cube[6]*(1000-50) +50   #Eta
    return cube


#----------------------Loading X-ray -------------
'''
d1=np.loadtxt('xrt_lc.dat')
nuobs=d1[:,0]
tsec=d1[:,1]*86400.
fdata=d1[:,2]
er1=fdata-d1[:,3]
er2=d1[:,4]-fdata
erdata = np.where(fdata==0.,d1[:,3],(er1+er2)*0.5)




'''
#----------------------Loading millimeter -------------
dta=np.loadtxt('ALMA97ghz_scint.txt')
nuobs = dta[:,0]*1.e9
tsec = 86400.*dta[:,1]
fdata = dta[:,2]
erdata = dta[:,3]+dta[:,4]





'''
#----------------------Loading X-ray and millimeter/radio together-------------
d1=np.loadtxt('xrt_lc.dat')
er1=d1[:,2]-d1[:,3]
er2=d1[:,4]-d1[:,2]
erdata1= np.where( d1[:,2]==0.,d1[:,3],(er1+er2)*0.5)


d2=np.loadtxt('radiomm_scint.txt')

nuobs=np.append(d1[:,0],d2[:,0]*1.e9)
tsec=86400.*np.append(d1[:,1],d2[:,1])
fdata=np.append(d1[:,2],d2[:,2])
erdata=np.append(erdata1,d2[:,3]+d2[:,4])
'''


#------------------------Initializing arrays----------------
Fnutot=np.ones(len(tsec))
FnuFS=np.ones(len(tsec))
FnuRS=np.ones(len(tsec))


def loglike(cube):
    chi2=0.0
#some default parameters. will get overwritten by cube depending on number of free parameters
    Eta = 500; thetaJ=0.25
    ppFS = 2.5
    lgeps_e =-1.0
    lgeps_B=-2.5
    lgE52= 1.2
    lgdens= -1.9
    lgRB = np.log10(1.)
    ppRS = 2.1
    parms=np.ones(ndim)

    for i in range(ndim):
          parms[i] = cube[i]

    ppFS,lgeps_e,lgeps_B,thetaJ,lgE52,lgdens,Eta = parms
    parmsin=np.array([ppFS,lgeps_e,lgeps_B,ppRS,lgRB,thetaJ,lgE52,lgdens,Eta])

    Fnutot=np.ones(len(tsec))
    
    for i in range(len(tsec)):
         FnuFS[i],FnuRS[i] = prog.aglightcurve(np.log10(nuobs[i]),tsec[i],parmsin)
         Fnutot[i]=FnuFS[i]+FnuRS[i]
    chi2= (((Fnutot- fdata)/erdata)**2).sum()
    mxlik =  -0.5*chi2-0.5*sum(np.log(2.*np.pi*erdata**2))
    return mxlik




'''
#--------------------------multinest------------------

nlive=1000
tol=0.3

#ppFS,lgeps_e,lgeps_B,thetaJ,lgE52,lgdens,Eta = parms
parameters=[r'$p_{FS}$',r'$log\/\epsilon_e$',r'$log\/\epsilon_B$',r'$\theta_j$',r'$log\/E_{52}$',r'$log\/n_0$',r'$\eta$']



solve(loglike, prior, n_dims=ndim, outputfiles_basename=datafile + '_', resume = True, verbose = False,n_live_points=nlive,sampling_efficiency=0.3)
json.dump(parameters, open(datafile + '_params.json', 'w')) # save parameter names

#-----------------------multinest OVER--------------------------
'''



#-----Plotting data 
plt.errorbar(tsec/86400.,fdata,yerr=erdata,ms=10,fmt='wo',mec='salmon',mfc='salmon',ecolor='salmon',label='97GHz')


'''
#---------------------------Running one solution


Eta = 500; thetaJ=0.25
ppFS = 2.5
lgeps_e =-1.0
lgeps_B=-2.5
lgE52= 1.2
lgdens= -1.9
lgRB = np.log10(1.)
ppRS = 2.1
parmsin=np.array([ppFS,lgeps_e,lgeps_B,ppRS,lgRB,thetaJ,lgE52,lgdens,Eta])

Fnutot=np.ones(len(tsec))

for i in range(len(tsec)):
     FnuFS[i],FnuRS[i] = prog.aglightcurve(np.log10(nuobs[i]),tsec[i],parmsin)
     Fnutot[i]=FnuFS[i]+FnuRS[i]

plt.plot(tsec/86400.,Fnutot,c='purple',ls='-',lw=2,alpha=0.4)
'''






#------------------------Running forest---------
nuobs0=nuobs[0]
tsec=np.logspace(-1.5,2.5,50)*86400.  #radio/mm
#tsec = np.logspace(-3.3,1.5,50)*86400.   #xray

Fnutot=np.ones(len(tsec))
FnuFS=np.ones(len(tsec))
FnuRS=np.ones(len(tsec))
col=['k','c','purple','forestgreen','y','maroon','darkorchid','crimson','yellowgreen','cadetblue','darkorange']

forest=datafile + "_ev.dat"
sm=np.loadtxt(forest)
Eta = 500; thetaJ=0.25
ppFS = 2.5
lgeps_e =-1.0
lgeps_B=-2.5
lgE52= 1.2
lgdens= -1.9
lgRB = np.log10(1.)
ppRS = 2.1
parms=np.ones(ndim)

sm0 =sm[:,0]
sm1 = sm[:,1]
sm2 =sm[:,2]
sm3 = sm[:,3]
sm4 = sm[:,4]
sm5 = sm[:,5]
sm6 = sm[:,6]

tot=len(sm0)

burn = int(tot*3./5.)  #Taking 60% as burn-in

nlines=100

for i in range(nlines):
   n=int(np.random.uniform(1,tot-burn))
   ppFS=sm0[n+burn]
   lgeps_e = sm1[n+burn]
   lgeps_B=sm2[n+burn]
   thetaJ = sm3[n+burn]
   lgE52 = sm4[n+burn]
   lgdens = sm5[n+burn]
   Eta = sm6[n+burn]

#   ppFS,lgeps_e,lgeps_B,thetaJ,lgE52,lgdens,Eta = parms    
#
   parms=np.array([ppFS,lgeps_e,lgeps_B,ppRS,lgRB,thetaJ,lgE52,lgdens,Eta])

   for i in range(len(tsec)):
         FnuFS[i],FnuRS[i] = prog.aglightcurve(np.log10(nuobs0),tsec[i],parms)
         Fnutot[i]=FnuFS[i]+FnuRS[i]
   plt.plot(tsec/86400.,FnuFS,c='gray',ls='--',lw=1,alpha=0.2)
   plt.plot(tsec/86400.,FnuRS,c='gray',ls='-.',lw=3,alpha=0.2)
   plt.plot(tsec/86400.,Fnutot,c='maroon',ls='-',lw=2,alpha=0.2)






#--------plot settings-----

#plt.legend(ncol=5,loc=0)
#plt.gca().set_xlim(5.e-4,50)  #for xray
#plt.gca().set_ylim(5.e-8,20)
plt.gca().set_xlim(1.e-2,300) #for radio/mm
plt.gca().set_ylim(0.001,100)
plt.gca().set_xscale('log',basex=10,size=30)  # define log axis
plt.gca().set_yscale('log',basex=10,size=30)
plt.xlabel('Time since burst (days)')
plt.ylabel(r'$F_{\nu}$ (mJy)')
plt.legend(loc=0,numpoints = 1)
plt.savefig("ALMA97-ism.pdf")
plt.show()
