# This code fits multiwavelength light curve data of GRBs using pymultinest
# Synchrotron Self Absorption is included. Reverse shock is currently not.
import numpy as np
from scipy.integrate import quad
from scipy import optimize
import par
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
'''
eE=par.eE 
phi=par.phi  # fluence
k_bol=par.k_bol   # for bat:5 and for fermi:1
eta=par.eta  # prompt radiation efficiency
eta_1=(1-eta)/eta
E1=(k_bol*4*(np.pi)*np.power(DL,2)*phi)/(1+z)
E=E1*eta_1
en=np.log10(E)
'''
# Smoothing parameters
sc=2.
sm=2.
s=2.

pi=np.pi
# tepoch=0.01*86400 #day

#Dimensionless paramater for p>2
def g(p):
    return (p-2.)/(p-1)

#Peak flux density (mJy) constant density CBM at jet break
def fm0(en,eB,n):
    return 1.6*(1.+z)*np.sqrt(eB/1.e-2)*(en/1.e52)*np.sqrt(n)*np.power(DL/1.e28,-2)
#Above, I changed 10**-2 to 1.e-2, 10**52 to 1.e52, and 10**28 to 1.e28. Is there any conflict due to this in a different machine? RESMI 07APR20
# NO.

#Return log value of flux density (mJy) at some time ts
def lgfmtd(en,eB,n,ts,tj):
    fm_0=fm0(en,eB,n)
    fmtd=fm_0*( 1 + (ts/tj)**s )**(-1./s)
    return np.log10(fmtd)

# nu_m (Hz) for constant density CBM at jet break
def num0(en,eB,p,eE,tj):
    return 3.3*(1.e14)*np.sqrt(1+z)*np.sqrt(eB/1.e-2)*((eE*g(p))**2)*np.sqrt(en/1.e52)*np.power(tj/86400.,-3./2.)
#Did similar 10**x to 1.ex here too. RESMI 07APR20

# Return log value of nu_m (Hz) at some time ts
def lgnumtd(en,eB,p,eE,ts,tj):
    num_0=num0(en,eB,p,eE,tj)
    numtd= num_0*( (ts/tj)**(1.5*s) + (ts/tj)**(2.*s) )**(-1./s)
    return np.log10(numtd)

# nu_c (Hz) for constant density CBM at jet break
def nuc0(en,eB,n,tj):
    test=6.3*(1.e15)*np.power(1+z,-1./2.)*np.power(eB/1.e-2,-3./2.)*np.power(en/1.e52,-1./2.)*(n**-1)*np.power(tj/86400,-1./2.)
    return test

# Return log value of nu_c (Hz) at some time ts
def lgnuctd(en,eB,n,ts,tj):
    nuc_0=nuc0(en,eB,n,tj)
    nuctd= nuc_0*(1+ (ts/tj)**(-0.5*s) )**(1./s)
    return np.log10(nuctd)

# Function 1 required to calculate nu_a
def lgtau_nu1(lgnu,x,tau_num,numtd): # x is the slope
    return np.log10(tau_num)+(x)*(lgnu-np.log10(numtd))

# Function 2 required to calculate nu_a
def lgtau_nu2(lgnu,x,tau_nuc,nuctd): # x is the slope
    return np.log10(tau_nuc)+(x)*(lgnu-np.log10(nuctd))

# Return log value of nu_c (Hz) at some time ts
def lgnuatd(ts,thetaj,en, n, p, eB, eE,tj,numtd,nuctd):
    Gamma_j=1/thetaj
    Rj = 2*(tj/(1+z))*(Gamma_j**2)*c
    if ts<tj:
        Gamma=np.power((3.*en)/((2**5)*pi*((ts/(1+z))**3)*n*mp*(c**5)),1./8.)
        R=2*(ts/(1+z))*(Gamma**2)*c
    else:
        Gamma=Gamma_j*((ts/tj)**(-1./2.))
        R=Rj

    if Gamma < 1.000001:
        lg_nuatd=-10
        return lg_nuatd
    k=0.
    u=-5./3.
    w=-(p+5.)/2.
    if numtd < nuctd:
         v=-(p+4.)/2.
         B=np.sqrt(32*pi*mp*(c**2)*eB*n*(Gamma-1)*Gamma)
         gamma_m= np.sqrt((numtd*(1+z))/(B*Gamma*(e/(2*pi*me*c))))
         tau_num= (5./(3.-k))*(e*n*R)/(B* gamma_m**5)
         if tau_num<1.:
              lg_nuatd=optimize.newton(lgtau_nu1, 15,args=(u,tau_num,numtd))
         else:
              tau_nuc=tau_num*(nuctd/numtd)**(v)
              if tau_nuc<1.:
                   lg_nuatd=optimize.newton(lgtau_nu1, 15, args=(v,tau_num,numtd))
              else:
                   lg_nuatd=optimize.newton(lgtau_nu2, 15, args=(w,tau_nuc,nuctd))

    else:
          v=-3
          B=np.sqrt(32*pi*mp*(c**2)*eB*n*(Gamma-1)*Gamma) # This is from Zhang not P
          gamma_c= np.sqrt((nuctd*(1+z))/(B*Gamma*(e/(2*pi*me*c))))
          tau_nuc= (5./(3.-k))*(e*n*R)/(B* gamma_c**5)
          if tau_nuc<1.:
              lg_nuatd=optimize.newton(lgtau_nu2, 15, args=(u,tau_nuc,nuctd))
          else:
              tau_num=tau_nuc*(numtd/nuctd)**(v)
              if tau_num<1.:
                  lg_nuatd=optimize.newton(lgtau_nu2, 15,args=(v,tau_nuc,nuctd))
              else:
                  lg_nuatd=optimize.newton(lgtau_nu1, 15,args=(w,tau_num,numtd))

    return lg_nuatd

# Returns an array of F_nu (mJy) for given arrays of tobs and nuobs and a given parmlist
def lcflux(tobs,nuobs, parmlist):
    thetaj, n, p, eB,eE, en = parmlist
    tj=(1+z)*np.power( ((3*(10**en)*(thetaj)**8)/((2**5)*np.pi*(10**n)*mp*(c**5))), (1./3.))
    lgfl=np.ones(len(tobs))
    for i in range(0,len(tobs)):
        lgnu=np.log10(nuobs[i])
        lgfmtd_=lgfmtd(10**en,10**eB,10**n,tobs[i],tj)
        lgnumtd_=lgnumtd(10**en,10**eB,p,10**eE,tobs[i],tj)
        lgnuctd_=lgnuctd(10**en,10**eB,10**n,tobs[i],tj)
        lgnuatd_=lgnuatd(tobs[i],thetaj,10**en, 10**n, p, 10**eB, 10**eE,tj,10**lgnumtd_,10**lgnuctd_)

        if lgnuatd_ == -10:
            return -999.

        if lgnumtd_<= lgnuctd_:  #****slow cooling****

               if lgnuatd_ <= lgnumtd_:
                   tp = 10**(-5.*(lgnumtd_-lgnuatd_)/3.0)
               elif lgnumtd_ < lgnuatd_<= lgnuctd_:
                   tp = 10**(-(p+4.)*0.5*(lgnumtd_-lgnuatd_))
               elif lgnuctd_ < lgnuatd_:
                   tp = 10**( (p+4.)*0.5*(lgnuctd_-lgnumtd_)+(p+5.)*0.5*(lgnuatd_-lgnuctd_) )


               taunu = tp*10**((lgnu-lgnumtd_)*(-5./3.))*( (1.+10.**(sm*((p/2.)+(1./3.))*(lgnu-lgnumtd_)))**(-1./sm) )*( (1.+10.**(sc*0.5*(lgnu-lgnuctd_)))**(-1./sc) )

       #================= calculating snu =================

               if lgnumtd_ >= lgnuatd_:
                    smt = lgfmtd_+(5.0*(lgnumtd_-lgnuatd_)/3.0)  #S_nu(nu = nu_p)
               elif lgnumtd_ < lgnuatd_ <= lgnuctd_:
                    smt = lgfmtd_+(p+4.)*0.5*(lgnumtd_-lgnuatd_)
               elif lgnuctd_ < lgnuatd_:
                    smt = lgfmtd_+(p+4.)*0.5*(lgnumtd_-lgnuctd_)+(p+5.)*0.5*(lgnuctd_-lgnuatd_)

               a1 = 10**(2.0*sm*(lgnu-lgnumtd_))
               a2= 10**(5.0*sm*(lgnu-lgnumtd_)*0.5)
               snu = smt+(np.log10(a1+a2))/sm


               if taunu <= 5e-8:
                     lgfl[i] = snu+np.log10(taunu)
               else:
                     lgfl[i] = snu+np.log10(1.-np.exp(-taunu))


        else:   #****fast cooling****

              if lgnuatd_ <= lgnuctd_:
                    tp = 10**(-5.*(lgnuctd_-lgnuatd_)/3.0)
              elif lgnuctd_ < lgnuatd_ <= lgnumtd_:
                    tp = 10**(-(2.+4.)*0.5*(lgnuctd_-lgnuatd_))
              elif lgnumtd_ < lgnuatd_:
                    tp = 10**( 3.*(lgnumtd_-lgnuctd_)+(p+5.)*0.5*(lgnuatd_-lgnumtd_) )


              taunu = tp*10**((lgnu-lgnuctd_)*(-5./3.))*( (1.+10.**(sc*((1.)+(1./3.))*(lgnu-lgnuctd_)))**(-1./sm) )*( (1.+10.**(sm*0.5*(p-1.)*(lgnu-lgnumtd_)))**(-1./sm) )

             #================= calculating snu =================

              if lgnuctd_ >= lgnuatd_:
                    smt = lgfmtd_+(5.0*(lgnuctd_-lgnuatd_)/3.0)  #S_nu(nu = nu_p)
              elif lgnuctd_ < lgnuatd_ <= lgnumtd_:
                    smt = lgfmtd_+(2+4.)*0.5*(lgnuctd_-lgnuatd_)
              elif lgnumtd_ < lgnuatd_:
                    smt = lgfmtd_+3.0*(lgnuctd_-lgnumtd_)+(p+5.)*0.5*(lgnumtd_-lgnuatd_)


              a1 = 10**(2.0*sm*(lgnu-lgnuctd_))
              a2= 10**(5.0*sm*(lgnu-lgnuctd_)*0.5)
              snu = smt+(np.log10(a1+a2))/sm

              #fnu = snu
              if taunu <= 5e-8:
                  lgfl[i] = snu+np.log10(taunu)
              else:
                  lgfl[i] = snu+np.log10(1.-np.exp(-taunu))

    return 10**lgfl
