


#----------------------------------------------------------------------------------------------------------------------------
#                                                  plotting.py
#----------------------------------------------------------------------------------------------------------------------------



import numpy as np
from pylab import *
from astropy.io import ascii
from astropy import units
import matplotlib
import matplotlib.pyplot as plt
from model_ism import lcflux


import matplotlib
font = {'family'         : 'serif',
	'weight'         : 'bold',
	'size'	         : 16}
matplotlib.rc('font',**font)
matplotlib.rc('grid',linewidth=1)
matplotlib.rc('xtick.major',width=2)
matplotlib.rc('xtick.major',size=7)
matplotlib.rc('xtick.minor',width=2)
matplotlib.rc('xtick.minor',size=4)
matplotlib.rc('ytick.major',width=2)
matplotlib.rc('ytick.major',size=7)
matplotlib.rc('ytick.minor',width=2)
matplotlib.rc('ytick.minor',size=4)
plt.style.use('classic')
fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(111)
ax.tick_params(direction="in")
#ax.yaxis.set_label_position("right")
#ax.yaxis.tick_right()
ax.tick_params(axis='y',which='minor',length=6,width=2,labelsize=18)
ax.tick_params(axis='y',which='major',length=10,width=2,labelsize=18)
ax.tick_params(axis='x',which='minor',length=6,width=2,labelsize=18)
ax.tick_params(axis='x',which='major',length=10,width=2,labelsize=18)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r"Time since burst $(\rm s)$",fontsize=20)
ax.set_ylabel(r'Flux Density $(\rm mJy)$',fontsize=20)


'''
ax2=ax.twinx()

ax2.yaxis.label.set_color('blue')
ax2.tick_params(axis='y',which='major', colors='blue',length=10,width=2,)
ax2.tick_params(axis='y',which='minor', colors='blue',length=6,width=2,)
ax2.spines['right'].set_color('blue')

ax2.tick_params(direction="in")
ax2.tick_params(axis='y',which='minor',length=6,width=2,labelsize=18)
ax2.tick_params(axis='y',which='major',length=10,width=2,labelsize=18)
#ax2.yaxis.set_ticks_position('both')
#ax2.xaxis.set_ticks_position('both')

ax2.set_xscale('log')
ax2.set_yscale('log')

ax2.set_ylabel(r'Frequency $(\rm Hz)$',fontsize=20)

ax2.set_xscale('log')
ax2.set_yscale('log')

'''

grb = input('enter grb name: ')

t1, tpos1, tneg1 ,  E1,flux1 ,fluxpos1 , fluxneg1 = np.loadtxt( 'fluxdensity@5keV_xrt_'+grb+'.txt' ,dtype=float,unpack=True,usecols=[0,1,2,3,4,5,6])
t, w, zp, m, merr,f,ferr = np.loadtxt( 'fluxdensity_op_'+grb+'.txt' ,dtype=float,unpack=True,usecols=[0,1,2,3,4,5,6])
#t0, w0, zp0, m0, merr0,f0,ferr0 = np.loadtxt( 'fluxdensity_op_ul_'+grb+'.txt' ,dtype=float,unpack=True,usecols=[0,1,2,3,4,5,6])
#tobs, nu_m, nu_c, nu_a = np.loadtxt( 'freq_evol_'+grb+'.txt' ,dtype=float,unpack=True,usecols=[0,1,2,3])


#ax2.plot(tobs,nu_m,c='cyan',ls='--',lw=5,label=r'$\nu_{m}$')            
#ax2.plot(tobs,nu_c,c='violet',ls='--',lw=5,label=r'$\nu_{c}$')    
#ax2.plot(tobs,nu_a,c='black',ls='--',lw=5,label=r'$\nu_{a}$')  
   


mask = [0.2684,0.2086,0.2246,0.3886,0.3520,0.5411,0.4346,  0.3600,0.4380,0.5450,0.6410,0.7980,1.220,1.630,2.190,  0.3560,0.4830,0.6260, 0.7670, 0.9100]


m0 = w == mask[0] # uvot-uvw1
m1 = w == mask[1] # uvot-uvw2  
m2 = w == mask[2] # uvot-uvm2
m3 = w == mask[3] # uvot-white
m4 = w == mask[4] # uvot-U
m5 = w == mask[5] # uvot-B
m6 = w == mask[6] # uvot-V


m7 = w == mask[7] # U band
m8 = w == mask[8] # B band
m9 = w == mask[9] # V band
m10 = w == mask[10] # R band
m11 = w == mask[11] # I band
m12 = w == mask[12] # J band
m13 = w == mask[13] # H band
m14 = w == mask[14] # K band


m15 = w == mask[15] # u band
m16 = w == mask[16] # g band
m17 = w == mask[17] # r band
m18 = w == mask[18] # i band
m19 = w == mask[19] # z band



t_uvot_uvw1 ,t_uvot_uvw2 ,t_uvot_uvm2 ,t_uvot_white ,t_uvot_u ,t_uvot_b ,t_uvot_v =t[m0],t[m1],t[m2],t[m3],t[m4],t[m5],t[m6]
t_U,t_B,t_V,t_R,t_I,t_J,t_H,t_K = t[m7],t[m8],t[m9],t[m10],t[m11],t[m12],t[m13],t[m14]
t_u,t_g,t_r,t_i,t_z = t[m15],t[m16],t[m17],t[m18],t[m19]

f_uvot_uvw1 ,f_uvot_uvw2 ,f_uvot_uvm2 ,f_uvot_white ,f_uvot_u ,f_uvot_b ,f_uvot_v =f[m0],f[m1],f[m2],f[m3],f[m4],f[m5],f[m6]
f_U,f_B,f_V,f_R,f_I,f_J,f_H,f_K = f[m7],f[m8],f[m9],f[m10],f[m11],f[m12],f[m13],f[m14]
f_u,f_g,f_r,f_i,f_z = f[m15],f[m16],f[m17],f[m18],f[m19]

ferr_uvot_uvw1 ,ferr_uvot_uvw2 ,ferr_uvot_uvm2 ,ferr_uvot_white ,ferr_uvot_u ,ferr_uvot_b ,ferr_uvot_v=ferr[m0],ferr[m1],ferr[m2],ferr[m3],ferr[m4],ferr[m5],ferr[m6]
ferr_U,ferr_B,ferr_V,ferr_R,ferr_I,ferr_J,ferr_H,ferr_K = ferr[m7],ferr[m8],ferr[m9],ferr[m10],ferr[m11],ferr[m12],ferr[m13],ferr[m14]
ferr_u,ferr_g,ferr_r,ferr_i,ferr_z = ferr[m15],ferr[m16],ferr[m17],ferr[m18],ferr[m19]
'''

mask = [0.2684,0.2086,0.2246,0.3886,0.3520,0.5411,0.4346,  0.3600,0.4380,0.5450,0.6410,0.7980,1.220,1.630,2.190,  0.3560,0.4830,0.6260, 0.7670, 0.9100]


m0 = w0 == mask[0] # uvot-uvw1
m1 = w0 == mask[1] # uvot-uvw2  
m2 = w0 == mask[2] # uvot-uvm2
m3 = w0 == mask[3] # uvot-white
m4 = w0 == mask[4] # uvot-U
m5 = w0 == mask[5] # uvot-B
m6 = w0 == mask[6] # uvot-V


m7 = w0 == mask[7] # U band
m8 = w0 == mask[8] # B band
m9 = w0 == mask[9] # V band
m10 = w0 == mask[10] # R band
m11 = w0 == mask[11] # I band
m12 = w0 == mask[12] # J band
m13 = w0 == mask[13] # H band
m14 = w0 == mask[14] # K band


m15 = w0 == mask[15] # u band
m16 = w0 == mask[16] # g band
m17 = w0 == mask[17] # r band
m18 = w0 == mask[18] # i band
m19 = w0 == mask[19] # z band


t_uvot_uvw1_ul,t_uvot_uvw2_ul,t_uvot_uvm2_ul,t_uvot_white_ul,t_uvot_u_ul,t_uvot_b_ul,t_uvot_v_ul=t0[m0],t0[m1],t0[m2],t0[m3],t0[m4],t0[m5],t0[m6]
t_U_ul,t_B_ul,t_V_ul,t_R_ul,t_I_ul,t_J_ul,t_H_ul,t_K_ul= t0[m7],t0[m8],t0[m9],t0[m10],t0[m11],t0[m12],t0[m13],t0[m14]
t_u_ul,t_g_ul,t_r_ul,t_i_ul,t_z_ul = t0[m15],t0[m16],t0[m17],t0[m18],t0[m19]

f_uvot_uvw1_ul,f_uvot_uvw2_ul,f_uvot_uvm2_ul,f_uvot_white_ul,f_uvot_u_ul,f_uvot_b_ul,f_uvot_v_ul=f0[m0],f0[m1],f0[m2],f0[m3],f0[m4],f0[m5],f0[m6]
f_U_ul,f_B_ul,f_V_ul,f_R_ul,f_I_ul,f_J_ul,f_H_ul,f_K_ul=f0[m7],f0[m8],f0[m9],f0[m10],f0[m11],f0[m12],f0[m13],f0[m14]
f_u_ul,f_g_ul,f_r_ul,f_i_ul,f_z_ul = f0[m15],f0[m16],f0[m17],f0[m18],f0[m19]

ferr_uvot_uvw1_ul,ferr_uvot_uvw2_ul,ferr_uvot_uvm2_ul,ferr_uvot_white_ul,ferr_uvot_u_ul,ferr_uvot_b_ul,ferr_uvot_v_ul=ferr0[m0],ferr0[m1],ferr0[m2],ferr0[m3],ferr0[m4],ferr0[m5],ferr0[m6]
ferr_U_ul,ferr_B_ul,ferr_V_ul,ferr_R_ul,ferr_I_ul,ferr_J_ul,ferr_H_ul,ferr_K_ul=ferr0[m7],ferr0[m8],ferr0[m9],ferr0[m10],ferr0[m11],ferr0[m12],ferr0[m13],ferr0[m14]
ferr_u_ul,ferr_g_ul,ferr_r_ul,ferr_i_ul,ferr_z_ul = ferr0[m15],ferr0[m16],ferr0[m17],ferr0[m18],ferr0[m19]

'''


parameters=[r'$\theta_{j}$',r'$log n_{0}$','p',r'$log \epsilon_{B}$']
col=['darkslategray','violet','darkslateblue','olive','turquoise','darkgoldenrod','deepskyblue','yellow','steelblue','darkorange','red','saddlebrown','maroon','lightcoral','hotpink','olivedrab','green','darksalmon','peru','brown','olive','grey']
#x-ray[21],uvot-uvw1[0],uvot-uvw2[1],uvot-uvm2[2],uvot-white[3],uvot-U[4],uvot-V[5],uvot-B[6],U[7],B[8],V[9],R[10],I[11],J[12],H[13],K[14],u[15],g[16],r[17],i[18],z[19],rad[20]
nuplt=np.array([1.11695231e+15,1.43715244e+15,1.33477293e+15,7.71461657e+14,8.51676136e+14,5.54038071e+14,6.89806719e+14,8.32750000e+14,6.84452055e+14,5.50073394e+14,4.67691108e+14,3.75676692e+14,2.45729508e+14,1.83920245e+14,1.36890411e+14,8.42106742e+14,6.20683230e+14,4.78897764e+14,3.90860495e+14,3.29439560e+14,8.46e9,1.20900000e+18])


nu_UVW1,nu_UVW2,nu_UVM2,nu_Wh,nu_U1,nu_V1,nu_B1,nu_U,nu_B,nu_V,nu_R,nu_I,nu_J,nu_H,nu_K,nu_u,nu_g,nu_r,nu_i,nu_z,nu_rad,nu_x=nuplt[0],nuplt[1],nuplt[2],nuplt[3],nuplt[4],nuplt[5],nuplt[6],nuplt[7],nuplt[8],nuplt[9],nuplt[10],nuplt[11],nuplt[12],nuplt[13],nuplt[14],nuplt[15],nuplt[16],nuplt[17],nuplt[18],nuplt[19],nuplt[20],nuplt[21]




#---------------------------------------------------------------------------------------------------------------------------

ax.errorbar(t1,flux1,  yerr = (-fluxneg1,fluxpos1), xerr = (-tneg1,tpos1),marker='D',ls='none',markersize='14',lw=2
                 ,color=col[21], ecolor=col[21],capsize=3,label='XRT@5KeV' )


if len(t_uvot_uvw1)!=0:
    ax.errorbar(t_uvot_uvw1,f_uvot_uvw1,marker='s', yerr = ferr_uvot_uvw1,ls='none',
                 markersize='14',color=col[0], ecolor=col[0],capsize=3,label='UVOT_UVW1' )
if len(t_uvot_uvw2)!=0:
    ax.errorbar(t_uvot_uvw2,f_uvot_uvw2,marker='s', yerr = ferr_uvot_uvw2,ls='none',
                 markersize='14',color=col[1], ecolor=col[1],capsize=3,label='UVOT_UVW2' )
if len(t_uvot_uvm2)!=0:
    ax.errorbar(t_uvot_uvm2,f_uvot_uvm2,marker='s', yerr = ferr_uvot_uvm2,ls='none',
                 markersize='14',color=col[2], ecolor=col[2],capsize=3,label='UVOT_UVM2' )    
if len(t_uvot_white)!=0:
    ax.errorbar(t_uvot_white,f_uvot_white,marker='o', yerr = ferr_uvot_white,ls='none',
                 markersize='14',color=col[3], ecolor=col[3],capsize=3,label=r'UVOT_White' )
if len(t_uvot_u)!=0:
    ax.errorbar(t_uvot_u,f_uvot_u, yerr = ferr_uvot_u,marker='s',ls='none',markersize='14',lw=2
                 ,color=col[4], ecolor=col[4],capsize=3,label=r'UVOT_U' )    
if len(t_uvot_v)!=0:
    ax.errorbar(t_uvot_v,f_uvot_v,marker='s', yerr = ferr_uvot_v,ls='none',
                 markersize='14',color=col[5], ecolor=col[5],capsize=3,label='UVOT_V' )

if len(t_uvot_b)!=0:
    ax.errorbar(t_uvot_b,f_uvot_b,marker='s', yerr = ferr_uvot_b,ls='none',
                 markersize='14',color=col[6], ecolor=col[6],capsize=3,label='UVOT_B' )

if len(t_U)!=0:
    ax.errorbar(t_U,f_U,marker='o', yerr = ferr_U,ls='none',
                 markersize='14',color=col[7], ecolor=col[7],capsize=3,label='U' ) 
if len(t_B)!=0:
    ax.errorbar(t_B,f_B,marker='o', yerr = ferr_B,ls='none',
                 markersize='14',color=col[8], ecolor=col[8],capsize=3,label='B' )    
if len(t_V)!=0:
    ax.errorbar(t_V,f_V,marker='o', yerr = ferr_V,ls='none',
                 markersize='14',color=col[9], ecolor=col[9],capsize=3,label='V' )      
if len(t_R)!=0:
    ax.errorbar(t_R,f_R,marker='o', yerr = ferr_R,ls='none',
                 markersize='14',color=col[10], ecolor=col[10],capsize=3,label='R' )     
if len(t_I)!=0:
    ax.errorbar(t_I,f_I,marker='o', yerr = ferr_I,ls='none',
                 markersize='14',color=col[11], ecolor=col[11],capsize=3,label='I' )
if len(t_J)!=0:
    ax.errorbar(t_J,f_J,marker='o', yerr = ferr_J,ls='none',
                 markersize='14',color=col[12], ecolor=col[12],capsize=3,label='J' )      
if len(t_H)!=0:
    ax.errorbar(t_H,f_H,marker='o', yerr = ferr_H,ls='none',
                 markersize='14',color=col[13], ecolor=col[13],capsize=3,label='H' )     
if len(t_K)!=0:
    ax.errorbar(t_K,f_K,marker='o', yerr = ferr_K,ls='none',
                 markersize='14',color=col[14], ecolor=col[14],capsize=3,label='K' )   
    
if len(t_u)!=0:
    ax.errorbar(t_u,f_u,marker='v', yerr = ferr_u,ls='none',
                 markersize='14',color=col[15], ecolor=col[15],capsize=3,label='u' ) 
if len(t_g)!=0:
    ax.errorbar(t_g,f_g,marker='v', yerr = ferr_g,ls='none',
                 markersize='14',color=col[16], ecolor=col[16],capsize=3,label='g' )    
if len(t_r)!=0:
    ax.errorbar(t_r,f_r,marker='v', yerr = ferr_r,ls='none',
                 markersize='14',color=col[17], ecolor=col[17],capsize=3,label='r' )      
if len(t_i)!=0:
    ax.errorbar(t_i,f_i,marker='v', yerr = ferr_i,ls='none',
                 markersize='14',color=col[18], ecolor=col[18],capsize=3,label='i' )     
if len(t_z)!=0:
    ax.errorbar(t_z,f_z,marker='v', yerr = ferr_z,ls='none',
                 markersize='14',color=col[19], ecolor=col[19],capsize=3,label='z' ) 

''' 
    
uplims = np.array([1], dtype=bool)    
    
    


if len(t_uvot_uvw1_ul)!=0:
    ax.errorbar(t_uvot_uvw1_ul,f_uvot_uvw1_ul,marker='s', yerr = ferr_uvot_uvw1_ul,ls='none',
                 markersize='14',color=col[0], ecolor=col[0],capsize=3,uplims=uplims,)
if len(t_uvot_uvw2_ul)!=0:
    ax.errorbar(t_uvot_uvw2_ul,f_uvot_uvw2_ul,marker='s', yerr = ferr_uvot_uvw2_ul,ls='none',
                 markersize='14',color=col[1], ecolor=col[1],capsize=3,uplims=uplims,label=r'UVOT_UVW2' )
if len(t_uvot_uvm2_ul)!=0:
    ax.errorbar(t_uvot_uvm2_ul,f_uvot_uvm2_ul,marker='s', yerr = ferr_uvot_uvm2_ul,ls='none',
                 markersize='14',color=col[2], ecolor=col[2],capsize=3,uplims=uplims,label=r'UVOT_UVM2' )    
if len(t_uvot_white_ul)!=0:
    ax.errorbar(t_uvot_white_ul,f_uvot_white_ul,marker='o', yerr = ferr_uvot_white_ul,ls='none',
                 markersize='14',color=col[3], ecolor=col[3],capsize=3,uplims=uplims,label=r'UVOT_White')
if len(t_uvot_u_ul)!=0:
    ax.errorbar(t_uvot_u_ul,f_uvot_u_ul, yerr = ferr_uvot_u_ul,marker='s',ls='none', markersize='14',lw=2
                 ,color=col[4], ecolor=col[4],capsize=3,uplims=uplims,label=r'UVOT_U' )    
if len(t_uvot_v_ul)!=0:
    ax.errorbar(t_uvot_v_ul,f_uvot_v_ul,marker='s', yerr = ferr_uvot_v_ul,ls='none',
                 markersize='14',color=col[5], ecolor=col[5],capsize=3,uplims=uplims,label=r'UVOT_V' )

if len(t_uvot_b_ul)!=0:
    ax.errorbar(t_uvot_b_ul,f_uvot_b_ul,marker='s', yerr = ferr_uvot_b_ul,ls='none',
                 markersize='14',color=col[6], ecolor=col[6],capsize=3,uplims=uplims,label=r'UVOT_B')

if len(t_U_ul)!=0:
    ax.errorbar(t_U_ul,f_U_ul,marker='o', yerr = ferr_U_ul,ls='none',
                 markersize='14',color=col[7], ecolor=col[7],capsize=3,uplims=uplims,label='U' ) 
if len(t_B_ul)!=0:
    ax.errorbar(t_B_ul,f_B_ul,marker='o', yerr = ferr_B_ul,ls='none',
                 markersize='14',color=col[8], ecolor=col[8],capsize=3,uplims=uplims,label='B')    
if len(t_V_ul)!=0:
    ax.errorbar(t_V_ul,f_V_ul,marker='o', yerr = ferr_V_ul,ls='none',
                 markersize='14',color=col[9], ecolor=col[9],capsize=3,uplims=uplims,label='V')      
if len(t_R_ul)!=0:
    ax.errorbar(t_R_ul,f_R_ul,marker='o', yerr = ferr_R_ul,ls='none',
                 markersize='14',color=col[10], ecolor=col[10],capsize=3,uplims=uplims)     
if len(t_I_ul)!=0:
    ax.errorbar(t_I_ul,f_I_ul,marker='o', yerr = ferr_I_ul,ls='none',
                 markersize='14',color=col[11], ecolor=col[11],capsize=3,uplims=uplims,)
if len(t_J_ul)!=0:
    ax.errorbar(t_J_ul,f_J_ul,marker='o', yerr = ferr_J_ul,ls='none',
                 markersize='14',color=col[12], ecolor=col[12],capsize=3,uplims=uplims)      
if len(t_H_ul)!=0:
    ax.errorbar(t_H_ul,f_H_ul,marker='o', yerr = ferr_H_ul,ls='none',
                 markersize='14',color=col[13], ecolor=col[13],capsize=3,uplims=uplims)     
if len(t_K_ul)!=0:
    ax.errorbar(t_K_ul,f_K_ul,marker='o', yerr = ferr_K_ul,ls='none',
                 markersize='14',color=col[14], ecolor=col[14],capsize=3,uplims=uplims,)    
   
    
if len(t_u_ul)!=0:
    ax.errorbar(t_u_ul,f_u_ul,marker='v', yerr = ferr_u_ul,ls='none',
                 markersize='14',color=col[15], ecolor=col[15],capsize=3,uplims=uplims,label='u') 
if len(t_g_ul)!=0:
    ax.errorbar(t_g_ul,f_g_ul,marker='v', yerr = ferr_g_ul,ls='none',
                 markersize='14',color=col[16], ecolor=col[16],capsize=3,uplims=uplims,label='g')    
if len(t_r_ul)!=0:
    ax.errorbar(t_r_ul,f_r_ul,marker='v', yerr = ferr_r_ul,ls='none',
                 markersize='14',color=col[17], ecolor=col[17],capsize=3,uplims=uplims,label='r' )      
if len(t_i_ul)!=0:
    ax.errorbar(t_i_ul,f_i_ul,marker='v', yerr = ferr_i_ul,ls='none',
                 markersize='14',color=col[18], ecolor=col[18],capsize=3,uplims=uplims,label='i' )     
if len(t_z_ul)!=0:
    ax.errorbar(t_z_ul,f_z_ul,marker='v', yerr = ferr_z_ul,ls='none',
                 markersize='14',color=col[19], ecolor=col[19],capsize=3,uplims=uplims,label='z')

'''
#------------------------Running forest---------
ndim=6 #Number of parameters in the model
datafile='mltnst_dim6' # Name of the 'PyMultiNest' datafile
forest=datafile + "_ev.dat"
smpl=np.loadtxt(grb+'/'+forest)

smpl0 =smpl[:,0]
smpl1 = smpl[:,1]
smpl2 =smpl[:,2]
smpl3 = smpl[:,3] #this is p
smpl4 = smpl[:,4]
smpl5 = smpl[:,5]
'''
fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
#samples = np.loadtxt('mltnst_dim6_ev.dat')
#labels = ['k','a']
for i in range(ndim):
    ax1 = axes[i]
    ax1.plot(smpl[:, i], "k", alpha=0.5)
    ax1.set_xlim(0, len(smpl))
    ax1.set_ylabel(parameters[i])
    ax1.yaxis.set_label_coords(-0.1, 0.5)

#axes[-1].set_xlabel("step number");
plt.savefig(datafile+'_walkers.pdf',format='pdf')
plt.close()
'''
tot=len(smpl0)
burn = int(tot*3./5.)  #Taking 60% as burn-in

nlines=500
tsec=np.logspace(1.,8,100)
 
if len(t_uvot_uvw1)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_UVW1,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[0],ls='-',lw=1,alpha=0.1)  
                                
if len(t_uvot_uvw2)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_UVW2,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[1],ls='-',lw=1,alpha=0.1)  

if len(t_uvot_uvm2)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_UVM2,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[2],ls='-',lw=1,alpha=0.1)  
                                
if len(t_uvot_white)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_Wh,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[3],ls='-',lw=1,alpha=0.1) 
                
if len(t_uvot_u)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_U1,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[4],ls='-',lw=1,alpha=0.1)  
                                
if len(t_uvot_v)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_V1,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[5],ls='-',lw=1,alpha=0.1) 
                
if len(t_uvot_b)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_B1,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[6],ls='-',lw=1,alpha=0.1)  

if len(t_U)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_U,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[7],ls='-',lw=1,alpha=0.1)     
                
if len(t_B)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_B,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[8],ls='-',lw=1,alpha=0.1)  
                                
if len(t_V)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_V,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[9],ls='-',lw=1,alpha=0.1)     
                
if len(t_R)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_R,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[10],ls='-',lw=1,alpha=0.1)  
                                
if len(t_I)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_I,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[11],ls='-',lw=1,alpha=0.1)  
                
if len(t_J)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_J,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[12],ls='-',lw=1,alpha=0.1)  
                                
if len(t_H)!=0 or len(t_H_ul)!=0:
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_H,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[13],ls='-',lw=1,alpha=0.1)  
                
if len(t_K)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_K,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[14],ls='-',lw=1,alpha=0.1)  
                                
if len(t_u)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_u,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[15],ls='-',lw=1,alpha=0.1)  
                
if len(t_g)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_g,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[16],ls='-',lw=1,alpha=0.1)  
                                
if len(t_r)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_r,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[17],ls='-',lw=1,alpha=0.1)                  
                
if len(t_i)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_i,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[18],ls='-',lw=1,alpha=0.1)  
                                
if len(t_z)!=0 :
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_z,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[19],ls='-',lw=1,alpha=0.1)  
'''                
if len(t_uvot_uvw1)!=0 or len(t_uvot_uvw1_ul)!=0:
    for i in range(nlines):
        n=int(np.random.uniform(1,tot-burn))
        parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
        Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_U1,parmlist)
        if len(np.asarray([Fnutot]).ravel())>1:
              ax.plot(tsec,Fnutot,c=col[20],ls='-',lw=1,alpha=0.1)  
'''                                
for i in range(nlines):
    n=int(np.random.uniform(1,tot-burn))
    parmlist=(smpl0[n+burn],smpl1[n+burn],smpl2[n+burn],smpl3[n+burn],smpl4[n+burn],smpl5[n+burn])
    Fnutot=lcflux(tsec,np.ones(len(tsec))*nu_x,parmlist)
    if len(np.asarray([Fnutot]).ravel())>1:
        ax.plot(tsec,Fnutot,c=col[21],ls='-',lw=1,alpha=0.1)    
                
#ax2.legend(numpoints=1,prop={'size':28,},loc='upper right',labelspacing=0.1)
ax.legend(numpoints=1,prop={'size':20,},loc='best',labelspacing=0.3)
ax.set_title('GRB'+ grb)
plt.savefig('GRB_'+grb+'_fitting.pdf', format='pdf')
plt.tight_layout()
show()






