# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 11:39:24 2022

@author: HP
"""

# Here we import import the python libraries we will use
import matplotlib.pyplot as plt
# import matplotlib.colors as colors
# from mpl_toolkits.mplot3d import Axes3D
import numpy as np
# from matplotlib.transforms import Bbox
# import matplotlib

def time2semichord(time):
    return 2*Uinf*time/chord

def semichord2time(s):
    return s/2/Uinf*chord

# define properties of the system
dt=.1
time =np.arange(0,500,dt)
Uinf=1

# properties of the airfoil
chord=1 # chord of the airfoil
dCn_dalpha=2*np.pi # lift slope
alpha0=0 # alpha for which normal load is zero in steady flow

# pitching motion of the airfoil
k=.0 # reduced frequency of the pitching motion
omega=k*2*chord/Uinf # frequency of the piching motion
Amplitude_alpha=10/180*np.pi # amplitude of the pitching motion
alpha_t0=15/180*np.pi # alpha at time=0
alpha=Amplitude_alpha*np.sin(omega*time)+alpha_t0 # calculate alpha
dalpha_dt=np.gradient(alpha,time) # calculate the time derivative of alpha


# plunge motion of the airfoil
k_plg=.3 # reduced frequency of the plunge motion
omega_plg=k_plg*2*chord/Uinf # frequency of the plunge motion
Amplitude_plg=.3 # amplitude of the plunge motion
hplg=Amplitude_plg*np.sin(omega_plg*time) #position
dhplg_dt=np.gradient(hplg,time) # plunge velocity



# define the array semi-chord time scale
sarray = time2semichord(time)


# calculate quasi-steady alpha
alpha0 = 0 # we define the \alpha_{0}, zero for the case of an uncambered plate/airfoil
alphaqs = alpha + dalpha_dt*(chord/2)/Uinf - dhplg_dt/Uinf
dalphaqs_dt=np.gradient(alphaqs,time) # calculate the time derivative of the quasi-steady alpha

# calculate the coefficient of normal force assuming quasi-steady flow asuming potential flow
Cnormal_quasisteady = 2*np.pi*(alphaqs-alpha0)

# we plot the effective quasi-steady angle of attack \alpha_{qs}

# plot figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
plt.rcParams["font.family"] = "serif" # define font
plt.rcParams["mathtext.fontset"] = "dejavuserif"  # define font
cmap = plt.get_cmap('BuGn')  # define colormap
fig,ax = plt.subplots(figsize=[6,6]) # define pointers for the figure and axes
ax.plot(alpha*180/np.pi, alphaqs*180/np.pi,color='black', linewidth=1) # plot equivalent quasi-steady angle of attack
ax.set_xlabel(r'$\alpha (^\circ)$') # set x-label
ax.set_ylabel(r'$\alpha_{qs} (^\circ)$') # set y-label
# add arrows to indicate the direction of the cycle
parr1=ax.annotate('', xy=(17.5, 20), xytext=(10,12.5), 
            arrowprops=dict(color='black', shrink=0.05, width=.5, headwidth=3,headlength=4, linewidth=.2))
parr1=ax.annotate('', xy=(10, 7.5), xytext=(17.7,15), 
            arrowprops=dict(color='black', shrink=0.05, width=.5, headwidth=3,headlength=4, linewidth=.2))
plt.grid() # add a grid 
ax.set_xlim(0,30) # define limits of the axis
ax.set_ylim(0,30) # define limits of the axis
plt.tight_layout() # all elements of figure inside plot area
plt.show() # show figure

#filename = 'figures/alpha_quasi_steady' # define name of the figure to be saved
#fig.savefig(filename+'.svg', pad_inches = 0) # save figure
#fig.savefig(filename+'.pdf',pad_inches = 0) # save figure
#fig.savefig(filename+'.png', pad_inches = 0, dpi=300) # save figure

# we plot the quasi-steady normal-force coefficient

# plot figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
plt.rcParams["font.family"] = "serif" # define font
plt.rcParams["mathtext.fontset"] = "dejavuserif"  # define font
cmap = plt.get_cmap('BuGn')  # define colormap
fig,ax = plt.subplots(figsize=[6,6]) # define pointers for the figure and axes
ax.plot(alpha*180/np.pi, Cnormal_quasisteady,color='black', linewidth=1) # plot equivalent quasi-steady angle of attack
ax.set_xlabel(r'$\alpha (^\circ)$') # set x-label
ax.set_ylabel(r'$Cn_{qs} $') # set y-label
# add arrows to indicate the direction of the cycle
parr1=ax.annotate('', xy=(17.5, 20/360*4*np.pi**2), xytext=(10,12.5/360*4*np.pi**2), 
            arrowprops=dict(color='black', shrink=0.05, width=.5, headwidth=3,headlength=4, linewidth=.2))
parr1=ax.annotate('', xy=(10, 7.5/360*4*np.pi**2), xytext=(17.7,15/360*4*np.pi**2), 
            arrowprops=dict(color='black', shrink=0.05, width=.5, headwidth=3,headlength=4, linewidth=.2))
plt.grid() # add a grid 
ax.set_xlim(0,30) # define limits of the axis
ax.set_ylim(0,3) # define limits of the axis
plt.tight_layout() # all elements of figure inside plot area
plt.show() # show figure

#filename = 'figures/Cnormal_quasisteady' # define name of the figure to be saved
#fig.savefig(filename+'.svg', pad_inches = 0) # save figure
#fig.savefig(filename+'.pdf',pad_inches = 0) # save figure
#fig.savefig(filename+'.png', pad_inches = 0, dpi=300) # save figure

# determining X and Y terms for recursive marching formula for approximation of Duhamel's integral 
def duhamel_approx(Xi,Yi,delta_s,delta_alpha,order=2,A1=0.3,A2=0.7,b1=0.14,b2=0.53):
    # A1=0.165,A2=0.335,b1=0.0455,b2=0.3
    # determine the next values of X and Y, named Xip1 and Yip1
    if order==1:
        Xip1= Xi*np.exp(-b1*delta_s)+A1*delta_alpha
        Yip1= Yi*np.exp(-b2*delta_s)+A2*delta_alpha
    elif order==2:
        Xip1= Xi*np.exp(-b1*delta_s)+A1*delta_alpha*np.exp(-b1*delta_s/2)
        Yip1= Yi*np.exp(-b2*delta_s)+A2*delta_alpha*np.exp(-b2*delta_s/2)        
    else:
        Xip1= Xi*np.exp(-b1*delta_s)+A1*delta_alpha*((1+4*np.exp(-b1*delta_s/2)+np.exp(-b1*delta_s))/6)
        Yip1= Yi*np.exp(-b2*delta_s)+A2*delta_alpha*((1+4*np.exp(-b2*delta_s/2)+np.exp(-b2*delta_s))/6)
    
    return Xip1,Yip1


# define function for circulatory force, potential flow
def circulatory_normal_force(dCn_dalpha,alpha_equivalent,alpha0):
    return dCn_dalpha*(alpha_equivalent-alpha0)


# define arrays for X,Y and alpha_equivalent
Xarray=np.zeros(np.shape(time))
Yarray=np.zeros(np.shape(time))

# define the array of alpha_equivalent
alpha_equivalent=np.zeros(np.shape(time))
alpha_equivalent[0]=alphaqs[0]

# march solution in time for alpha_E
for i,val in enumerate(time[:-1]):
    Xarray[i+1],Yarray[i+1]=duhamel_approx(Xarray[i],Yarray[i],sarray[i+1]-sarray[i],alphaqs[i+1]-alphaqs[i])

alpha_equivalent=alphaqs-Xarray-Yarray


# plot solutions of test of duhamel_approx

# plot figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
plt.rcParams["font.family"] = "serif" # define font
plt.rcParams["mathtext.fontset"] = "dejavuserif" # define font
cmap = plt.get_cmap('BuGn') # define colormap
fig,ax = plt.subplots(figsize=[6,6]) # define pointers for figure and axes


#we will only plot the last cycle
Ncycles = np.floor(time[-1]*omega/(2*np.pi)) # determine number of cycles
n_of_cycle = time*omega/(2*np.pi) # calculate the phase of the different points of the cycle
i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1))) # index of start of cycle plotted
i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5))) # index of 180 degrees
i3=np.argmin(np.abs(n_of_cycle-(Ncycles))) # index of 360 degrees

# plot last cycle of the simulation, steady, quasi-steady and unsteady equivalent angle of attack 
ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1]), alpha[i1:i3]*180/np.pi,color='blue',linestyle='--', label=r'$\alpha$')
ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1]), alphaqs[i1:i3]*180/np.pi,color='red',linestyle='-.', label=r'$\alpha_{qs}$')
ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1]), alpha_equivalent[i1:i3]*180/np.pi,color='green',linestyle='-', label=r'$\alpha_{eq}$')
ax.set_xlabel('s semichords') # set x-label
ax.set_ylabel(r'$(^\circ)$') # set y-label
ax.set_xlim(0,2) # define limits of the axis
ax.set_ylim(0,30) # define limits of the axis
ax.grid() # add grid
ax.legend(loc='lower left')
plt.tight_layout() # all elements of figure inside plot area
plt.show() # show figure



#filename = 'figures/comparison_alpha_st_qs_circ' # define name of the figure to be saved
#fig.savefig(filename+'.svg', pad_inches = 0) # save figure
#fig.savefig(filename+'.pdf',pad_inches = 0) # save figure
#fig.savefig(filename+'.png', pad_inches = 0, dpi=300) # save figure


# plot figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
cmap = plt.get_cmap('BuGn')
fig,ax = plt.subplots(figsize=[6,6])

#we will only plot the last cycle
Ncycles = np.floor(time[-1]*omega/(2*np.pi))
n_of_cycle = time*omega/(2*np.pi) # calculate the phase of the different points of the cycle
i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1))) # index of start of cycle plotted
i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5))) # index of 180 degrees
i3=np.argmin(np.abs(n_of_cycle-(Ncycles))) # index of 360 degrees


ax.plot(alpha[i1:i3]*180/np.pi, alphaqs[i1:i3]*180/np.pi,color='blue',linestyle='--', label=r'$\alpha_{qs}$')
ax.plot(alpha[i1:i3]*180/np.pi, alpha_equivalent[i1:i3]*180/np.pi,color='red',linestyle='dashdot', label=r'$\alpha_{eq}$')


# we will plot arrows to see the direction of the cycle
scale_arrow=3 # scale od arrow
dx = (alpha[i1]-alpha[i1-1]) # dx of arrow
dy = (alphaqs[i1]-alphaqs[i1-1])  # dy of arrow   
ax.arrow(alpha[i1]*180/np.pi, alphaqs[i1]*180/np.pi, 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='blue', width=scale_arrow*.1, shape='left') # plot arrow at 0 degrees of cycle
dx = (alpha[i2]-alpha[i2-1]) # dx of arrow
dy = (alphaqs[i2]-alphaqs[i2-1])  # dy of arrow   
ax.arrow(alpha[i2]*180/np.pi, alphaqs[i2]*180/np.pi, 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='blue', width=scale_arrow*.1, shape='left') # plot arrow at 0 degrees of cycle


dx = (alpha[i1]-alpha[i1-1]) # dx of arrow
dy = (alpha_equivalent[i1]-alpha_equivalent[i1-1])  # dy of arrow   
ax.arrow(alpha[i1]*180/np.pi, alpha_equivalent[i1]*180/np.pi, 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='red', width=scale_arrow*.1, shape='left') # plot arrow at 0 degrees of cycle
dx = (alpha[i2]-alpha[i2-1]) # dx of arrow
dy = (alpha_equivalent[i2]-alpha_equivalent[i2-1])  # dy of arrow   
ax.arrow(alpha[i2]*180/np.pi, alpha_equivalent[i2]*180/np.pi, 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='red', width=scale_arrow*.1, shape='left') # plot arrow at 0 degrees of cycle

# ax.set_aspect(aspect=40.0)
ax.set_xlabel(r'$\alpha (^\circ)$')
ax.set_ylabel(r'$ (^\circ)$')
# ax.set_xlim(0,time.max())
ax.legend(loc='lower right')
plt.grid()
plt.tight_layout() # all elements of figure inside plot area

plt.show()

#filename = 'figures/comparison_cycle_alpha_qs_circ' # define name of the figure to be saved
#fig.savefig(filename+'.svg', pad_inches = 0) # save figure
#fig.savefig(filename+'.pdf',pad_inches = 0) # save figure
#fig.savefig(filename+'.png', pad_inches = 0, dpi=300) # save figure


# plot solutions of test of duhamel_approx

# plot figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
plt.rcParams["font.family"] = "serif" # define font
plt.rcParams["mathtext.fontset"] = "dejavuserif" # define font
cmap = plt.get_cmap('BuGn') # define colormap
fig,ax = plt.subplots(figsize=[6,6]) # define pointers for figure and axes


#we will only plot the last cycle
Ncycles = np.floor(time[-1]*omega/(2*np.pi)) # determine number of cycles
n_of_cycle = time*omega/(2*np.pi) # calculate the phase of the different points of the cycle
i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1))) # index of start of cycle plotted
i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5))) # index of 180 degrees
i3=np.argmin(np.abs(n_of_cycle-(Ncycles))) # index of 360 degrees

# plot last cycle of the simulation, steady, quasi-steady and unsteady normal force coefficient  
ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1]), 
        circulatory_normal_force(2*np.pi,alpha[i1:i3],0),color='blue',linestyle='--', label=r'$Cn_{st}$')
ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1]), 
        circulatory_normal_force(2*np.pi,alphaqs[i1:i3],0),color='red',linestyle='-.', label=r'$Cn_{qs}$')
ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1]), 
        circulatory_normal_force(2*np.pi,alpha_equivalent[i1:i3],0),color='green',linestyle='-', label=r'$Cn_{c}$')
ax.set_xlabel('s semichords') # set x-label
ax.set_ylabel(r'$Cn$') # set y-label
ax.set_xlim(0,2) # define limits of the axis
ax.set_ylim(0,3) # define limits of the axis
ax.grid() # add grid
ax.legend(loc='lower left')
plt.tight_layout() # all elements of figure inside plot area
plt.show() # show figure



#filename = 'figures/comparison_Cn_st_qs_circ' # define name of the figure to be saved
#fig.savefig(filename+'.svg', pad_inches = .0) # save figure
#fig.savefig(filename+'.pdf',pad_inches = .0) # save figure
#fig.savefig(filename+'.png', pad_inches = .0, dpi=300) # save figure


# plot figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
cmap = plt.get_cmap('BuGn')
fig,ax = plt.subplots(figsize=[6,6])

#we will only plot the last cycle
Ncycles = np.floor(time[-1]*omega/(2*np.pi))
n_of_cycle = time*omega/(2*np.pi) # calculate the phase of the different points of the cycle
i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1))) # index of start of cycle plotted
i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5))) # index of 180 degrees
i3=np.argmin(np.abs(n_of_cycle-(Ncycles))) # index of 360 degrees


ax.plot(alpha[i1:i3]*180/np.pi, 
        circulatory_normal_force(2*np.pi,alphaqs[i1:i3],0),color='blue',linestyle='--', label=r'$Cn_{qs}$')
ax.plot(alpha[i1:i3]*180/np.pi, 
        circulatory_normal_force(2*np.pi,alpha_equivalent[i1:i3],0),color='red',linestyle='dashdot', label=r'$Cn_{c}$')


# we will plot arrows to see the direction of the cycle
scale_arrow=3 # scale od arrow
dx = (alpha[i1]-alpha[i1-1])*180/np.pi # dx of arrow
dy = circulatory_normal_force(2*np.pi,(alphaqs[i1]-alphaqs[i1-1]),0)  # dy of arrow   
ax.arrow(alpha[i1]*180/np.pi, circulatory_normal_force(2*np.pi,alphaqs[i1],0), 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='blue', width=scale_arrow*.02,shape='left') # plot arrow at 0 degrees of cycle
dx = (alpha[i2]-alpha[i2-1])*180/np.pi # dx of arrow
dy = circulatory_normal_force(2*np.pi,(alphaqs[i2]-alphaqs[i2-1]),0)  # dy of arrow   
ax.arrow(alpha[i2]*180/np.pi, circulatory_normal_force(2*np.pi,alphaqs[i2],0), 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='blue', width=scale_arrow*.02,shape='left') # plot arrow at 0 degrees of cycle


dx = (alpha[i1]-alpha[i1-1])*180/np.pi # dx of arrow
dy = circulatory_normal_force(2*np.pi,(alpha_equivalent[i1]-alpha_equivalent[i1-1]),0)  # dy of arrow   
ax.arrow(alpha[i1]*180/np.pi, circulatory_normal_force(2*np.pi,alpha_equivalent[i1],0), 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='red', width=scale_arrow*.02, shape='left') # plot arrow at 0 degrees of cycle
dx = (alpha[i2]-alpha[i2-1])*180/np.pi # dx of arrow
dy = circulatory_normal_force(2*np.pi,(alpha_equivalent[i2]-alpha_equivalent[i2-1]),0)  # dy of arrow   
ax.arrow(alpha[i2]*180/np.pi, circulatory_normal_force(2*np.pi,alpha_equivalent[i2],0), 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='red', width=scale_arrow*.02, shape='left') # plot arrow at 0 degrees of cycle

# ax.set_aspect(aspect=40.0)
ax.set_xlabel(r'$\alpha (^\circ)$')
ax.set_ylabel(r'$Cn$')
# ax.set_xlim(0,time.max())
ax.set_ylim(0,3)
ax.legend(loc='lower right')
plt.grid()
plt.tight_layout() # all elements of figure inside plot area

plt.show()

#filename = 'figures/comparison_cycle_Cn_qs_circ' # define name of the figure to be saved
#fig.savefig(filename+'.svg', pad_inches = 0.) # save figure
#fig.savefig(filename+'.pdf',pad_inches = 0.) # save figure
#fig.savefig(filename+'.png', pad_inches = 0., dpi=300) # save figure



# deficiency function for non-circulatory normal force

def deficiency_function(Dnoncirc_i,delta_dalpha_dt,delta_t,chord,asound=343,kalpha=0.75):
    # a sound is the speed of sound
    TI=chord/asound
    Dnoncirc_ip1 = Dnoncirc_i*np.exp(-delta_t/(kalpha*TI))+delta_dalpha_dt*np.exp(-delta_t/(2*kalpha*TI))
    return Dnoncirc_ip1

# non-circulatory normal force
def non_circulatory_normal_force(dalpha_dt,chord,Uinf,Dnoncirc,kalpha=0.75):
    return 4*kalpha*chord/Uinf*(dalpha_dt-Dnoncirc)
    


# define arrays for Dnoncirc, the deficiency function for non-circulatory loading
Dnoncirc=np.zeros(np.shape(time))

# march solution in time
for i,val in enumerate(time[:-1]):
    Dnoncirc[i+1]=deficiency_function(Dnoncirc[i],dalphaqs_dt[i+1]-dalphaqs_dt[i],dt,chord)
    
    
Cnormal_circ = circulatory_normal_force(dCn_dalpha,alpha_equivalent,alpha0)
Cnormal_noncirc =  non_circulatory_normal_force(dalphaqs_dt,chord,Uinf,Dnoncirc)
Cnormal_p = Cnormal_circ+Cnormal_noncirc

# plot figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
plt.rcParams["font.family"] = "serif" # define font
plt.rcParams["mathtext.fontset"] = "dejavuserif" # define font
cmap = plt.get_cmap('BuGn') # define colormap
fig,ax = plt.subplots(figsize=[6,6]) # define pointers for figure and axes


#we will only plot the last cycle
Ncycles = np.floor(time[-1]*omega/(2*np.pi)) # determine number of cycles
n_of_cycle = time*omega/(2*np.pi) # calculate the phase of the different points of the cycle
i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1))) # index of start of cycle plotted
i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5))) # index of 180 degrees
i3=np.argmin(np.abs(n_of_cycle-(Ncycles))) # index of 360 degrees

# plot last cycle of the simulation, normal force coefficient  
ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1]), 
        Cnormal_p[i1:i3],color='blue',linestyle='--', label=r'$Cn_{p}$')
ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1]), 
        Cnormal_circ[i1:i3],color='red',linestyle='-.', label=r'$Cn_{c}$')
ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1]), 
        Cnormal_noncirc[i1:i3],color='green',linestyle='-', label=r'$Cn_{nc}$')
ax.set_xlabel('s semichords') # set x-label
ax.set_ylabel(r'$Cn$') # set y-label
ax.set_xlim(0,2) # define limits of the axis
ax.set_ylim(-1,3) # define limits of the axis
ax.grid() # add grid
ax.legend(loc='lower left')
plt.tight_layout() # all elements of figure inside plot area

plt.show() # show figure



#filename = 'figures/comparison_Cn_p_circ_noncirc' # define name of the figure to be saved
#fig.savefig(filename+'.svg', pad_inches = 0.) # save figure
#fig.savefig(filename+'.pdf',pad_inches = 0.) # save figure
#fig.savefig(filename+'.png', pad_inches = 0., dpi=300) # save figure

# plot figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
cmap = plt.get_cmap('BuGn')
fig,ax = plt.subplots(figsize=[6,6])

#we will only plot the last cycle
Ncycles = np.floor(time[-1]*omega/(2*np.pi))
n_of_cycle = time*omega/(2*np.pi) # calculate the phase of the different points of the cycle
i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1))) # index of start of cycle plotted
i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5))) # index of 180 degrees
i3=np.argmin(np.abs(n_of_cycle-(Ncycles))) # index of 360 degrees


ax.plot(alpha[i1:i3]*180/np.pi, 
        Cnormal_p[i1:i3],color='blue',linestyle='--', label=r'$Cn_{p}$')
ax.plot(alpha[i1:i3]*180/np.pi, 
        Cnormal_circ[i1:i3],color='red',linestyle='dashdot', label=r'$Cn_{c}$')
ax.plot(alpha[i1:i3]*180/np.pi, 
        Cnormal_noncirc[i1:i3],color='green',linestyle='dashdot', label=r'$Cn_{nc}$')

# we will plot arrows to see the direction of the cycle
scale_arrow=3 # scale od arrow
dx = (alpha[i1]-alpha[i1-1])*180/np.pi # dx of arrow
dy = Cnormal_p[i1]-Cnormal_p[i1-1]  # dy of arrow   
ax.arrow(alpha[i1]*180/np.pi, Cnormal_p[i1], 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='blue', width=scale_arrow*.02,shape='left') # plot arrow at 0 degrees of cycle
dx = (alpha[i2]-alpha[i2-1])*180/np.pi # dx of arrow
dy = Cnormal_p[i2]-Cnormal_p[i2-1]  # dy of arrow   
ax.arrow(alpha[i2]*180/np.pi, Cnormal_p[i2], 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='blue', width=scale_arrow*.02,shape='left') # plot arrow at 0 degrees of cycle


dx = (alpha[i1]-alpha[i1-1])*180/np.pi # dx of arrow
dy = Cnormal_circ[i1]-Cnormal_circ[i1-1]  # dy of arrow      
ax.arrow(alpha[i1]*180/np.pi, Cnormal_circ[i1], 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='red', width=scale_arrow*.02, shape='left') # plot arrow at 0 degrees of cycle
dx = (alpha[i2]-alpha[i2-1])*180/np.pi # dx of arrow
dy = Cnormal_circ[i2]-Cnormal_circ[i2-1]  # dy of arrow   
ax.arrow(alpha[i2]*180/np.pi, Cnormal_circ[i2], 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='red', width=scale_arrow*.02, shape='left') # plot arrow at 0 degrees of cycle

dx = (alpha[i1]-alpha[i1-1])*180/np.pi # dx of arrow
dy = Cnormal_noncirc[i1]-Cnormal_noncirc[i1-1]  # dy of arrow      
ax.arrow(alpha[i1]*180/np.pi, Cnormal_noncirc[i1], 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='green', width=scale_arrow*.02, shape='left') # plot arrow at 0 degrees of cycle
dx = (alpha[i2]-alpha[i2-1])*180/np.pi # dx of arrow
dy = Cnormal_noncirc[i2]-Cnormal_noncirc[i2-1]  # dy of arrow   
ax.arrow(alpha[i2]*180/np.pi, Cnormal_noncirc[i2], 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='green', width=scale_arrow*.02, shape='left') # plot arrow at 0 degrees of cycle

# ax.set_aspect(aspect=40.0)
ax.set_xlabel(r'$\alpha (^\circ)$')
ax.set_ylabel(r'$Cn$')
ax.set_ylim(-1,3)
ax.set_xlim(0,25)

ax.legend(loc='upper left')
plt.grid()

plt.tight_layout() # all elements of figure inside plot area

plt.show()

#filename = 'figures/comparison_cycle_Cn_p_circ_noncirc' # define name of the figure to be saved
#fig.savefig(filename+'.svg', pad_inches = 0.) # save figure
#fig.savefig(filename+'.pdf',pad_inches = 0.) # save figure
#fig.savefig(filename+'.png', pad_inches = 0., dpi=300) # save figure


# definition of a function for the trailing edge separation point "f"

def f_trailing_edge_separation_point(alpha, a1=7,a2=15,a3=21):
    # receives alpha in radians, converts to degrees
    alphadeg = alpha*180/np.pi
    if alphadeg<=a1:
        f=1
    elif ((alphadeg>a1) and (alphadeg<=a2)):
        f= 1 - .8*((alphadeg-a1)/(a2-a1))
    elif ((alphadeg>a2) and (alphadeg<a3)):
        f= .2 *(1- ((alphadeg-a2)/(a3-a2))**.3)
    else:
        f=0
    return f


# test polar with trailing edge separation point "f"
alpha_polar=np.arange(-5,30,.1)
CNsep=np.zeros(np.shape(alpha_polar))
for i,val in enumerate(alpha_polar):
    CNsep[i]= dCn_dalpha*(((1+np.sqrt(f_trailing_edge_separation_point(alpha_polar[i]*np.pi/180)))/2)**2)*(alpha_polar[i]*np.pi/180-alpha0)
# we will also calculate the steady value of Cn to compare with later results
CNsteady=np.zeros(np.shape(time))
for i,val in enumerate(time):
    CNsteady[i]= dCn_dalpha*(((1+np.sqrt(f_trailing_edge_separation_point(alpha[i])))/2)**2)*(alpha[i]-alpha0)



# plot figure
#plt.rcParams["font.family"] = "serif"
#plt.rcParams["mathtext.fontset"] = "dejavuserif"
#cmap = plt.get_cmap('BuGn')
#fig,ax = plt.subplots(figsize=[6,6])
#ax.plot(alpha_polar, dCn_dalpha*(alpha_polar*np.pi/180-alpha0),color='green',linestyle='dashdot', label=r'potential flow')
#ax.plot(alpha_polar, CNsep,color='red',linestyle='-', label=r'viscous flow')
# ax.set_aspect(aspect=40.0)
#ax.set_xlabel(r'$\alpha (^\circ)$')
#ax.set_ylabel(r'$C_n$')
#ax.set_xlim(0,25)
#ax.legend(loc='lower right')

plt.tight_layout() # all elements of figure inside plot area

plt.grid()


plt.show()

# we will now determine the effect of the pressure lag in terms of onset of the separation point

def pressure_lag_deficiency(Dpress_i,delta_s, delta_CNpot, Tp=1.7):
    return Dpress_i*np.exp(-delta_s/Tp)+ delta_CNpot*np.exp(-delta_s/2/Tp) 
    
# we need to define an array for the pressure lag deficiency function
Dpress = np.zeros(np.shape(time))

# we will now do the time marching to solve for the pressure lag deficiency function 
for i,val in enumerate(time[:-1]):
    Dpress[i+1] = pressure_lag_deficiency(Dpress[i],sarray[i+1]-sarray[i], Cnormal_p[i+1]-Cnormal_p[i])

# we now determine the normal force coefficient due to the pressure lag
Cnormal_prime = Cnormal_p-Dpress

# an based on this Cnormal_prime, we determine a new equivalent angle of attack 
# to determine the onset of trailing edge separation
alpha_f = Cnormal_prime/dCn_dalpha+alpha0

# we use this equivalent angle of attack alpha_f to determine a new trailign edge separation point effect f_prime
fprime = np.zeros(np.shape(time))
for i,val in enumerate(time):
    fprime[i]  = f_trailing_edge_separation_point(alpha_f[i])
    
    
# plot figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
cmap = plt.get_cmap('BuGn')
fig,ax = plt.subplots(figsize=[6,6])

#we will only plot the last cycle
Ncycles = np.floor(time[-1]*omega/(2*np.pi))
n_of_cycle = time*omega/(2*np.pi) # calculate the phase of the different points of the cycle
i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1))) # index of start of cycle plotted
i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5))) # index of 180 degrees
i3=np.argmin(np.abs(n_of_cycle-(Ncycles))) # index of 360 degrees


ax.plot(alpha[i1:i3]*180/np.pi, alpha_f[i1:i3]*180/np.pi,color='blue',linestyle='--', label=r'$\alpha_{f}$')
ax.plot(alpha[i1:i3]*180/np.pi, alpha_equivalent[i1:i3]*180/np.pi,color='red',linestyle='dashdot', label=r'$\alpha_{eq}$')


# we will plot arrows to see the direction of the cycle
scale_arrow=3 # scale od arrow
dx = (alpha[i1]-alpha[i1-1]) # dx of arrow
dy = (alphaqs[i1]-alphaqs[i1-1])  # dy of arrow   
ax.arrow(alpha[i1]*180/np.pi, alpha_f[i1]*180/np.pi, 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='blue', width=scale_arrow*.1, shape='left') # plot arrow at 0 degrees of cycle
dx = (alpha[i2]-alpha[i2-1]) # dx of arrow
dy = (alphaqs[i2]-alphaqs[i2-1])  # dy of arrow   
ax.arrow(alpha[i2]*180/np.pi, alpha_f[i2]*180/np.pi, 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='blue', width=scale_arrow*.1, shape='left') # plot arrow at 0 degrees of cycle


dx = (alpha[i1]-alpha[i1-1]) # dx of arrow
dy = (alpha_equivalent[i1]-alpha_equivalent[i1-1])  # dy of arrow   
ax.arrow(alpha[i1]*180/np.pi, alpha_equivalent[i1]*180/np.pi, 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='red', width=scale_arrow*.1, shape='left') # plot arrow at 0 degrees of cycle
dx = (alpha[i2]-alpha[i2-1]) # dx of arrow
dy = (alpha_equivalent[i2]-alpha_equivalent[i2-1])  # dy of arrow   
ax.arrow(alpha[i2]*180/np.pi, alpha_equivalent[i2]*180/np.pi, 
             scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
              color='red', width=scale_arrow*.1, shape='left') # plot arrow at 0 degrees of cycle

# ax.set_aspect(aspect=40.0)
ax.set_xlabel(r'$\alpha (^\circ)$')
ax.set_ylabel(r'$ (^\circ)$')
# ax.set_xlim(0,time.max())
ax.legend(loc='lower right')
plt.grid()

plt.tight_layout() # all elements of figure inside plot area

plt.show()

#filename = 'figures/comparison_cycle_alpha_circ_f' # define name of the figure to be saved
#fig.savefig(filename+'.svg', pad_inches = 0.) # save figure
#fig.savefig(filename+'.pdf',pad_inches = 0.) # save figure
#fig.savefig(filename+'.png', pad_inches = 0., dpi=300) # save figure


## we will now implement a delay function for separation point for boundary layer lag effects

def boundary_layer_lag_deficiency(Dbl_i,delta_s, delta_fprime, Tf=3.0):
    return Dbl_i*np.exp(-delta_s/Tf)+ delta_fprime*np.exp(-delta_s/2/Tf) 
    
# we need to define an array for the boundary layer lag deficiency function
Dbl = np.zeros(np.shape(time))

# we will now do the time marching to solve for the boundary layer lag deficiency function 
for i,val in enumerate(time[:-1]):
    Dbl[i+1] = boundary_layer_lag_deficiency(Dbl[i],sarray[i+1]-sarray[i], fprime[i+1]-fprime[i])

# we now determine the a new expression of fprimeprime due to the boundary layer lag
fprimeprime = fprime-Dbl


## we can now determine the normal force due to trailing edge boundary layer separation

Cnormal_f = dCn_dalpha*((1+np.sqrt(fprimeprime))/2)**2*(alpha_equivalent-alpha0)+ Cnormal_noncirc



# plot figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
cmap = plt.get_cmap('BuGn')


#we will only plot the last cycle
Ncycles = np.floor(time[-1]*omega/(2*np.pi))
n_of_cycle = time*omega/(2*np.pi) # calculate the phase of the different points of the cycle
i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1))) # index of start of cycle plotted
i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5))) # index of 180 degrees
i3=np.argmin(np.abs(n_of_cycle-(Ncycles))) # index of 360 degrees


fig,ax = plt.subplots(figsize=[10,6])
ax.plot(alpha[i1:i3]*180/np.pi, Cnormal_f[i1:i3],color='black',linestyle='-', label=r'$Cn_f$')
ax.plot(alpha[i1:i3]*180/np.pi, Cnormal_p[i1:i3],color='red',linestyle='dashdot', label=r'$Cn_p$')
ax.plot(alpha[i1:i3]*180/np.pi, Cnormal_circ[i1:i3],color='blue',linestyle='-.', label=r'$Cn_{c}$')
ax.plot(alpha[i1:i3]*180/np.pi, Cnormal_noncirc[i1:i3],color='green',linestyle='--', label=r'$Cn_{nc}$')
ax.plot(alpha[i1:i3]*180/np.pi, CNsteady[i1:i3],color='grey',linestyle='--', label=r'$Cn_{st}$')
# ax.set_aspect(aspect=40.0)
ax.set_xlabel(r'$\alpha (^\circ)$')
ax.set_ylabel(r'$Cn$')
ax.set_xlim(0,40)
ax.legend(loc='lower right')
plt.grid()

plt.tight_layout() # all elements of figure inside plot area

plt.show()

#filename = 'figures/comparison_cycle_Cn_f_p_circ_noncirc_st' # define name of the figure to be saved
#fig.savefig(filename+'.svg', pad_inches = 0.) # save figure
#fig.savefig(filename+'.pdf',pad_inches = 0.) # save figure
#fig.savefig(filename+'.png', pad_inches = 0., dpi=300) # save figure


# we will now setup the non-dimensional vortex-time parameter vortime
# we will setup an equation for vortime, integrating in time

def vortime_function(vortime_i,delta_s,delta_alphaqs, Cnormal_prime, CN1=1.0093):
    if Cnormal_prime>CN1:
        vortime_ip1 = vortime_i + 0.45*delta_s
    else:
        if (delta_alphaqs<0 and vortime_i>0):
            vortime_ip1 = vortime_i + 0.45*delta_s
        else:
            vortime_ip1 = 0 
    
    return vortime_ip1



# we need to define an array for the non-dimensional vortex-time parameter vortime
vortime = np.zeros(np.shape(time))

# we will now do the time marching to solve for the non-dimensional vortex-time parameter vortime
for i,val in enumerate(time[:-1]):
    vortime[i+1] = vortime_function(vortime[i],sarray[i+1]-sarray[i],dalphaqs_dt[i+1]-dalphaqs_dt[i], Cnormal_prime[i])


# # plot figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
cmap = plt.get_cmap('BuGn')
fig,ax = plt.subplots(figsize=[6,6])
ax.plot(time, vortime,color='black',linestyle='-', label=r'$\tau_v$')
# # ax.set_aspect(aspect=40.0)
secax = ax.secondary_xaxis('top', functions=(time2semichord, semichord2time))
secax.set_xlabel('s semichords')
ax.set_xlabel(r'$t_{(s)}$')
ax.set_ylabel(r'$\tau_v$')
# # ax.set_xlim(0,time.max())
ax.legend(loc='lower right')
plt.grid()


# plt.show()


## we will now define the added normal force due to the presence of the leading edge vortex

# first, we determine the increment in normal force due to the presence of a vortex 'Cvortex' as a function 
# of the circulatory normal force and the separation location

Cvortex=Cnormal_circ*(1 - ( ((1+np.sqrt(fprimeprime))/2)**2  ) )

# we will now define the function for decay of the cumulative normal force due to the presence of the leading edge vortex

def leading_edge_vortex_normal_force(Cnormal_vortex_i,delta_s,delta_Cvortex,vortime,TVL=11,TV=6):
    if (vortime>0.001 and vortime<TVL):
        Cnormal_vortex_ip1=Cnormal_vortex_i*np.exp(-delta_s/TV)+ delta_Cvortex*np.exp(-delta_s/2/TV)      
    else:    
        Cnormal_vortex_ip1=Cnormal_vortex_i*np.exp(-delta_s/TV) 
    return Cnormal_vortex_ip1

# We will now solve for the cumulative normal force due to the leadign edge vortex by marching in time.
# First, we will define the array

Cnormal_vortex = np.zeros(np.shape(time))
Cnormal_vortex[0]=Cvortex[0]

for i,val in enumerate(time[:-1]):
    Cnormal_vortex[i+1] = leading_edge_vortex_normal_force(Cnormal_vortex[i],sarray[i+1]-sarray[i], 
                                                           Cvortex[i+1]-Cvortex[i],vortime[i])



# calculating total load

Cnormal_total = Cnormal_f + Cnormal_vortex

# plot figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
cmap = plt.get_cmap('BuGn')


#we will only plot the last cycle
Ncycles = np.floor(time[-1]*omega/(2*np.pi))
n_of_cycle = time*omega/(2*np.pi) # calculate the phase of the different points of the cycle
i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1))) # index of start of cycle plotted
i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5))) # index of 180 degrees
i3=np.argmin(np.abs(n_of_cycle-(Ncycles))) # index of 360 degrees


fig,ax = plt.subplots(figsize=[6,6])
ax.plot(alpha[i1:i3]*180/np.pi, Cnormal_total[i1:i3],color='black',linestyle='-', linewidth=1.5, label=r'$Cn_t$')
ax.plot(alpha[i1:i3]*180/np.pi, Cnormal_vortex[i1:i3],color='cyan',linestyle='-',linewidth=.75, label=r'$Cn_v$')
ax.plot(alpha[i1:i3]*180/np.pi, Cnormal_f[i1:i3],color='magenta',linestyle='--',linewidth=.75, label=r'$Cn_f$')
ax.plot(alpha[i1:i3]*180/np.pi, Cnormal_p[i1:i3],color='red',linestyle='dashdot',linewidth=.75, label=r'$Cn_p$')
ax.plot(alpha[i1:i3]*180/np.pi, Cnormal_circ[i1:i3],color='blue',linestyle='-.',linewidth=.75, label=r'$Cn_{c}$')
ax.plot(alpha[i1:i3]*180/np.pi, Cnormal_noncirc[i1:i3],color='green',linestyle='--',linewidth=.75, label=r'$Cn_{nc}$')
ax.plot(alpha[i1:i3]*180/np.pi, CNsteady[i1:i3],color='grey',linestyle='--',linewidth=.75, label=r'$Cn_{st}$')
# ax.set_aspect(aspect=40.0)
ax.set_xlabel(r'$\alpha (^\circ)$')
ax.set_ylabel(r'$Cn$')
ax.set_xlim(0,30)
ax.legend(loc='upper left')
plt.grid()

plt.tight_layout() # all elements of figure inside plot area

plt.show()

#filename = 'figures/comparison_cycle_Cn_t_v_f_p_circ_noncirc_st' # define name of the figure to be saved
#fig.savefig(filename+'.svg', pad_inches = 0) # save figure
#fig.savefig(filename+'.pdf',pad_inches = 0) # save figure
#fig.savefig(filename+'.png', pad_inches = 0, dpi=300) # save figure



Ncycles = np.floor(time[-1]*omega/2/np.pi)
# print('Number of cycles ',Ncycles)
ind1=np.argmin(np.abs((time*omega/2/np.pi)-(Ncycles-1)))
ind2=np.argmin(np.abs((time*omega/2/np.pi)-Ncycles))+1


alphacycle=alpha[ind1:ind2]
Cncycle=Cnormal_total[ind1:ind2]
timecycle=time[ind1:ind2]
fprimeprimecycle=fprimeprime[ind1:ind2]

# determine point of trailing edge separation
i1= np.where(fprimeprime[ind1:ind2]>.99)
i1= np.where(alphacycle == np.max(alphacycle[i1]))
p1 = [alphacycle[i1]*180/np.pi, Cncycle[i1]]
#print(p1)

# determine point leading edge separation
i2= np.where(vortime[ind1:ind2]<.01)
i2= np.where(alphacycle == np.max(alphacycle[i2]))
p2 = [alphacycle[i2]*180/np.pi, Cncycle[i2]]
# print(p2)


# determine point of leading edge vortex shed in the wake and contribution due to the leading-edge vortex starts to decay
i3= np.where(Cnormal_vortex[ind1:ind2]==np.max(Cnormal_vortex[ind1:ind2]))
i3= np.where(alphacycle == np.max(alphacycle[i3]))
p3 = [alphacycle[i3]*180/np.pi, Cncycle[i3]]
# print(p2)


# determine point where contribution by leading edge vortex is small, the flow is dominated by the trailing edge separation
i4= np.where(Cnormal_vortex[ind1:ind2]<0.05)
i4= np.where(alphacycle == np.max(alphacycle[i4]))
p4 = [alphacycle[i4]*180/np.pi, Cncycle[i4]]
# print(p4)


# determine point of reattachment of the flow
i5= np.where(np.gradient(alpha[ind1:ind2])<0)
# print(i6[0][-1])
i5= np.where(fprimeprimecycle[0:i5[0][-1]]<0.95)
# print(i6)
p5 = [alphacycle[i5[0][-1]]*180/np.pi, Cncycle[i5[0][-1]]]
# print(p5)



# plot figure
plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
cmap = plt.get_cmap('BuGn')
fig,ax = plt.subplots(figsize=[6,6])
ax.plot(alpha[ind1:ind2]*180/np.pi, Cnormal_total[ind1:ind2],color='blue',linestyle='-', label=r'Unsteady')
ax.plot(alpha_polar, CNsep,color='green',linestyle='--', label=r'Steady')
# ax.set_aspect(aspect=40.0)
ax.set_xlabel(r'$\alpha (^\circ)$')
ax.set_ylabel(r'$C_n$')
ax.set_xlim(0,25)
ax.set_ylim(0,2.0)
ax.legend(loc='lower right')
plt.grid()

# annotate dynamic stall points
bbox = dict(boxstyle="round", fc="0.9")
p1an=ax.annotate('1', xy=(p1[0], p1[1]), xytext=(p1[0]-5, p1[1]+.1),bbox=bbox,
            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=6))
p2an=ax.annotate('2', xy=(p2[0], p2[1]), xytext=(p2[0]-4, p2[1]+.3),bbox=bbox,
            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=6))
p3an=ax.annotate('3', xy=(p3[0], p3[1]), xytext=(p3[0]+2, p3[1]+.3),bbox=bbox,
            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=6))
p4an=ax.annotate('4', xy=(p4[0], p4[1]), xytext=(p4[0]-5, p4[1]+.1),bbox=bbox,
            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=6))
p5an=ax.annotate('5', xy=(p5[0], p5[1]), xytext=(p5[0]-2, p5[1]+.1),bbox=bbox,
            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=6))

parr1=ax.annotate('', xy=(15, 1.3), xytext=(11,1),
            arrowprops=dict(color='grey', shrink=0.05, width=.5, headwidth=3,headlength=4))
parr2=ax.annotate('', xy=(23, 0.8), xytext=(21,1.2),
            arrowprops=dict(color='grey', shrink=0.05, width=.5, headwidth=3,headlength=4))
parr3=ax.annotate('', xy=(12, 0.75), xytext=(20,0.6),
            arrowprops=dict(color='grey', shrink=0.05, width=.5, headwidth=3,headlength=4))

plt.tight_layout() # all elements of figure inside plot area
plt.show()



#filename = 'figures/dynamic_stall_cycle'
#fig.savefig(filename+'.svg', pad_inches = 0) # save figure
#fig.savefig(filename+'.pdf',pad_inches = 0) # save figure
#fig.savefig(filename+'.png', pad_inches = 0, dpi=300) # save figure


