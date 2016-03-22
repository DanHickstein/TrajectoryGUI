# This program plots the trajectories from two-color circularly polarized fields.
# You need to install python. If you're using a Mac with OS 10.9, all you need is already installed
# DanHickstein@gmail.com

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons

def plot_color_line(ax,x,y,colors,linewidth=2):
    global cbar
    #This function plots a colorful line into axes "ax"
    #using data x and y and the colors from "colors"
    from matplotlib.collections import LineCollection
    points = np.array([x,y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments,cmap=plt.get_cmap('jet'),norm=plt.Normalize(np.min(colors), np.max(colors)))
    lc.set_array(colors[:-1])
    lc.set_linewidth(linewidth)
    ax.add_collection(lc)
    ax.plot(x,y,alpha=0)
    
    if cbar=='':
        cbar = fig.colorbar(lc,ax=ax,fraction=.1,shrink=0.9,pad=0.03)
    else:
        cbar.set_clim(np.min(colors),np.max(colors))
        cbar.draw_all()
    

q = -1.602e-19       #Coulombs   Charge of electron
h_bar = 1.054e-34   #J*s        Plank's Constant div by 2Pi
c = 3.0e8           #m/s        Speed of light
eo = 8.8541e-12     #C^2/(Nm^2) Permittivity of vacuum
me = 9.109e-31      #kg         Mass of electron

def EVD(t,wavelength,intensity=5e14,tb=0,phase=0,scale=1):
    #returns the position and velocity of the electron at a time or range of times
    #wavelength is in meters
    #I is in W/cm^2
    #tb is in seconds
    #phase is used to switch between a sin field (phase=0) or a cos field (phase = pi/2)
    
    Eo = np.sqrt(2*intensity*10**4/(3e8*8.85e-12))
    w = c/wavelength * 2 * np.pi
    p = phase
    
    E_field  = scale * Eo*np.sin(w*t+p)
    E_birth  = scale * Eo*np.sin(w*tb+p)
    distance = scale * q*Eo/(me*w) * ( np.cos(w*tb+p)*(t-tb) - np.sin(w*t+p)/w + np.sin(w*tb+p)/w ) 
    velocity = scale * q*Eo/(me*w) * ( np.cos(w*tb+p) - np.cos(w*t+p)) 
    
    vector_potential = scale * q*Eo/w * np.cos(w*t+p)
    
    return E_birth, E_field,velocity,distance

###################### Setting up the Figure #####################
fig = plt.figure(figsize=(12,10))
ax1 = plt.subplot(221)
ax2 = plt.subplot(223)

###Slider axes###
sl = 0.60
sw = 0.35
sh = 0.02
sv = 0.97 #start vertical
ss = -0.045 #spacing

axL1    = plt.axes([sl, sv+1*ss, sw, sh])
axL2    = plt.axes([sl, sv+2*ss, sw, sh])

axI1    = plt.axes([sl, sv+4*ss, sw, sh])
axI2    = plt.axes([sl, sv+5*ss, sw, sh])

axtbnum = plt.axes([sl, sv+7*ss, sw, sh])
axtbmin = plt.axes([sl, sv+8*ss, sw, sh])
axtbmax = plt.axes([sl, sv+9*ss, sw, sh])

axv0long = plt.axes([sl, sv+11*ss, sw, sh])
axv0tran = plt.axes([sl, sv+12*ss, sw, sh])

axtnum = plt.axes([sl, sv+14*ss, sw, sh])
axtsta = plt.axes([sl, sv+15*ss, sw, sh])
axtend = plt.axes([sl, sv+16*ss, sw, sh])

axcolor = plt.axes([sl, sv+18*ss, sw/2.5, sh*4])
axrot   = plt.axes([sl+0.2, sv+18*ss, sw/2.5, sh*4])


###Sliders###
sL1 = Slider(axL1, 'Wavelength 1 \n(nm)', 0,2000, valinit=400)
sL2 = Slider(axL2, 'Wavelength 2 \n(nm)', 0,2000, valinit=800)

sI1 = Slider(axI1, 'Intensity 1 \n(1e14 W/cm2)', 0,10, valinit=5)
sI2 = Slider(axI2, 'Intensity 2 \n(1e14 W/cm2)', 0,10, valinit=5)

stbnum = Slider(axtbnum, 'num tb',      0,100,  valinit=10)
stbmin = Slider(axtbmin, 'tb min\n(fs)', 0,2.66, valinit=.01)
stbmax = Slider(axtbmax, 'tb max\n(fs)', 0,2.66, valinit=2.66)

sv0long = Slider(axv0long, 'v0 tunnel dir\n(1e6 m/s)', -10,10, valinit=0)
sv0tran = Slider(axv0tran, 'v0 transverse dir\n(1e6 m/s)',-10,10, valinit=0)

stnum = Slider(axtnum, 'Num time steps', 0,1000, valinit=200)
ststa = Slider(axtsta, 'Start time\n(fs)', 0,10, valinit=.0001)
stend = Slider(axtend, 'Time duration\n(fs)', 0,10, valinit=6)

colorbutton = RadioButtons(axcolor, ('Time color', 'Energy color'))
global color_status; color_status = 'Time color'

rotbutton = RadioButtons(axrot, ('Counter rotating', 'Co rotating'))
global rotation; rotation = -1

plt.subplots_adjust(left=0.08,bottom=0.05,right=0.95,top=0.96,hspace=0.14)
###################### Done setting up the figure #####################

global cbar ;cbar=''


###################### This activates everytime something is clicked ###############
def update(val):
    global cbar
    print 'Updating with val = '+str(val)
    
    L1 = sL1.val*1e-9
    L2 = sL2.val*1e-9

    I1 = sI1.val*1e14
    I2 = sI2.val*1e14
    
    tbnum = stbnum.val
    if tbnum<1: tbnum=1
    tbmin = stbmin.val*1e-15
    tbmax = stbmax.val*1e-15
    
    v0long = sv0long.val*1e6
    v0tran = sv0tran.val*1e6
    
    tnum = stnum.val 
    tsta = ststa.val * 1e-15
    tend = stend.val * 1e-15
        
    if not val == 'init':
        xlim1 = ax1.get_xlim()
        ylim1 = ax1.get_ylim()
        xlim2 = ax2.get_xlim()
        ylim2 = ax2.get_ylim()
    
    ax1.clear()
    ax2.clear()
    for ax in (ax1,ax2): ax.axvline(0,alpha=0.2);ax.axhline(0,alpha=0.2)
    
    ax1.set_xlabel('Ex (V/m)')
    ax1.set_ylabel('Ey (V/m)')
    ax2.set_xlabel('x-direction (meters)')
    ax2.set_ylabel('y-direction (meters)')
    
###################### THIS IS WHERE THE MAGIC HAPPENS #####################
    
    for tb in np.linspace(tbmin,tbmax,tbnum):
        times = np.linspace(tb,tb+tend,tnum)
    
        Bx1,Ex1,vx1,dx1 = EVD(times,L1,intensity=I1,phase=np.pi*0 ,tb=tb)
        Bx2,Ex2,vx2,dx2 = EVD(times,L2,intensity=I2,phase=rotation*np.pi*0,tb=tb)

        By1,Ey1,vy1,dy1 = EVD(times,L1,intensity=I1,phase=np.pi*0.5 ,tb=tb)
        By2,Ey2,vy2,dy2 = EVD(times,L2,intensity=I2,phase=rotation*np.pi*0.5,tb=tb)

        Bx=Bx1+Bx2; By=By1+By2; 
        B = np.sqrt(Bx**2+By**2)
        
        v0x = Bx/B * v0long + By/B * v0tran
        v0y = By/B * v0long - Bx/B * v0tran

        Ex=Ex1+Ex2; vx=vx1+vx2+v0x; dx=dx1+dx2+v0x*(times-tb)
        Ey=Ey1+Ey2; vy=vy1+vy2+v0y; dy=dy1+dy2+v0y*(times-tb)

        ax1.plot(Bx,By,'bo',mec=None)
        sc = 1.0e5 #This is an arbitrary scale to make ax1 first figure look nice
        ax1.plot((Bx,Bx+v0x*sc),(By,By+v0y*sc),color='b')
        
        if color_status=='Energy color':
            energy = 0.5*me*(vx**2+vy**2)/1.602e-19
            plot_color_line(ax2,dx[times>tb],dy[times>tb],colors=energy)
        else:
            plot_color_line(ax2,dx[times>tb],dy[times>tb],colors=(times-tb)[times>tb])
        
    
    ax1.plot(Ex,Ey,color='r')
####################### NO MORE MAGIC BELOW HERE ############################################            
    
    if not val == 'init':
        ax1.set_xlim(*xlim1)  
        ax1.set_ylim(*ylim1) 
        ax2.set_xlim(*xlim2)  
        ax2.set_ylim(*ylim2)
     
    fig.canvas.draw_idle()

def color_update(val):
    global color_status
    color_status = val
    update('Updating color')

def rotation_update(val):
    global rotation
    if val == 'Counter rotating': rotation = -1
    else: rotation = 1
    update('Updating rotation')
    

for s in (sL1,sL2,stbmin,stbmax,stbnum,sI1,sI2,sv0long,sv0tran,stnum,ststa,stend):
    s.on_changed(update)

colorbutton.on_clicked(color_update)
rotbutton.on_clicked(rotation_update)
    
update('init') #make the plots of the first time

plt.show()  