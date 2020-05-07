import random
import math
import numpy as np
kB=1.38#Boltzmann Constant
T=2#Temperature
m=5#mass of a atom
l=7#length of the box
rc=2.5#cutoff distance
x=np.zeros(512)#r(t)
y=np.zeros(512)
z=np.zeros(512)
px=np.zeros(512)#r(t-dt)
py=np.zeros(512)
pz=np.zeros(512)
nx=np.zeros(512)#r(t+dt)
ny=np.zeros(512)
nz=np.zeros(512)
u=np.zeros(512)#potential energy
ke=np.zeros(512)#kinetic energy
e=np.zeros(512)#total energy
c=0
#Initial Coordinate Assignment
for i in range(0,7):
    for j in range(0,7):
        for k in range(0,7):
            x[c]=i
            y[c]=j
            z[c]=k
            c=c+1
vx=np.zeros(512)#v(t)
vy=np.zeros(512)
vz=np.zeros(512)
#Initial Velocity Assignment
for i in range(0,511):
    vx[i]=random.gauss(0,1)*math.sqrt(kB*T/m)
    vy[i]=random.gauss(0,1)*math.sqrt(kB*T/m)
    vz[i]=random.gauss(0,1)*math.sqrt(kB*T/m)
dt=0.001 #timestep
#Declaration of Force field
fx=np.zeros(512)
fy=np.zeros(512)
fz=np.zeros(512)
#Simulation Begins for 100 time steps
for i in range(1,10):
	#Lennard Jones Force Field Calculation
	for i in range(0,511):
        	for j in range(0,511):
            		if(i!=j):
				dx=x[j]-x[i]
				dy=y[j]-y[i]
				dz=z[j]-z[i]
				# Using Minimum Image Convention
				if(dx>l/2):
				    dx=dx-l
				elif(dx<-l/2):
				    dx=dx+l
				if(dy>l/2):
				    dy=dy-l
				elif(dy<-l/2):
				    dy=dy+l
				if(dz>l/2):
				    dz=dz-l
				elif(dz<-l/2):
				    dz=dz+l    
				r=math.sqrt(dx*dx+dy*dy+dz*dz)
				if(r>rc):
				    fx[i]=fx[i]+(48/(r*r))*(((1/r)**12)-(0.5*((1/r)**6)))*dx
				    fy[i]=fy[i]+(48/(r*r))*(((1/r)**12)-(0.5*((1/r)**6)))*dy
				    fz[i]=fz[i]+(48/(r*r))*(((1/r)**12)-(0.5*((1/r)**6)))*dz
    	#Calculating r(t+dt) using Verlet method
    	for i in range(0,511):
		nx[i]=(2*x[i])-px[i]+(fx[i]*((dt**2)/m))
		ny[i]=(2*y[i])-py[i]+(fy[i]*((dt**2)/m))
		nz[i]=(2*z[i])-pz[i]+(fz[i]*((dt**2)/m))
		#Using Periodic Boundary Conditions
		if(nx[i]>l):
		    nx[i]=nx[i]-l
		elif(nx[i]<0):
		    nx[i]=nx[i]+l
		if(ny[i]>l):
		    ny[i]=ny[i]-l
		elif(ny[i]<0):
		    ny[i]=ny[i]+l
		if(nz[i]>l):
		    nz[i]=nz[i]-l
		elif(nz[i]<0):
		    nz[i]=nz[i]+l    
		vx[i]=(nx[i]-px[i])/(2*dt)
		vy[i]=(ny[i]-py[i])/(2*dt)
		vz[i]=(nz[i]-pz[i])/(2*dt)
    	#Changing coordinates
    	for i in range(0,511):
		px[i]=x[i]
		py[i]=y[i]
		pz[i]=z[i]
		x[i]=nx[i]
		y[i]=ny[i]
		z[i]=nz[i]
		#Calculating Potential Energy
		for j in range(0,511):
			if(i!=j):
				dx=x[j]-x[i]
				dy=y[j]-y[i]
				dz=z[j]-z[i]
				#Reusing Minimum Image Convention
				if(dx>l/2):
				    dx=dx-l
				elif(dx<-l/2):
				    dx=dx+l
				if(dy>l/2):
				    dy=dy-l
				elif(dy<-l/2):
				    dy=dy+l
				if(dz>l/2):
				    dz=dz-l
				elif(dz<-l/2):
				    dz=dz+l    
				r=math.sqrt(dx*dx+dy*dy+dz*dz)
				if(r>rc):
				    u[i]=4*(((1/r)**12)-((1/r)**6))    
		#Calculating Kinetic Energy
		ke[i]=(0.5*m*((vx[i]*vx[i])+(vy[i]*vy[i])+(vz[i]*vz[i])))
		#Calculating Total Energy
		e[i]=u[i]+ke[i]
    e[i]=e[i]/512
    print e[i]

#plot Avg PE vs MD Step
#plot Avg KE vs MD Step
#plot Avg Eng vs MD Step
    
    

                
            


    
            
