#!usr/bin/python
"""
*********************************************
***   Program for the computation of the  ***
***          radiative inhibition         ***         
*********************************************
This code is part of the LIFELINE program.

Copyright (C) 2020-2021 University of Liege (Belgium)
Enmanuelle Mossoux (STAR Institute)
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.
You should have received a copy of the GNU General Public License along with this program.
If not, see <http://www.gnu.org/licenses/>.

Run: from radiative_inhibition.py import radiative_inhibition
     radiative_inhibition_comp(Mdot1, Mdot2, Per, mass_ratio, R1, R2, Teff_1, Teff_2, d, direct, mu, H_mass_frac, beta1, beta2, ipar)

# Parameters 
# ==========
Mdoti - The mass-loss rate [Msun/yr]
Per - Orbital period [days]
mass_ratio - M1/M2
Ri - Stellar radius [Rsun]
Teff_i - Effective temperature [K]
d - Separation [cm]
direct - Directory where to save the results
mu - Mean molecular weight. Default value: 0.62 (totally ionized gas)
H_mass_frac - Fractional mass abundance of H. Default value: 0.7381 (Sun)
betai - Parameter of the beta law of the wind velocity for each star. Default value: 1.
ipar - Index of the binary system

# Versions
# ========
v1 - 1/3/2019 - Application of radiative inhibition_models.py for a particular binary.
"""      

# Subroutine: to plot stellar surface
# ===================================
def plot_star(centerx,centery,radius,x):
	import numpy as np
	return np.sqrt(radius**2-(x-centerx)**2)+centery

# Subroutine: velocity at the critical point
# ==========================================
def vz(M1,Gamma1,R1,M2,Gamma2,R2,grav,sound_vel,d,uc,hc,ud,alpha1,alpha2):
	import math
	Ac=(ud/(uc-ud))**2
	Cz=2.*grav*M1*(1.-Gamma1)/sound_vel**2
	factor=R1*uc/Cz
	K1=(1.-(1.-factor**2)**(1.+alpha1))/((1.+alpha1)*factor**2)
	dK1du=2.*R1*((1.-factor**2)**alpha1*(1.+alpha1)*factor**2-1.+(1-factor**2)**(1.+alpha1))/(Cz*(1.+alpha1)*factor**3)
	factor=(R2/(d+Cz/uc))**2
	K2=(1.-(1-factor)**(1.+alpha2))/((1.+alpha2)*factor)
	dK2du=(2.*Cz/((d+Cz/uc)*uc**2))*((1.-factor)**alpha2*(1.+alpha2)*factor-1.+(1-factor)**(1.+alpha2))/((1.+alpha2)*factor)
	dAdu=-2.*ud/(uc-ud)**3
	Bc=abs(K1-M2*Gamma2/(M1*Gamma1)*K2*Ac)
	dBdu=dK1du-M2*Gamma2*(dK2du*Ac+dAdu*K2)/(M1*Gamma1)
	dhdu=4.*M2*(1.-Gamma2)*(dAdu-Ac/uc)/(uc*M1*(1.-Gamma1))
	omegac=1.+alpha1*hc/((1.-alpha1)*(hc*dBdu/((1.-alpha1)*Bc)-dhdu)**0.5)
	return sound_vel*math.sqrt(omegac)

# Subroutine: change of coordinate
# ================================
def uz(z,M,Gamma,grav,sound_vel):
	return -2.*grav*M*(1.-Gamma)/(z*sound_vel**2)

# Subroutines: compute velocity
# =============================
def dv_ligne_centre(x, urr, r, dMdotdsigma1, sound_vel, deltar, alpha1, rad, gravity):
	import math
	rho1=dMdotdsigma1/(abs(x)*4.*math.pi*r**2) #[g/cm^3]
	eq=(1.-sound_vel**2/x**2)*urr*(x-urr)/deltar-2.*sound_vel**2/r+gravity-rad*(abs(x-urr)/(rho1*deltar))**alpha1
	return eq**2

def dv_hors_ligne_centre(p, uthetar,urr, uthetat, urt, r, dMdotdsigma1, theta, theta2, sound_vel, deltar, deltatheta, alpha1, rad1, rad2, grav2, gravity, ir):
	import numpy as np
	import math

    	x, y = p
	vhere=np.sqrt(x**2+y**2)
	rho1=dMdotdsigma1/(vhere*4.*math.pi*r**2) #[g/cm^3]
	sound=sound_vel**2/vhere**2

	angler=math.cos(theta)*math.cos(theta2)-math.sin(theta)*math.sin(theta2)
	eq1=(1.-sound)*x*(x-urr)/deltar - 2.*sound_vel**2/r + y*(x-urt)/(r*deltatheta) - sound*y*(y-uthetar)/deltar - y**2/r - grav2*angler + gravity - (rad1-rad2*angler)*abs((x*(x-urr)+y*(y-uthetar))/(vhere*rho1*deltar))**alpha1

	angler=math.cos(theta)*math.sin(theta2)+math.sin(theta)*math.cos(theta2)
	eq2=(1.-sound)*y*(y-uthetat)/(r*deltatheta) - sound*x*(x-urt)/(r*deltatheta) + y*x/r + x*(y-uthetar)/deltar - grav2*angler +rad2*angler*abs((x*(x-urr)+y*(y-uthetar))/(vhere*rho1*deltar))**alpha1
	if (abs(x*(x-urr)+y*(y-uthetar)) == 0.):
		eq1=(1.-sound)*x*(x-urr)/deltar - 2.*sound_vel**2/r + y*(x-urt)/(r*deltatheta) - sound*y*(y-uthetar)/deltar - y**2/r - grav2*angler + gravity
		eq2=(1.-sound)*y*(y-uthetat)/(r*deltatheta) - sound*x*(x-urt)/(r*deltatheta) + y*x/r + x*(y-uthetar)/deltar - grav2*angler
	return math.sqrt(eq1**2+eq2**2)

def jac(p, uthetar,urr, uthetat, urt, r, dMdotdsigma1, theta, theta2, sound_vel, deltar, deltatheta, alpha1, rad1, rad2, grav2, gravity, ir):
	import numpy as np
	import math

	x,y=p
	vhere=np.sqrt(x**2+y**2)
	rho1=dMdotdsigma1/(vhere*4.*math.pi*r**2) #[cm^-3]
	sound=sound_vel**2/vhere**2
	const2=abs(x*(x-urr)+y*(y-uthetar))/vhere

	angler=math.cos(theta)*math.cos(theta2)-math.sin(theta)*math.sin(theta2)
	const1=alpha1*(rad1-rad2*angler)*abs((x*(x-urr)+y*(y-uthetar))/(vhere*rho1*deltar))**(alpha1-1)/(vhere**2*rho1*deltar)
	if (abs(x*(x-urr)+y*(y-uthetar)) == 0.):
		const1=0.
		const2=0.
	eq1x=((y**2-sound_vel**2)*(2.*x-urr)+x**2*(4.*x-3.*urr))/(deltar*vhere**2)+2.*x*(sound_vel**2*y*(y-uthetar)-x*(x-urr)*(vhere**2-sound_vel**2))/(deltar*vhere**4)+y/(r*deltatheta)-const1*(2.*x-urr)
	eq1y=(-sound_vel**2*(y*2.-uthetar)+2.*y*x*(x-urr))/(deltar*vhere**2)+2.*y*(sound_vel**2*y*(y-uthetar)-x*(x-urr)*(vhere**2-sound_vel**2))/(deltar*vhere**4)-2.*y/r+(x-urr)/(r*deltatheta)-const1*(2.*y-uthetar)
	eq1=(1.-sound)*x*(x-urr)/deltar - 2.*sound_vel**2/r + y*(x-urt)/(r*deltatheta) - sound*y*(y-uthetar)/deltar - y**2/r - grav2*angler + gravity - (rad1-rad2*angler)*abs((x*(x-urr)+y*(y-uthetar))/(vhere*rho1*deltar))**alpha1

	angler=math.cos(theta)*math.sin(theta2)+math.sin(theta)*math.cos(theta2)
	const1=alpha1*rad2*angler*abs((x*(x-urr)+y*(y-uthetar))/(vhere*rho1*deltar))**(alpha1-1)/(vhere**2*rho1*deltar)
	if (abs(x*(x-urr)+y*(y-uthetar)) == 0.):
		const1=0.
	eq2y=((x**2-sound_vel**2)*(2.*y-uthetat)+y**2*(4.*y-3.*uthetat))/(r*deltatheta*vhere**2)+2.*y*(sound_vel**2*x*(x-urt)-y*(y-uthetat)*(vhere**2-sound_vel**2))/(r*deltatheta*vhere**4)+x/r+x/deltar-const1*(2.*y-uthetar)
	eq2x=(-sound_vel**2*(x*2.-urt)+2.*y*x*(y-uthetat))/(r*deltatheta*vhere**2)+2.*x*(sound_vel**2*x*(x-urt)-y*(y-uthetat)*(vhere**2-sound_vel**2))/(r*deltatheta*vhere**4)+y/r+(x-uthetar)/deltar-const1*(2.*x-urr)
	eq2=(1.-sound)*y*(y-uthetat)/(r*deltatheta) - sound*x*(x-urt)/(r*deltatheta) + y*x/r + x*(y-uthetar)/deltar - grav2*angler +rad2*angler*abs((x*(x-urr)+y*(y-uthetar))/(vhere*rho1*deltar))**alpha1

	if (math.isnan(eq1x) or math.isnan(eq2x) or math.isnan(eq1y) or math.isnan(eq2y)):
		angler=math.cos(theta)*math.cos(theta2)-math.sin(theta)*math.sin(theta2)
		eq1x=((y**2-sound_vel**2)*(2.*x-urr)+x**2*(4.*x-3.*urr))/(deltar*vhere**2)+2.*x*(sound_vel**2*y*(y-uthetar)-x*(x-urr)*(vhere**2-sound_vel**2))/(deltar*vhere**4)+y/(r*deltatheta)
		eq1y=(-sound_vel**2*(y*2.-uthetar)+2.*y*x*(x-urr))/(deltar*vhere**2)+2.*y*(sound_vel**2*y*(y-uthetar)-x*(x-urr)*(vhere**2-sound_vel**2))/(deltar*vhere**4)-2.*y/r+(x-urr)/(r*deltatheta)
		eq1y=((x**2-sound_vel**2)*(2.*y-uthetat)+y**2*(4.*y-3.*uthetat))/(r*deltatheta*vhere**2)+2.*y*(sound_vel**2*x*(x-urt)-y*(y-uthetat)*(vhere**2-sound_vel**2))/(r*deltatheta*vhere**4)+x/r+x/deltar
		eq1x=(-sound_vel**2*(x*2.-urt)+2.*y*x*(y-uthetat))/(r*deltatheta*vhere**2)+2.*x*(sound_vel**2*x*(x-urt)-y*(y-uthetat)*(vhere**2-sound_vel**2))/(r*deltatheta*vhere**4)+y/r+(x-uthetar)/deltar
		eq1=(1.-sound)*x*(x-urr)/deltar - 2.*sound_vel**2/r + y*(x-urt)/(r*deltatheta) - sound*y*(y-uthetar)/deltar - y**2/r - grav2*angler + gravity
		eq2=(1.-sound)*y*(y-uthetat)/(r*deltatheta) - sound*x*(x-urt)/(r*deltatheta) + y*x/r + x*(y-uthetar)/deltar - grav2*angler
	return np.array([(2.*eq1*eq1x+2.*eq2*eq2x),(2.*eq1*eq1y+2.*eq2*eq2y)])

# Main
# ====
def radiative_inhibition_comp(type1, type2, Mdot1, Mdot2, Per, mass_ratio, R1, R2, Teff_1, Teff_2, d, direct, mu, H_mass_frac, beta1, beta2, ipar):
	import sys
	import math
	import numpy as np
	from constantes import constante
	import matplotlib.pyplot as plt
	import matplotlib.mlab as mlab
	import matplotlib
	from scipy.optimize import minimize,minimize_scalar
	import os
	import pandas as pd
	import time

	time1=time.time()

	Rsun=constante('Rsun') #[cm]
	Msun=constante('Msun') #[g]
	grav=constante('G')*1000. #[cm^3/g/s^2]
	light_vel=constante('c')*100. #[cm/s]
	sigma_t=constante('sigma_t')*10000. #[cm^2]
	sigma_k=constante('sigma_k') #[W/m^2/K^4]
	mh=(constante('mp')+constante('me'))*1000. #[g]
	mp=constante('mp')*1000. #[g]
	kb=constante('kb') #Boltzmann constant [erg/K]
	sigma_e=(1+H_mass_frac)*sigma_t/(2.*mh)	#scattering electron opacity [cm^2/g]       
	pi= math.pi

	abbott82_teff=np.array([6000., 6000., 6000., 8000., 8000., 8000., 10000., 10000., 10000., 15000., 15000., 15000., 20000., 20000., 20000., 30000., 30000., 30000., 40000., 40000., 40000., 50000., 50000., 50000.])
	abbott82_dens=np.array([6.8e6, 6.8e8, 6.8e11, 1.7e6, 1.7e9, 1.7e12, 3.2e6, 3.2e9, 3.9e12, 1.3e7, 1.3e10, 1.3e13, 3.0e7, 3.0e10, 3.0e13, 1.0e8, 1.0e11, 1.0e14, 1.8e8, 1.8e11, 1.8e14, 3.1e8, 3.1e11, 3.1e14])
	k_vec=np.array([0.018, 0.029, 0.272, 0.11, 0.105, 0.108, 0.288, 0.362, 0.37, 0.189, 0.253, 0.945, 0.14, 0.477, 0.617, 0.093, 0.156, 0.571, 0.051, 0.174, 0.533, 0.089, 0.178, 0.472])
	alpha_vec=[0.502, 0.465, 0.444, 0.521, 0.542, 0.555, 0.499, 0.538, 0.54, 0.505, 0.511, 0.517, 0.559, 0.506, 0.523, 0.576, 0.609, 0.545, 0.684, 0.606, 0.571, 0.64, 0.606, 0.582]
	Tuniq=np.unique(abbott82_teff)

	mass_sum=4.*pi**2*d**3/(grav*(Per*86400.)**2) #g
	M2=mass_sum/(1.+mass_ratio)
	M1=mass_sum-M2
	R1=R1*Rsun #cm
	Mdot1=Mdot1*Msun/(86400.*365.25) #g/s
	R2=R2*Rsun #cm
	Mdot2=Mdot2*Msun/(86400.*365.25) #g/s

	nbr_star=2-(type1 == type2)
	for istar in range(nbr_star):
		if (istar == 1):
			Mdot1, Mdot2 = Mdot2, Mdot1
			R1, R2 = R2, R1
			M1, M2 = M2, M1
			Teff_1, Teff_2 = Teff_2, Teff_1
			beta1, beta2 = beta2, beta1
			vinf1, vinf2 = vinf2, vinf1

		v1=2.6*math.sqrt(2.*grav*M1/R1) #vinf=2.6*vesc
		v2=2.6*math.sqrt(2.*grav*M2/R2)
		vinf1=v1
		vinf2=v2
		vlocal=[v1*(1.-0.99*R1/(d/2.))**beta1,v2*(1.-0.99*R2/(d/2.))**beta2] #cm/s
		if ((min(vlocal)/1.0E8)**4*((d/2.)/1.0e12)/(max([Mdot1,Mdot2])*365.24*86400./(1.0e-5*Msun))<1):
			v1=v1*(1.-0.99*R1/(d/2.))**beta1
			v2=v2*(1.-0.99*R2/(d/2.))**beta2

		L1=4.*sigma_k*pi*Teff_1**4*(R1/100.)**2 # [W = 10^7 g*cm*cm/s*s*s]
		L2=4.*sigma_k*pi*Teff_2**4*(R2/100.)**2 # [W]

		sound_vel=math.sqrt(kb*Teff_1/(mh*mu)) # isothermal flow [cm/s]
		Gamma1=sigma_t*L1*1e7/(M1*4*pi*grav*light_vel*mp) # Eddington ratio (L/Le)
		Gamma2=sigma_t*L2*1e7/(M2*4*pi*grav*light_vel*mp)
		rhomean=Mdot1/(mh*v1*4.*pi*(d/2.)**2)
		ind0=np.argmin(abs(Tuniq-Teff_1))
		ind=np.argmin(abs(abbott82_dens[int(ind0*3.):int(ind0*3.+3.)]-rhomean))
		alpha1=alpha_vec[int(ind0*3.+ind)]
		k1=k_vec[int(ind0*3.+ind)]
		rhomean=Mdot2/(mh*v2*4.*pi*(d/2.)**2)
		ind0=np.argmin(abs(Tuniq-Teff_2))
		ind=np.argmin(abs(abbott82_dens[int(ind0*3.):int(ind0*3.+3.)]-rhomean))
		alpha2=alpha_vec[int(ind0*3.+ind)]
		alpha2=np.mean(alpha_vec)
		K2c=(1.-(1-(R2/(d-1.05*R2))**2)**(1.+alpha2))/((1.+alpha2)*(R2/(d-1.05*R2))**2)
		if (R2/(d-1.05*R2)>1):
			K2c=(1.-(1-(R2/(1.05*R2))**2)**(1.+alpha2))/((1.+alpha2)*(R2/(1.05*R2))**2)
		K1c=(1.-(1-(1./1.05)**2)**(1.+alpha1))/((1.+alpha1)*(1./1.05)**2)
		c1=M1*(1.-Gamma1)/(M2*(1.-Gamma2))
		c2=L1*K1c/(L2*K2c)
		ud=uz(d,M1,Gamma1,grav,sound_vel)
		uc=uz(1.05*R1,M1,Gamma1,grav,sound_vel) #Where the opacity is equal to 1
		hc=1.+4.*M2*(1.-Gamma2)*(ud/(uc-ud))**2/(uc*M1*(1.-Gamma1))

		# Equation 24 of Stevens et al. 1994
		dMdotdsigma1=Mdot1*(1.-alpha1*(ud/(uc-ud))**2*(c2-(1.-alpha1)*c1))

		# Velocity at the critical point
		Cz=2.*grav*M1*(1.-Gamma1)/sound_vel**2 #cm
		zc=-Cz/uc
		vc=vz(M1,Gamma1,R1,M2,Gamma2,R2,grav,sound_vel,d,uc,hc,ud,alpha1,alpha2)

		Vth1=math.sqrt(2.*kb*Teff_1/mh) #[cm/s]

		# along the line of center (theta=0)
		urr=vc
		r=zc
		deltar=d/1000.
		nur=(2.*d-R2-zc)/deltar
		ur_vec=[urr]
		utheta_vec=[0.]
		rvec=[r]
		x_vec=[r]
		y_vec=[0.]
		ux=[urr]
		uy=[0.]

		for ir in range(int(nur-1)):
			r=r+deltar
			r2=d-r
			K1=(1.-(1-(R1/r)**2)**(1.+alpha1))/((1.+alpha1)*(R1/r)**2)
			flux=L1/(4.*math.pi*r**2) # [W/cm^2]
			flux2=L2/(4.*pi*r2**2) # [W/cm^2]
			K2=0. #Inside star 2
			if (abs(r2) > R2):
				K2=(1.-(1.-(R2/r2)**2)**(1.+alpha2))/((1.+alpha2)*(R2/r2)**2)
			rad1=1.0e7*sigma_e**(1.-alpha1)*k1*(flux*K1-flux2*K2)/(light_vel*Vth1**alpha1)
			gravity=-grav*M2*(1.-Gamma2)/r2**2+grav*M1*(1.-Gamma1)/r**2
			sol1 = minimize_scalar(dv_ligne_centre, bounds=(urr,v1), method='bounded',args=(urr, r, dMdotdsigma1, sound_vel, deltar, alpha1, rad1, gravity))
			urr=sol1.x
			ur_vec.append(urr)
			ux.append(urr)
			utheta_vec.append(0.)
			uy.append(0.)
			rvec.append(r)
			x_vec.append(r)
			y_vec.append(0.)
		nur=len(rvec)
		rvec=np.array(rvec)

		# Outside the line of center (theta>0)
		theta=0.
		deltatheta= pi/179.
		nutheta= pi/deltatheta

		r=rvec[1:]
		K1=(1.-(1-(R1/r)**2)**(1.+alpha1))/((1.+alpha1)*(R1/r)**2)
		flux=L1/(4.*pi*r**2) # [W/cm^2]
		rad1=1.0e7*sigma_e**(1.-alpha1)*k1*flux*K1/(light_vel*Vth1**alpha1)
		gravity=grav*M1*(1.-Gamma1)/r**2 #[cm/s^2]
		for _ in range(int(nutheta)):
			urr=ur_vec[0]
			uthetar=0.
			utheta=utheta_vec[1]
			ur=ur_vec[1]
			theta=theta+deltatheta
			print("Theta: "+str(round(theta*180./pi))+" on 180 degree")
			ux.append(urr*math.cos(theta)-uthetar*np.sin(theta))
			uy.append(urr*math.sin(theta)+uthetar*np.cos(theta))
			r2=np.sqrt(d**2+r**2-2.*r*d*math.cos(theta))
			ind=np.where((abs(r2) < R2))[0]
			K2=(1.-(1-(R2/r2)**2)**(1.+alpha2))/((1.+alpha2)*(R2/r2)**2)
			flux2=L2/(4.*pi*r2**2) # [W/cm^2]
			if (len(ind) > 0):
				K2[ind]=0.
			rad2=1.0e7*sigma_e**(1.-alpha1)*k1*flux2*K2/(light_vel*Vth1**alpha1) #[cm/s^2]
			grav2=grav*M2*(1.-Gamma2)/r2**2
			theta2=np.arccos((d-r*math.cos(theta))/r2)
			tole=0.1

			try:
				for ir in range(int(nur-1)):
					# Occultation by star 1	
					if (theta > pi/2. and theta2[ir] < math.atan(R1/d)):
						rad2=r2*0.
					cons = ({'type': 'ineq', 'fun': lambda x:  x[0]**2 + x[1]**2 - urr**2 - uthetar**2})
					if (utheta == 0):
						bnds = ((ur/1.1, ur*1.1), (0., 1.e7))
					else:
						bnds = ((ur/1.1, ur*1.1), (utheta/1.1, utheta*1.1))
					sol = minimize(dv_hors_ligne_centre, [ur*1.05,utheta*1.05], bounds=bnds,args=(uthetar,urr, utheta, ur, r[ir], dMdotdsigma1, theta, theta2[ir], sound_vel, deltar, deltatheta, alpha1, rad1[ir], rad2[ir], grav2[ir], gravity[ir], ir),jac=jac,constraints=cons, tol=tole,options={'maxiter':100})
					urr=ur
					uthetar=utheta
					if (sol.success == True):
						tole=sol.fun/50.
						uthetar=sol['x'][1]
						urr=sol['x'][0]

					ur_vec[ir+1]=urr
					utheta_vec[ir+1]=uthetar
					if (ir < int(nur-2)):
						ur=ur_vec[ir+2]
						utheta=utheta_vec[ir+2]
					ux.append(urr*math.cos(theta)-uthetar*np.sin(theta))
					uy.append(urr*math.sin(theta)+uthetar*np.cos(theta))
					if (ux[-1]**2+uy[-1]**2 > vinf1**2):
						if ux[-1] < 0:
							angle=pi-math.atan(-uy[-1]/ux[-1])
						else:
							angle=math.atan(uy[-1]/ux[-1])
						ux[-1]=vinf1*math.cos(angle)
						uy[-1]=vinf1*math.sin(angle)

				x_vec=np.concatenate((x_vec,rvec*math.cos(theta)))
				y_vec=np.concatenate((y_vec,rvec*math.sin(theta)))
			except KeyboardInterrupt:
				ux=ux[:len(x_vec)]
				uy=uy[:len(x_vec)]
				break

		ux=np.array(ux)
		uy=np.array(uy)
		x_vec=np.array(x_vec)
		y_vec=np.array(y_vec)

		# Occultation by star 2
		ux[((np.sqrt(x_vec*x_vec+y_vec*y_vec) > d-R2) & (y_vec/x_vec < R2/d)) & (x_vec > 0.)]=0.
		uy[((np.sqrt(x_vec*x_vec+y_vec*y_vec) > d-R2) & (y_vec/x_vec < R2/d)) & (x_vec > 0.)]=0.

		############
		### Save ###
		############
		os.chdir(direct)
		os.remove("wind_star"+str(istar+1)+"_param"+str(ipar)+".h5") if os.path.exists("wind_star"+str(istar+1)+"_param"+str(ipar)+".h5") else None
		df_xyz = pd.DataFrame(dict(x=np.array(x_vec)/d, y=np.array(y_vec)/d, ux=ux, uy=uy))
		df_xyz.to_hdf(direct+"/winds/wind_star"+str(istar+1)+"_param"+str(ipar)+".h5",'wind',format='t')
		if (nbr_star==1):
			os.remove("wind_star2_param"+str(ipar)+".h5") if os.path.exists("wind_star2_param"+str(ipar)+".h5") else None
			df_xyz = pd.DataFrame(dict(x=np.array(x_vec)/d, y=np.array(y_vec)/d, ux=ux, uy=uy))
			df_xyz.to_hdf(direct+"/winds/wind_star2_param"+str(ipar)+".h5",'wind',format='t')
		del df_xyz

		############
		### Plot ###
		############
		from plots_winds import plots
		nothing=plots(R1,R2,vinf1,d,nbr_star,istar,ipar,direct+'/winds')

	time3=time.time()
	print("I have worked during "+str((time3-time1)/60.)+" min or "+str((time3-time1)/3600.)+" hours.")

