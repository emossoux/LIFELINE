#!usr/bin/python
"""
***************************************************        
***   Program for the computation of the shock  ***
***     in radiative colliding wind binaries    ***         
***************************************************
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

Run: from rshock import rshock
     rien=rshock(Mdot1, Mdot2, vinf1, vinf2, mass_ratio, a, per, R1, R2, beta1, beta2, Teff_1, Teff_2, nbr_points_shock_2D, nbr_points_shock_3D, nstep_width, nbr_bin_profile, direct, cut_lim, wind_prim, wind_sec, ex, omega, phase, incl, ipar, mu)

Use the Antokhin et al. (2004) formalism for the shock creation

# Parameters 
# ==========
Mdoti - Mass loss rate [Msun/yr]
vinfi - terminal velocity of the winds [km/s]
mass_ratio - mass ratio
a - semi-major axis [Rsun]
per  - orbital period [days]
Ri - stellar radius [Rsun]
betai - index of the beta velocity law
Teff_i - efective temperature [K]
nbr_points_shock_2D - number of points to create the 2D shape of the shock
nbr_points_shock_3D - number of points to create the 3D shape of the shock
nstep_width - number of points to discretize the inside of the shock
direct - directory where to save files
cut_lim - Limit of the tangential velocity to see the emission
wind_prim - h5 file containing the wind distribution of the primary 
wind_sec - h5 file containing the wind distribution of the secondary 
ex - eccentricity
omega - argument of the periapsis [radian]
phase - phase of the system
incl - inclination of the orbit
ipar - index on the set of stellar parameter
mu - Mean molecular weight

# Versions
# ========
v1 - 16/03/18
v2 - 26/03/19 - Taking the radiative inhibition into account
"""      
# Subroutine: compute eta
# =======================
def comp_eta(x, y, R1, R2, Mdot1, Mdot2, d, tree_prim_wind, tree_sec_wind, dwind_prim, dwind_sec):
	import numpy as np
	if (x**2+y**2 < R1**2):
		vv1=100.e5
		# closest point
		rien, ind_close = tree_sec_wind.query([[1.-R1/d,y/d]])
		ind_close=ind_close[0]
		vv2=np.sqrt(dwind_sec.iloc[ind_close]['ux']**2+dwind_sec.iloc[ind_close]['uy']**2)
		if (dwind_sec.iloc[ind_close]['ux'] == 0 and dwind_sec.iloc[ind_close]['uy']==0):
			ind_closep=ind_close+1
			ind_closem=ind_close-1
			stade="up"
			dw1=dwind_sec.iloc[ind_closep]['ux']
			dw2=dwind_sec.iloc[ind_closep]['uy']
			ind_closep=ind_close+1
			while (dw1 == 0 and dw2==0):
				if (stade=="up"):
					dw1=dwind_sec.iloc[ind_closem]['ux']
					dw2=dwind_sec.iloc[ind_closem]['uy']
					ind_closem=ind_closem-1
					ind_close=ind_closem
					stade="down"
				else:
					stade="up"
					dw1=dwind_sec.iloc[ind_closep]['ux']
					dw2=dwind_sec.iloc[ind_closep]['uy']
					ind_closep=ind_closep+1
					ind_close=ind_closep
			vv2=np.sqrt(dwind_sec.iloc[ind_close]['ux']**2+dwind_sec.iloc[ind_close]['uy']**2)
	if ((d-x)**2+y**2 < R2**2):
		vv2=100.e5
		# closest point
		rien, ind_close = tree_prim_wind.query([[R2/d,y/d]])
		ind_close=ind_close[0]
		vv1=np.sqrt(dwind_prim.iloc[ind_close]['ux']**2+dwind_prim.iloc[ind_close]['uy']**2)
		if (dwind_prim.iloc[ind_close]['ux'] == 0 and dwind_prim.iloc[ind_close]['uy']==0):
			ind_closep=ind_close+1
			ind_closem=ind_close-1
			stade="up"
			dw1=dwind_prim.iloc[ind_closep]['ux']
			dw2=dwind_prim.iloc[ind_closep]['uy']
			ind_closep=ind_close+1
			while (dw1 == 0 and dw2==0):
				if (stade=="up"):
					dw1=dwind_prim.iloc[ind_closem]['ux']
					dw2=dwind_prim.iloc[ind_closem]['uy']
					ind_closem=ind_closem-1
					ind_close=ind_closem
					stade="down"
				else:
					stade="up"
					dw1=dwind_prim.iloc[ind_closep]['ux']
					dw2=dwind_prim.iloc[ind_closep]['uy']
					ind_closep=ind_closep+1
					ind_close=ind_closep
			vv1=np.sqrt(dwind_prim.iloc[ind_close]['ux']**2+dwind_prim.iloc[ind_close]['uy']**2)
	if ((x**2+y**2 >= R1**2) and ((d-x)**2+y**2 >= R2**2)):
		# closest point
		rien, ind_close = tree_prim_wind.query([[x/d,y/d]])
		ind_close=ind_close[0]
		vv1=np.sqrt(dwind_prim.iloc[ind_close]['ux']**2+dwind_prim.iloc[ind_close]['uy']**2)
		if (dwind_prim.iloc[ind_close]['ux'] == 0 and dwind_prim.iloc[ind_close]['uy']==0):
			ind_closep=ind_close+1
			ind_closem=ind_close-1
			stade="up"
			dw1=dwind_prim.iloc[ind_closep]['ux']
			dw2=dwind_prim.iloc[ind_closep]['uy']
			ind_closep=ind_close+1
			while (dw1 == 0 and dw2==0):
				if (stade=="up"):
					dw1=dwind_prim.iloc[ind_closem]['ux']
					dw2=dwind_prim.iloc[ind_closem]['uy']
					ind_closem=ind_closem-1
					ind_close=ind_closem
					stade="down"
				else:
					stade="up"
					dw1=dwind_prim.iloc[ind_closep]['ux']
					dw2=dwind_prim.iloc[ind_closep]['uy']
					ind_closep=ind_closep+1
					ind_close=ind_closep
			vv1=np.sqrt(dwind_prim.iloc[ind_close]['ux']**2+dwind_prim.iloc[ind_close]['uy']**2)
		rien, ind_close = tree_sec_wind.query([[1.-x/d,y/d]])
		ind_close=ind_close[0]
		vv2=np.sqrt(dwind_sec.iloc[ind_close]['ux']**2+dwind_sec.iloc[ind_close]['uy']**2)
		if (dwind_sec.iloc[ind_close]['ux'] == 0 and dwind_sec.iloc[ind_close]['uy']==0):
			ind_closep=ind_close+1
			ind_closem=ind_close-1
			stade="up"
			dw1=dwind_sec.iloc[ind_closep]['ux']
			dw2=dwind_sec.iloc[ind_closep]['uy']
			ind_closep=ind_close+1
			while (dw1 == 0 and dw2==0):
				if (stade=="up"):
					dw1=dwind_sec.iloc[ind_closem]['ux']
					dw2=dwind_sec.iloc[ind_closem]['uy']
					ind_closem=ind_closem-1
					ind_close=ind_closem
					stade="down"
				else:
					stade="up"
					dw1=dwind_sec.iloc[ind_closep]['ux']
					dw2=dwind_sec.iloc[ind_closep]['uy']
					ind_closep=ind_closep+1
					ind_close=ind_closep
			vv2=np.sqrt(dwind_sec.iloc[ind_close]['ux']**2+dwind_sec.iloc[ind_close]['uy']**2)
	return (Mdot2*vv2)/(Mdot1*vv1)

# Subroutine: find x0
# ===================
def find_x0(x, d, R1, R2, v1, v2, Mdot1, Mdot2, tree_prim_wind, tree_sec_wind, dwind_prim, dwind_sec):
	import math
	eta=comp_eta(x, 0., R1, R2, Mdot1, Mdot2, d, tree_prim_wind, tree_sec_wind, dwind_prim, dwind_sec)
	return x-d/(1.+math.sqrt(eta))

# Subroutine: find dzeta
# ======================
def ddzeta_dy(x, yy,ddelta,aangle,sslope,MMdot):
	import numpy as np
	import math
	try:
		ret=MMdot*math.sin(ddelta)*yy/(4.*math.pi*(yy/math.sin(aangle))**2*math.sin(sslope))
	except:
		ret=0.
	return ret

# Subroutine: mass column
# =======================
def sigmas(y_new,z_new,xstag, y_vec, anglesign, deltasign, slopesign, Mdot, vsign, etasign, rhosign, beta, cote_choc, dzeta_old, y_old, R1, R2, d, radius, crashing):
	from scipy.integrate import odeint
	import numpy as np
	import math
	# Surface density of the cooled material within the front interaction layer
	if (y_new > 0.2*xstag):
		dzeta=0.
		if (math.sin(anglesign)*math.sin(slopesign) != 0.):
			dzeta = odeint(ddzeta_dy, dzeta_old, [y_old,y_new], args=(deltasign,anglesign,slopesign,Mdot)) # Equation A7 of Antokhin 2004
		sigma=dzeta/(y_new*vsign*np.cos(deltasign))
		dzeta_old=-1.0
		y_old=-1.0
	else:
		if (crashing=='yes'):
			c1=0.
			c2=beta/((d-R1)**2*((d-R1/R2)-1.))
			z0=(4.*np.sqrt(1.-np.sqrt(1./etasign))/(np.sqrt(1./etasign)*(d-R1))+(d-R1)*np.sqrt(1./etasign)*(c1-c2)/(1.+np.sqrt(1./etasign)))/(6.-(d-R1)**2*np.sqrt(1./etasign)*(c1+c2*np.sqrt(1./etasign))/(1.+np.sqrt(1./etasign)))
			if (z0 <0):
				z0=0.
			# Equation A4 of Antokhin 2004
			# Caution: eta of antokhin = our 1/eta -> change x0 to d-x0, r1_vec to r2_vec
			sigma_0=rhosign*(d-radius)/(2.*(1.+z0*(radius-d))) # Equation A3 of Antokhin 2004
			sigma=sigma_0*(1.-(y_new+z_new)**2*(1.+2.*(d-radius)*z0)/(6.*(d-radius)**2))
			if (cote_choc == 2):
				sigma_0=rhosign*radius/(2.*(1.+z0*radius)) # Equation A2 of Antokhin 2004
				sigma=sigma_0*(1.-(y_new+z_new)**2*(1.+2.*radius*z0)/(6.*radius**2))
		else:
			c1=beta/((xstag)**2*((xstag/R1)-1.))
			c2=beta/((d-xstag)**2*((d-xstag/R2)-1.))
			z0=(4.*np.sqrt(1.-np.sqrt(1./etasign))/(np.sqrt(1./etasign)*(d-xstag))+(d-xstag)*np.sqrt(1./etasign)*(c1-c2)/(1.+np.sqrt(1./etasign)))/(6.-(d-xstag)**2*np.sqrt(1./etasign)*(c1+c2*np.sqrt(1./etasign))/(1.+np.sqrt(1./etasign)))
			if (z0 <0):
				z0=0.
			# Equation A4 of Antokhin 2004
			sigma_0=rhosign*xstag/(2.*(1.+z0*xstag)) # Equation A2 of Antokhin 2004
			sigma=sigma_0*(1.-(y_new+z_new)**2*(1.+2.*(d-xstag)*z0)/(6.*(d-xstag)**2))
			if (cote_choc == 2):
				sigma_0=rhosign*(d-xstag)/(2.*(1.+z0*(xstag-d))) # Equation A3 of Antokhin 2004
				sigma=sigma_0*(1.-(y_new+z_new)**2*(1.+2.*xstag*z0)/(6.*xstag**2))
		dzeta_old=sigma*(y_new*vsign*np.cos(deltasign))
		y_old=y_new
	return [sigma,dzeta_old,y_old]

# Subroutine: evolution T
# =======================
def integrand(kT, kT_tab, cool_tab, direct, kT_ion, ion_frac, sunabund):
	import numpy as np
	import math
	from sun_abund import abund
	import pyatomdb
	from constantes import constante
	import os

	os.environ["ATOMDB"] = 'https://hea-www.cfa.harvard.edu/AtomDB/'
	# write user_data (where filemap and APED)
	fuser = open(direct+"/userdata", 'w')
	fuser.write("filemap="+direct+"/filemap\n")
	fuser.write("atomdbroot=https://hea-www.cfa.harvard.edu/AtomDB/\n")
	fuser.close()
	# read user_data
	setting=pyatomdb.util.load_user_prefs(adbroot=direct)

	ind=np.where(kT_tab <= kT*8.61732e-8)[0]
	if (kT*8.61732e-8 < min(kT_tab)):
		ind=[0]
	if (kT*8.61732e-8 > max(kT_tab)):
		ind=[int(len(kT_tab))-1]
	cool_tab_here=cool_tab[ind[-1]]/1.e-23

	nz=25
	sum_vec=0.
	number_fraction=abund(ref=sunabund)
	number_fraction=np.array(number_fraction)
	mass_fraction=[]
	for zz in range(nz+1):
		Z=zz+1
		mass_fraction.append(number_fraction[zz]*pyatomdb.atomic.Z_to_mass(Z)/6.022e23)
	frac_tot=sum(mass_fraction)
	mass_fraction=np.array(mass_fraction)/frac_tot

	ind=np.argmin(abs(kT_ion-kT))
	ind=np.where((kT_ion == kT_ion[ind]))[0]
	ion_frac_here=np.take(ion_frac, ind)

	Z_vec=np.arange(2)
	Z_vec=Z_vec[::-1]
	Az=mass_fraction[0]/pyatomdb.atomic.Z_to_mass(1)
	Azn=mass_fraction[0]/pyatomdb.atomic.Z_to_mass(1)
	for zz in range(nz-1):
		Z=zz+2
		Z_vec=np.arange(Z+1)+1
		Z_vec=Z_vec[::-1]
		Xz=ion_frac_here[Z-1]
		if (math.isnan(Xz[0]) == False):
			sum_vec=sum_vec+sum(Xz*Z_vec)/Z
			atomic_weight=Z*1.007276466879+(math.floor(pyatomdb.atomic.Z_to_mass(Z))-Z)*1.00866491588+sum(Xz*Z_vec)*5.48579909070e-4
			Az=Az+mass_fraction[Z-1]/atomic_weight # Mean atomic weight of ions
			Azn=Azn+sum(Xz*Z_vec)*mass_fraction[Z-1]/atomic_weight #Mean atomic weight of electrons
	mu=1./(Az+Azn)
	mubar=1./Az
	return 1./((mu**3/(mubar**2))*cool_tab_here/kT**2) #kT [K]

# Subroutine: evolution rho
# =========================
def var_rho(kT, kT_tab, cool_tab, direct, kT_ion, ion_frac, sunabund):
	import numpy as np
	import math
	from sun_abund import abund
	import pyatomdb
	from constantes import constante
	import os

	os.environ["ATOMDB"] = 'https://hea-www.cfa.harvard.edu/AtomDB/'
	# write user_data (where filemap and APED)
	fuser = open(direct+"/userdata", 'w')
	fuser.write("filemap="+direct+"/filemap\n")
	fuser.write("atomdbroot=https://hea-www.cfa.harvard.edu/AtomDB/\n")
	fuser.close()
	# read user_data
	setting=pyatomdb.util.load_user_prefs(adbroot=direct)

	ind=np.where(kT_tab <= kT)[0]
	if (kT < min(kT_tab)):
		ind=[0]
	if (kT > max(kT_tab)):
		ind=[int(len(kT_tab))-1]
	cool_tab_here=cool_tab[ind[-1]]/1.e-23

	nz=25
	sum_vec=0.
	number_fraction=abund(ref=sunabund)
	number_fraction=np.array(number_fraction)
	mass_fraction=[]
	for zz in range(nz+1):
		Z=zz+1
		mass_fraction.append(number_fraction[zz]*pyatomdb.atomic.Z_to_mass(Z)/6.022e23)
	frac_tot=sum(mass_fraction)
	mass_fraction=np.array(mass_fraction)/frac_tot

	ind=np.argmin(abs(kT_ion-kT/8.61732e-8))
	ind=np.where((kT_ion == kT_ion[ind]))[0]
	ion_frac_here=np.take(ion_frac, ind)

	Z_vec=np.arange(2)
	Z_vec=Z_vec[::-1]
	Xh=ion_frac_here[0]
	atomic_weight=1.007276466879+sum(Xh*Z_vec)*5.48579909070e-4
	Az=mass_fraction[0]/pyatomdb.atomic.Z_to_mass(1)
	Azn=mass_fraction[0]/pyatomdb.atomic.Z_to_mass(1)
	for zz in range(nz-1):
		Z=zz+2
		Z_vec=np.arange(Z+1)+1
		Z_vec=Z_vec[::-1]
		Xz=ion_frac_here[Z-1]
		if (math.isnan(Xz[0]) == False):
			sum_vec=sum_vec+sum(Xz*Z_vec)/Z
			atomic_weight=Z*1.007276466879+(math.floor(pyatomdb.atomic.Z_to_mass(Z))-Z)*1.00866491588+sum(Xz*Z_vec)*5.48579909070e-4
			Az=Az+mass_fraction[Z-1]/atomic_weight # Mean atomic weight of ions
			Azn=Azn+sum(Xz*Z_vec)*mass_fraction[Z-1]/atomic_weight #Mean atomic weight of electrons
	mu=1./(Az+Azn)
	mubar=1./Az
	return (mu**3/(mubar**2))*cool_tab_here/kT**3

# Subroutine: plot the star radius
# ================================
def plot_star(centerx,centery,radius,x):
	import numpy as np
	return np.sqrt(radius**2-(x-centerx)**2)+centery

def angleyz(y,z):
	import numpy as np
	opp=z
	adj=y
	hypo=np.sqrt(opp**2+adj**2)
	return np.arctan2(opp/hypo,adj/hypo)

def projx(vec,angle):
	import math
	return vec*math.cos(angle)
def projy(vec,angle1,angle2):
	import math
	return vec*math.sin(angle1)*math.cos(angle2)
def projz(vec,angle1,angle2):
	import math
	return vec*math.sin(angle1)*math.sin(angle2)

def apply_mask(triang, x, y, alpha=0.4):
	import numpy as np
	# Mask triangles with sidelength bigger some alpha
	triangles = triang.triangles
	# Mask off unwanted triangles.
	xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
	ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
	maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)
	# apply masking
	triang.set_mask(maxi > alpha)

def find_xy(p1, p2, x):
	import numpy as np
	x1, y1, z1 = p1
	x2, y2, z2 = p2
	resz=[]
	resy=[]
	for i in range(int(len(x))):
		resy.append(np.interp(x[i], (x1, x2), (y1, y2)))
		resz.append(np.interp(x[i], (x1, x2), (z1, z2)))

	return resy, resz

# Main
# ====
def rshock(liste_param):
	from scipy.spatial import KDTree
	from scipy.optimize import fsolve
	from scipy.interpolate import interp1d
	import numpy as np
	import math
	from constantes import constante
	import matplotlib.pyplot as plt
	import matplotlib.tri as tri
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib 
	import os
	import pandas as pd
	import sys
	from astropy.io import ascii
	from raytracing import RT

	Mdot1,Mdot2,vinf1,vinf2, mass_ratio,ex, omega, a, per,R1,R2,beta1,beta2,Teff_1, Teff_2, nbr_points_shock_2D, nbr_points_shock_3D, nstep_width, nbr_bin_profile, direct, cut_lim, wind_prim, wind_sec, phase, incl,ipar, sunabund, mu, M_conj, atom, ion=liste_param

	sys.setrecursionlimit(10000)
        
	pi=math.pi
	fcool = direct+"/cooling_functions/cooling_function_"+str(atom)+str(ion)+".tab"
		                  
	dwind_prim=pd.read_hdf(wind_prim).reset_index(drop=True)
	dwind_sec=pd.read_hdf(wind_sec).reset_index(drop=True)
	tree_prim_wind=KDTree(np.reshape(np.concatenate((dwind_prim.iloc[:]['x'],dwind_prim.iloc[:]['y'])),(-1,int(len(dwind_prim.iloc[:]['x'])))).T,leafsize=1000)
	tree_sec_wind=KDTree(np.reshape(np.concatenate((dwind_sec.iloc[:]['x'],dwind_sec.iloc[:]['y'])),(-1,int(len(dwind_sec.iloc[:]['x'])))).T,leafsize=1000)

	# Model parameters                                             
	# ================                                             
	beta=(Mdot2*vinf2)/(Mdot1*vinf1)
	SumXzZXh=0.99
	mubar=1.3

	Msun=constante('Msun') # solar mass in g
	Rsun=constante('Rsun') # solar radius in cm
	mp=constante('mp') # proton mass in kg
	kb=constante('kb') # Boltzmann constant in erg/K

	R1=R1*Rsun #cm
	R2=R2*Rsun #cm
	a=a*Rsun #cm
	Mdot1=Mdot1*Msun/(365.24*86400.) #g/s
	Mdot2=Mdot2*Msun/(365.24*86400.) #g/s

	vmax=max([vinf1,vinf2]) #km/s
	vec=np.arange(nbr_bin_profile, dtype=np.float)+1.
	vvv=(vec-(math.floor(float(nbr_bin_profile)/2.)+1.))*vmax/math.floor(float(nbr_bin_profile)/2.)

	# Cooling function (lambda(T))
	# ============================
	# computed as cloudy cooling_func.in (lykins13)
	cool_tab=[]
	kT_tab=[]
	with open(fcool) as f:
		for line in f:
			if (line.split()[0] != "#depth"):
				line2=line.split()
				kT_tab.append(float(line2[1]))
				cool_tab.append(float(line2[3]))

	kT_tab=np.array(kT_tab)*8.61732e-8
	plt.plot(kT_tab,np.array(cool_tab)/1.e-23)
	plt.yscale('log')
	plt.xscale('log')
	plt.ylabel("lambda(T) (10^-23 erg cm^3/s)")
	plt.xlabel("Temperature (keV)")
	plt.savefig(direct+"/cooling_func_"+str(atom)+str(ion)+".pdf")
	plt.close()

	# Initialisation                                             
	# ==============

	# Use of the Antokhin's formalism
	# ===============================
	# Equation 28 of Canto96, theta_inf
	thetainf=pi/2.
	if (beta != 1.):
		asymp=pi/(1.-1./beta)
		thetainf=pi/2.
		deltaf=1.-(pi/2.-asymp)/math.tan(pi/2.)
		deltatheta=deltaf/(asymp-thetainf)
		thetainf=thetainf+deltatheta
		while (abs(deltatheta) >= 1.0E-8 and abs(deltaf) >= 1.0E-8):
			deltaf=1.-(thetainf-asymp)/math.tan(thetainf)
			deltatheta=deltaf/(asymp-thetainf)
			thetainf=thetainf+deltatheta
		if (thetainf < pi/2.):
			thetainf=pi/2.+thetainf

	kT_ion=[]
	Z_ion=[]
	ion_frac=[]
	with open(direct+"/ionisation_fraction.tab") as fion:
		for line in fion:
			ligne=line.split()
			if (ligne[0][0] != "#"):
				kT_ion.append(float(ligne[0]))
				Z_ion.append(float(ligne[1]))
				ion_frac_here=[]
				for jion in range(int(ligne[1])+1):
					ion_frac_here.append(float(ligne[jion+2]))
				ion_frac.append(ion_frac_here)
	kT_ion=np.array(kT_ion)

	# Compute the phase we are interested in
	M=phase
	E=M
	dE=(M-E+ex*math.sin(E))/(1.-ex*math.cos(E))
	while (abs(dE) >= 1.0E-8):
		E=E+dE
		dE=(M-E+ex*math.sin(E))/(1.-ex*math.cos(E))
	coi=(math.cos(E)-ex)/(1.-ex*math.cos(E))
	sii=math.sqrt(1.-ex*ex)*math.sin(E)/(1.-ex*math.cos(E))
	# True anomaly 
	phase2=math.atan2(sii,coi)
	if (phase2 == 2.*pi):
		phase2 = 0.
	d=a*(1.-ex**2)/(1.+ex*math.cos(phase2)) #cm
	ouverture_min=math.atan(R1/d)

	x_cdw1, y_cdw1, z_cdw1, x_cdw2, y_cdw2, z_cdw2, x_cdw1g, y_cdw1g, z_cdw1g, x_cdw2g, y_cdw2g, z_cdw2g, kT_all, rho_all, sigma_all=([] for _ in range(15))
	x_cd, y_cd, z_cd, vol_cd, vx_cd, vy_cd, vz_cd, x_skew, y_skew, z_skew, xx_all, yy_all, zz_all, cote_cd=([] for _ in range(14))
	x_cdw1g1_all, x_cdw2g1_all, x_cdw1g1a_all, x_cdw2g1a_all, y_cdw1g1_all, y_cdw2g1_all, y_cdw1g1a_all, y_cdw2g1a_all =([] for _ in range(8))
	vol_cdo, kT_allo, rho_allo, vt1_veco, vt2_veco, sigma_allo =([] for _ in range(6))

	phi_vec=np.arange(nbr_points_shock_3D)*360.*pi/(180.*nbr_points_shock_3D)
	dphi=360.*pi/(180.*nbr_points_shock_3D)
	dzeta_o_vec=np.arange(2.*nstep_width)*0.
	y_o_vec=np.arange(2.*nstep_width)*0.
	crashing='no'
	lmax=0.

	xstar1=np.linspace(-R1/d,R1/d,40)
	ystar1=plot_star(0.0,0.0,R1/d,xstar1)
	xstar2=np.linspace(-R2/d,R2/d,40)
	ystar2=plot_star(0.0,0.0,R2/d,xstar2)
	xstar2=xstar2+1.

	xstag=fsolve(find_x0, (R1+d-R2)/2., args=(d, R1, R2, vinf1, vinf2, Mdot1, Mdot2, tree_prim_wind,tree_sec_wind, dwind_prim, dwind_sec)) #cm 
	xstag=xstag[0]
	fxstag = open(direct+"/histograms/xstag_par"+str(ipar), 'w')
	fxstag.write(str(xstag)+" cm\n")
	fxstag.write(str(xstag/R1)+" R1\n")
	fxstag.write(str(xstag/d)+" d\n")
	fxstag.close()
	xstag=max([xstag,R1])

	if (xstag <= R1 or xstag >= d-R2):
		if (ibal == 0):
			crashing='yes'
			print("")
			print('Crashing wind')
			print("")

	# Equation 24 of Canto96 for a good sampling of the shock
	theta_vec=np.arange(nbr_points_shock_2D)*thetainf/nbr_points_shock_2D+0.0001*pi/2.
	theta1_vec=(np.arange(nbr_points_shock_2D*500.)+1.)*(pi-thetainf)/(nbr_points_shock_2D*500.-1)
	ThetaCotTheta=theta1_vec/np.tan(theta1_vec)
	point=1.+(1./beta)*(theta_vec/np.tan(theta_vec)-1.)
	the1_vec=[]
	for pointhere in point:
		ind=np.where((abs(ThetaCotTheta-pointhere) == np.nanmin(abs(ThetaCotTheta-pointhere))))[0]
		the1_vec.append(theta1_vec[ind[0]])
	the1_vec=np.array(the1_vec)
	r_test=np.sin(the1_vec)/np.sin(the1_vec+theta_vec) #unit of d
	y_vec=r_test*np.sin(theta_vec)*d #cm

	# Equation 6 of Antokhin et al. 2004
	x_old=xstag
	y_old=0.
	x_vec=[]
	eta_vec=[]
	for y_vechere in y_vec:
		r1=math.sqrt(x_old**2+y_old**2)
		r2=math.sqrt((d-x_old)**2+y_old**2)
		eta=comp_eta(x_old, y_old, R1, R2, Mdot1, Mdot2, d, tree_prim_wind, tree_sec_wind, dwind_prim, dwind_sec)
		if (eta != 1):
			pente=(x_old-(d*r1**2*math.sqrt(eta))/(r1**2*math.sqrt(eta)+r2**2))/y_vechere
			x_old=pente*(y_vechere-y_old)+x_old
		y_old=y_vechere
		x_vec.append(x_old/d)	
		eta_vec.append(eta)
	eta_vec=np.array(eta_vec)
	x_vec=np.array(x_vec)

	r_vec=np.sqrt(x_vec**2+(y_vec/d)**2)
	sin_theta=(y_vec/d)/r_vec
	cos_theta=x_vec/r_vec
	theta_vec=np.arctan2(sin_theta,cos_theta) # Equivalent to the delta1 of Antokhin 2004

	adj=x_vec[1:]-x_vec[:-1]
	opp=y_vec[1:]/d-y_vec[:-1]/d
	hypo=np.sqrt(opp**2+adj**2)
	slope_vec=np.arctan2(opp/hypo,adj/hypo) # Equivalent to phi in Antokhin 2004
	first=np.array([slope_vec[0]])
	if (xstag > R1):
		first=np.array([pi/2.])
	slope_vec=np.concatenate((first,slope_vec))
	delta1_vec=slope_vec-theta_vec # Equivalent to theta1 in Antokhin 2004
	delta1_vec[0]=pi/2. # Equivalent to theta1 in Antokhin 2004

	xaxis=[1.,0.,0.]
	x_vec=x_vec*d

	for cut in range(int(len(x_vec))):
		print((str(cut)+" "+str(int(len(x_vec)))))
		if (x_vec[cut]<0. and math.atan(y_vec[cut]/(d-x_vec[cut]))<ouverture_min):
			Rprime=math.tan(ouverture_min)*(d-x_vec[cut])
			y_vec[cut]=abs(Rprime)*float(np.sign(y_vec[cut]))

		thetah=np.arccos(np.dot([x_vec[cut],y_vec[cut],0.],xaxis)/(np.linalg.norm([x_vec[cut],y_vec[cut],0.], axis=0)*np.linalg.norm(xaxis,axis=0)))
		the1_vech=pi-np.arccos(np.dot([x_vec[cut]-d,y_vec[cut],0.],xaxis)/(np.linalg.norm([x_vec[cut]-d,y_vec[cut],0.], axis=0)*np.linalg.norm(xaxis,axis=0)))
		r1_vech=np.sqrt(x_vec[cut]**2+y_vec[cut]**2) #cm
		r2_vech=np.sqrt((d-x_vec[cut])**2+y_vec[cut]**2)
		thetahm1=0.
		the1_vechm1=0.
		slope_vech=pi/2.
		if (cut > 0):
			thetahm1=np.arccos(np.dot([x_vec[cut-1],y_vec[cut-1],0.],xaxis)/(np.linalg.norm([x_vec[cut-1],y_vec[cut-1],0.], axis=0)*np.linalg.norm(xaxis,axis=0)))
			the1_vechm1=pi-np.arccos(np.dot([x_vec[cut-1]-d,y_vec[cut-1],0.],xaxis)/(np.linalg.norm([x_vec[cut-1]-d,y_vec[cut-1],0.], axis=0)*np.linalg.norm(xaxis,axis=0)))
			slope_vech=np.arccos(np.dot([x_vec[cut]-x_vec[cut-1],y_vec[cut]-y_vec[cut-1],0.],xaxis)/(np.linalg.norm([x_vec[cut]-x_vec[cut-1],y_vec[cut]-y_vec[cut-1],0.], axis=0)* np.linalg.norm(xaxis,axis=0)))
		dtheta=abs(thetah-thetahm1)
		dthe1=abs(the1_vech-the1_vechm1)

		# Shock thickness
		l0_1=0. 
		l_1_width=np.arange(nstep_width)*0.
		if (r1_vech>R1):
			l0=np.linspace(0.,(r1_vech/d-R1/d)*0.99,200)
			l0_1=[]
			ind_half=np.where(((dwind_prim['x']**2+dwind_prim['y']**2)**0.5 < r1_vech))[0]
			dwind_prim_choc=dwind_prim.iloc[ind_half]
			ind_half=np.where(((dwind_prim_choc['x']**2+dwind_prim_choc['y']**2)**0.5 > R1/d))[0]
			dwind_prim_choc=dwind_prim_choc.iloc[ind_half]
			ind_half=np.where((dwind_prim_choc['y'] < y_vec[cut]))[0]
			dwind_prim_choc=dwind_prim_choc.iloc[ind_half]
			for il0 in range (200):
				# closest point
				xh=projx(r1_vech-l0[il0]*d,thetah)
				yh=projy(r1_vech-l0[il0]*d,thetah,0.)
				rien, ind_close = tree_prim_wind.query([[xh/d,yh/d]])
				ind_close=ind_close[0]
				v=max(vinf1,np.sqrt(dwind_prim.iloc[ind_close]['ux']**2+dwind_prim.iloc[ind_close]['uy']**2))
				if (dwind_prim.iloc[ind_close]['ux'] == 0 and dwind_prim.iloc[ind_close]['uy']==0):
					ind_closep=ind_close+1
					ind_closem=ind_close-1
					stade="up"
					dw1=dwind_prim.iloc[ind_closep]['ux']
					dw2=dwind_prim.iloc[ind_closep]['uy']
					ind_closep=ind_close+1
					while (dw1 == 0 and dw2==0):
						if (stade=="up"):
							dw1=dwind_prim.iloc[ind_closem]['ux']
							dw2=dwind_prim.iloc[ind_closem]['uy']
							ind_closem=ind_closem-1
							ind_close=ind_closem
							stade="down"
						else:
							stade="up"
							dw1=dwind_prim.iloc[ind_closep]['ux']
							dw2=dwind_prim.iloc[ind_closep]['uy']
							ind_closep=ind_closep+1
							ind_close=ind_closep
					v=max(vinf1,np.sqrt(dwind_prim.iloc[ind_close]['ux']**2+dwind_prim.iloc[ind_close]['uy']**2))
				# Temperature of the postschock gas (eq 18 of Antokhin 2004)
				kT0=1.21*(mu/0.62)*(v*math.sin(slope_vech-thetah)/1.E8)**2 #keV
				# Cooling function
				ind=np.where(kT_tab <= kT0)[0]
				if (kT0 < min(kT_tab)):
					ind=[0]
				if (kT0 > max(kT_tab)):
					ind=[int(len(kT_tab))-1]
				cool_tab_here=cool_tab[ind[-1]]
				rho=4.*Mdot1/(v*4.*pi*(r1_vech-l0[il0]*d)**2) #g/cm^3
				# Width of the shock
				l0_1.append(15.*mubar**2*(mp*1000.)**2*abs(v*math.sin(slope_vech-thetah))**3/(d*512.*SumXzZXh*rho*cool_tab_here)) #units of d
			ind=np.where(abs(l0_1-l0) == min(abs(l0_1-l0)))[0]
			l0_1=l0_1[ind[0]]*d
			l_1_width=(np.arange(nstep_width))*((2.*nstep_width-1)*l0_1/(2.*nstep_width))/(nstep_width-1)+l0_1/(2.*nstep_width)

		l0_2=0.
		l_2_width=np.arange(nstep_width)*0.
		if (abs(math.sin(pi-slope_vech-the1_vech)) > 1.0e-8):
			l0=np.linspace(0.,(r2_vech/d-R2/d)*0.99,200)
			l0_2=[]
			ind_half=np.where(((dwind_sec['x']**2+dwind_sec['y']**2)**0.5 < r2_vech))[0]
			dwind_sec_choc=dwind_sec.iloc[ind_half]
			ind_half=np.where(((dwind_sec_choc['x']**2+dwind_sec_choc['y']**2)**0.5 > R2/d))[0]
			dwind_sec_choc=dwind_sec_choc.iloc[ind_half]
			ind_half=np.where((dwind_sec_choc['y'] < y_vec[cut]))[0]
			dwind_sec_choc=dwind_sec_choc.iloc[ind_half]
			for il0 in range (200):
				# closest point
				xh=projx(r2_vech-l0[il0]*d,pi-the1_vech)+d
				yh=projy(r2_vech-l0[il0]*d,pi-the1_vech,0.)
				rien, ind_close = tree_sec_wind.query([[(d-xh)/d,yh/d]])
				ind_close=ind_close[0]
				v=np.sqrt(dwind_sec_choc.iloc[ind_close]['ux']**2+dwind_sec_choc.iloc[ind_close]['uy']**2)
				if (dwind_sec.iloc[ind_close]['ux'] == 0 and dwind_sec.iloc[ind_close]['uy']==0):
					ind_closep=ind_close+1
					ind_closem=ind_close-1
					stade="up"
					dw1=dwind_sec.iloc[ind_closep]['ux']
					dw2=dwind_sec.iloc[ind_closep]['uy']
					ind_closep=ind_close+1
					while (dw1 == 0 and dw2==0):
						if (stade=="up"):
							dw1=dwind_sec.iloc[ind_closem]['ux']
							dw2=dwind_sec.iloc[ind_closem]['uy']
							ind_closem=ind_closem-1
							ind_close=ind_closem
							stade="down"
						else:
							stade="up"
							dw1=dwind_sec.iloc[ind_closep]['ux']
							dw2=dwind_sec.iloc[ind_closep]['uy']
							ind_closep=ind_closep+1
							ind_close=ind_closep
					v=np.sqrt(dwind_sec.iloc[ind_close]['ux']**2+dwind_sec.iloc[ind_close]['uy']**2)

				if (v> vinf2):
					v=vinf2
				# Temperature of the postschock gas (eq 18 of Antokhin 2004)
				kT0=1.21*(mu/0.62)*(v*math.sin(pi-slope_vech-the1_vech)/1.E8)**2 #keV
				# Cooling function
				ind=np.where(kT_tab <= kT0)[0]
				if (kT0 < min(kT_tab)):
					ind=[0]
				if (kT0 > max(kT_tab)):
					ind=[int(len(kT_tab))-1]
				cool_tab_here=cool_tab[ind[-1]]
				rho=4.*Mdot2/(v*4.*pi*(r2_vech-l0[il0]*d)**2) #g/cm^3
				# Width of the shock
				l0_2.append(15.*mubar**2*(mp*1000.)**2*abs(v*math.sin(pi-slope_vech-the1_vech))**3/(d*512.*SumXzZXh*rho*cool_tab_here)) #units of d
			ind=np.where(abs(l0_2-l0) == min(abs(l0_2-l0)))[0]
			l0_2=l0_2[ind[0]]*d
			l_2_width=(np.arange(nstep_width))*((2.*nstep_width-1)*l0_2/(2.*nstep_width))/(nstep_width-1)+l0_2/(2.*nstep_width)
		
		if (cut >= 1):
			if (r1_vech-l0_1<rl1):
				l0_1=l011
			if (r2_vech-l0_2<rl2):
				l0_2=l022
		lmax=max(lmax,l0_1,l0_2)

		xh=projx(r1_vech-l0[il0],thetah)
		yh=projy(r1_vech-l0[il0],thetah,0.)
		rien, ind_close = tree_prim_wind.query([[xh/d,yh/d]])
		ind_close=ind_close[0]
		v1=min(vinf1,np.sqrt(dwind_prim.iloc[ind_close]['ux']**2+dwind_prim.iloc[ind_close]['uy']**2))
		if (dwind_prim.iloc[ind_close]['ux'] == 0 and dwind_prim.iloc[ind_close]['uy']==0):
			ind_closep=ind_close+1
			ind_closem=ind_close-1
			stade="up"
			dw1=dwind_prim.iloc[ind_closep]['ux']
			dw2=dwind_prim.iloc[ind_closep]['uy']
			ind_closep=ind_close+1
			while (dw1 == 0 and dw2==0):
				if (stade=="up"):
					dw1=dwind_prim.iloc[ind_closem]['ux']
					dw2=dwind_prim.iloc[ind_closem]['uy']
					ind_closem=ind_closem-1
					ind_close=ind_closem
					stade="down"
				else:
					stade="up"
					dw1=dwind_prim.iloc[ind_closep]['ux']
					dw2=dwind_prim.iloc[ind_closep]['uy']
					ind_closep=ind_closep+1
					ind_close=ind_closep
			v1=min(vinf1,np.sqrt(dwind_prim.iloc[ind_close]['ux']**2+dwind_prim.iloc[ind_close]['uy']**2))
		v1=max(100.e5,v1)

		xh=projx(r2_vech-l0[il0],pi-the1_vech)+d
		yh=projy(r2_vech-l0[il0],pi-the1_vech,0.)
		rien, ind_close = tree_sec_wind.query([[(d-xh)/d,yh/d]])
		ind_close=ind_close[0]
		v2=min(vinf2,np.sqrt(dwind_sec.iloc[ind_close]['ux']**2+dwind_sec.iloc[ind_close]['uy']**2))
		if (dwind_sec.iloc[ind_close]['ux'] == 0 and dwind_sec.iloc[ind_close]['uy']==0):
			ind_closep=ind_close+1
			ind_closem=ind_close-1
			stade="up"
			dw1=dwind_sec.iloc[ind_closep]['ux']
			dw2=dwind_sec.iloc[ind_closep]['uy']
			ind_closep=ind_close+1
			while (dw1 == 0 and dw2==0):
				if (stade=="up"):
					dw1=dwind_sec.iloc[ind_closem]['ux']
					dw2=dwind_sec.iloc[ind_closem]['uy']
					ind_closem=ind_closem-1
					ind_close=ind_closem
					stade="down"
				else:
					stade="up"
					dw1=dwind_sec.iloc[ind_closep]['ux']
					dw2=dwind_sec.iloc[ind_closep]['uy']
					ind_closep=ind_closep+1
					ind_close=ind_closep
			v2=min(vinf2,np.sqrt(dwind_sec.iloc[ind_close]['ux']**2+dwind_sec.iloc[ind_close]['uy']**2))
		v2=max(100.e5,v2)
				
		increasevol=np.arange(nstep_width)
		vol_cdo=np.concatenate((vol_cdo,l0_1*(r1_vech-l0_1*increasevol/nstep_width)**2*math.sin(pi/2.-thetah)*dtheta*dphi/(abs(math.sin(slope_vech-thetah))*nstep_width)))
		vol_cdo=np.concatenate((vol_cdo,l0_2*(r2_vech-l0_2*increasevol/nstep_width)**2*math.sin(pi/2.-the1_vec)*dthe1*dphi/(abs(math.sin(slope_vech-the1_vech))*nstep_width)))

		x_cdw1g1=projx(r1_vech-l0_1,thetah)
		y_cdw1g1=projy(r1_vech-l0_1,thetah,0.)
		x_cdw2g1=projx(r2_vech-l0_2,pi-the1_vech)+d
		y_cdw2g1=projy(r2_vech-l0_2,pi-the1_vech,0.)
		x_cdw1g1a=projx(r1_vech-l_1_width,thetah)
		y_cdw1g1a=projy(r1_vech-l_1_width,thetah,0.)
		x_cdw2g1a=projx(r2_vech-l_2_width,pi-the1_vech)+d
		y_cdw2g1a=projy(r2_vech-l_2_width,pi-the1_vech,0.)

		slope_vech=pi/2.
		vt1=0.
		vt2=0.
		if (cut >= 1):
			slope_vech=np.arccos(np.dot([x_cdw1g1-x_cdw1g11,y_cdw1g1-y_cdw1g11,0.],xaxis)/(np.linalg.norm([x_cdw1g1-x_cdw1g11,y_cdw1g1-y_cdw1g11,0.], axis=0)* np.linalg.norm(xaxis,axis=0)))
			vt1=v1*math.cos(slope_vech-thetah)
			slope_vech=np.arccos(np.dot([x_cdw2g1-x_cdw2g11,y_cdw2g1-y_cdw2g11,0.],xaxis)/(np.linalg.norm([x_cdw2g1-x_cdw2g11,y_cdw2g1-y_cdw2g11,0.], axis=0)* np.linalg.norm(xaxis,axis=0)))
			vt2=v2*math.cos(pi-slope_vech-the1_vech)
		x_cdw1g11=x_cdw1g1
		y_cdw1g11=y_cdw1g1
		x_cdw2g11=x_cdw2g1
		y_cdw2g11=y_cdw2g1
		rl1=r1_vech-l0_1
		rl2=r2_vech-l0_2
		l011=l0_1
		l022=l0_2 

		vp1=math.sqrt(v1**2-vt1**2)
		vp2=math.sqrt(v2**2-vt2**2)

		vt1_veco=np.concatenate((vt1_veco,[vt1]))
		vt2_veco=np.concatenate((vt2_veco,[vt2]))

		#Evolution of the temperature and density inside the shock
		kT0_1=0.
		l_1_width=np.arange(nstep_width)*0.
		kT1=np.arange(nstep_width)*0. #keV
		rho1_vec=4.*Mdot1/(v1*4.*pi*r1_vech**2)+np.arange(nstep_width)*0. #g/cm^3
		if (l0_1>0.):
			kT0_1=1.21*(mu/0.62)*(vp1/1.E8)**2 #keV
			l_1_width=np.linspace(l0_1/(2.*nstep_width),l0_1,nstep_width) #cm
			rho0=4.*Mdot1/(v1*4.*pi*(r1_vech-l0_1)**2) #g/cm^3
			kT_width=np.linspace(1.,kT0_1/8.61732e-8,nstep_width)
			const=1./(9.*0.99*mp*1000.*v1**3*rho0/(40.*kb**3))
			test=[]
			for i in kT_width:
				test.append(1.e23*integrand(i,kT_tab, cool_tab, direct, kT_ion/8.61732e-8, ion_frac, sunabund))
			test=(const*(np.cumsum(test)*np.mean(kT_width[1:]-kT_width[:-1])))/d
			test=test*l0_1/max(test)
			f = interp1d(test,kT_width*8.61732e-8, kind='cubic',assume_sorted=True)
			whereto=l_1_width*0.+l_1_width
			ind=np.where((whereto<min(test)))[0]
			if (len(ind)>0):
				whereto[ind]=[min(test)]*len(ind)
			ind=np.where((whereto>max(test)))[0]
			if (len(ind)>0):
				whereto[ind]=[max(test)]*len(ind)
			kT1=f(whereto)

			const3=[]
			rho1_vec=[rho0]
			for i in range(int(len(kT1))):
				if (i > 0):
					const2=1.e23*var_rho(kT1[i]/8.61732e-8, kT_tab, cool_tab, direct, kT_ion/8.61732e-8, ion_frac, sunabund)
					const3.append(const2*abs(l_1_width[i-1]-l_1_width[i])/const)
			const3=np.array(const3)/(max(const3)*1.1)
			for i in range(int(len(kT1))):
				if (i > 0):
					drho=const3[i-1]*rho1_vec[-1]/(1.-const3[i-1])
					rho1_vec.append(rho1_vec[-1]+drho)
			#rho_width=np.linspace(rho0,rho0*1e10,nstep_width)
			#const=1.2**cool_tab_here/(rho0**3*v1**3*1.3*1.6734E-24)
			#rho1_vec=1/(4.*const*l_1_width)**0.25-1/(4.*const*l_1_width[-1])**0.25+4.*rho0

		kT0_2=0.
		l_2_width=np.arange(nstep_width)*0.
		kT2=np.arange(nstep_width)*0. #keV
		rho2_vec=4.*Mdot2/(v2*4.*pi*r2_vech**2)+np.arange(nstep_width)*0. #g/cm^3
		if (l0_2>0.):
			kT0_2=1.21*(mu/0.62)*(vp2/1.E8)**2 #keV
			l_2_width=np.linspace(l0_2/(2.*nstep_width),l0_2,nstep_width) #cm
			rho0=4.*Mdot2/(v2*4.*pi*(r2_vech-l0_2)**2) #g/cm^3
			kT_width=np.linspace(1.,kT0_2/8.61732e-8,nstep_width)
			const=1./(9.*0.99*mp*1000.*v2**3*rho0/(40.*kb**3))
			test=[]
			for i in kT_width:
				test.append(1.e23*integrand(i,kT_tab, cool_tab, direct, kT_ion/8.61732e-8, ion_frac, sunabund))
			test=(const*(np.cumsum(test)*np.mean(kT_width[1:]-kT_width[:-1])))/d
			test=test*l0_2/max(test)
			f = interp1d(test,kT_width*8.61732e-8, kind='cubic',assume_sorted=True)
			whereto=l_2_width*0.+l_2_width
			ind=np.where((whereto<min(test)))[0]
			if (len(ind)>0):
				whereto[ind]=[min(test)]*len(ind)
			ind=np.where((whereto>max(test)))[0]
			if (len(ind)>0):
				whereto[ind]=[max(test)]*len(ind)
			kT2=f(whereto)

			const3=[]
			rho2_vec=[rho0]
			for i in range(int(len(kT2))):
				if (i > 0):
					const2=1.e23*var_rho(kT2[i]/8.61732e-8, kT_tab, cool_tab, direct, kT_ion/8.61732e-8, ion_frac, sunabund)
					const3.append(const2*abs(l_2_width[i-1]-l_2_width[i])/const)
			const3=np.array(const3)/(max(const3)*1.1)
			for i in range(int(len(kT2))):
				if (i > 0):
					drho=const3[i-1]*rho2_vec[-1]/(1.-const3[i-1])
					rho2_vec.append(rho2_vec[-1]+drho)

		eta_vec[cut]=max([eta_vec[cut], 1./eta_vec[cut]])
		#Evolution of the surface density
		if (kT0_1 > 0):
			for isigma in range(int(len(l_1_width))):
				sigma=sigmas(projy((r1_vech-l_1_width[isigma]),thetah,0.), 0., xstag, y_vec, thetah, slope_vech-thetah, slope_vech, Mdot1, v1, eta_vec[cut], rho1_vec[isigma], beta1, 1, dzeta_o_vec[isigma], y_o_vec[isigma],R1,R2,d,R1,crashing)
				if (sigma[1] != -1.0):
					if (len(sigma) == 3):
						dzeta_o_vec[isigma]=sigma[1]
						y_o_vec[isigma]=sigma[2]
					else:
						dzeta_o_vec[isigma]=sigma[1][0]
						y_o_vec[isigma]=sigma[2]
				while (not isinstance(sigma,float)):
					sigma=sigma[0]
				sigma_allo=np.concatenate((sigma_allo,[sigma]))
		else:
			sigma=4.*Mdot1/(v1*4.*pi*(r1_vech-l0_1)**2) #g/cm^3
			for isigma in range(int(len(l_1_width))):
				sigma_allo=np.concatenate((sigma_allo,[sigma]))
				dzeta_o_vec[isigma]=sigma*projy(R1,thetah,0.)*v1*np.cos(slope_vech-thetah)
				y_o_vec[isigma]=projy(R1,thetah,0.)


		if (kT0_2 > 0):
			for isigma in range(int(len(l_2_width))):
				sigma=sigmas(projy((r2_vech-l_2_width[isigma]),the1_vech,0.), 0., xstag, y_vec, the1_vech, pi-slope_vech-the1_vech, slope_vech, Mdot2, v2, eta_vec[cut], rho2_vec[isigma], beta2, 2, dzeta_o_vec[isigma+int(len(l_1_width))], y_o_vec[isigma+int(len(l_1_width))],R1,R2,d,R2,crashing)
				if (sigma[1] != -1.0):
					if (len(sigma) == 3):
						dzeta_o_vec[isigma+nstep_width]=sigma[1]
						y_o_vec[isigma+nstep_width]=sigma[2]
					else:
						dzeta_o_vec[isigma+nstep_width]=sigma[1][0]
						y_o_vec[isigma+nstep_width]=sigma[2]
				while (not isinstance(sigma,float)):
					sigma=sigma[0]
				sigma_allo=np.concatenate((sigma_allo,[sigma]))
		else:
			sigma=4.*Mdot2/(v2*4.*pi*(r2_vech-l0_2)**2) #g/cm^3
			for isigma in range(int(len(l_2_width))):
				sigma_allo=np.concatenate((sigma_allo,[sigma]))
				dzeta_o_vec[isigma+nstep_width]=sigma*projy(R2,the1_vech,0.)*v2*np.cos(pi-slope_vech-the1_vech)
				y_o_vec[isigma+nstep_width]=projy(R2,the1_vech,0.)

		x_cdw1g1_all=np.concatenate((x_cdw1g1_all,[x_cdw1g1]))
		y_cdw1g1_all=np.concatenate((y_cdw1g1_all,[y_cdw1g1]))
		x_cdw2g1_all=np.concatenate((x_cdw2g1_all,[x_cdw2g1]))
		y_cdw2g1_all=np.concatenate((y_cdw2g1_all,[y_cdw2g1]))
		x_cdw1g1a_all=np.concatenate((x_cdw1g1a_all,x_cdw1g1a))
		y_cdw1g1a_all=np.concatenate((y_cdw1g1a_all,y_cdw1g1a))
		x_cdw2g1a_all=np.concatenate((x_cdw2g1a_all,x_cdw2g1a))
		y_cdw2g1a_all=np.concatenate((y_cdw2g1a_all,y_cdw2g1a))
		kT_allo=np.concatenate((kT_allo,kT1))
		rho_allo=np.concatenate((rho_allo,rho1_vec))
		kT_allo=np.concatenate((kT_allo,kT2))
		rho_allo=np.concatenate((rho_allo,rho2_vec))


	for icut in range(nbr_points_shock_3D):
		if (float(icut+1.)%100. == 0.):
			sentence="* "+str(icut+1)+"/"+str(nbr_points_shock_3D)+" *"
			print(sentence)
		xcut=x_vec
		ycut=y_vec*np.cos(phi_vec[icut])
		zcut=y_vec*np.sin(phi_vec[icut])
		y_cdw1g1=y_cdw1g1_all*np.cos(phi_vec[icut])
		z_cdw1g1=y_cdw1g1_all*np.sin(phi_vec[icut])
		y_cdw2g1=y_cdw2g1_all*np.cos(phi_vec[icut])
		z_cdw2g1=y_cdw2g1_all*np.sin(phi_vec[icut])
		y_cdw1g1a=y_cdw1g1a_all*np.cos(phi_vec[icut])
		z_cdw1g1a=y_cdw1g1a_all*np.sin(phi_vec[icut])
		y_cdw2g1a=y_cdw2g1a_all*np.cos(phi_vec[icut])
		z_cdw2g1a=y_cdw2g1a_all*np.sin(phi_vec[icut])

		# Discontinuity surface
		# named _skew by historical reasons by not skewing applied
		x_skew=np.concatenate((x_skew,xcut)) #cm
		y_skew=np.concatenate((y_skew,ycut))
		z_skew=np.concatenate((z_skew,zcut))

		# shock surface
		x_cdw1g=np.concatenate((x_cdw1g,x_cdw1g1_all))
		y_cdw1g=np.concatenate((y_cdw1g,y_cdw1g1))
		z_cdw1g=np.concatenate((z_cdw1g,z_cdw1g1))
		x_cdw2g=np.concatenate((x_cdw2g,x_cdw2g1_all))
		y_cdw2g=np.concatenate((y_cdw2g,y_cdw2g1))
		z_cdw2g=np.concatenate((z_cdw2g,z_cdw2g1))

		# Overall points inside the shock
		xx_all=np.concatenate((xx_all,x_cdw1g1a_all)) #cm
		yy_all=np.concatenate((yy_all,y_cdw1g1a))
		zz_all=np.concatenate((zz_all,z_cdw1g1a))
		xx_all=np.concatenate((xx_all,x_cdw2g1a_all)) #cm
		yy_all=np.concatenate((yy_all,y_cdw2g1a))
		zz_all=np.concatenate((zz_all,z_cdw2g1a))


		for cut in range(int(len(x_vec))):
			slope_vech=pi/2.
			if (cut >= 1):
				slope_vech=np.arccos(np.dot([x_cdw1g1_all[cut]-x_cdw1g1_all[cut-1],y_cdw1g1[cut]-y_cdw1g1[cut-1],z_cdw1g1[cut]-z_cdw1g1[cut-1]],xaxis)/(np.linalg.norm([x_cdw1g1_all[cut]-x_cdw1g1_all[cut-1],y_cdw1g1[cut]-y_cdw1g1[cut-1],z_cdw1g1[cut]-z_cdw1g1[cut-1]], axis=0)* np.linalg.norm(xaxis,axis=0)))
			lat=angleyz(y_cdw1g1[cut],z_cdw1g1[cut])

			vx_cd=np.concatenate((vx_cd,projx(np.arange(nstep_width)*0.+vt1_veco[cut],slope_vech))) #cm
			vy_cd=np.concatenate((vx_cd,projy(np.arange(nstep_width)*0.+vt1_veco[cut],slope_vech,lat)))
			vz_cd=np.concatenate((vx_cd,projz(np.arange(nstep_width)*0.+vt1_veco[cut],slope_vech,lat)))

			slope_vech=pi/2.
			if (cut >= 1):
				slope_vech=np.arccos(np.dot([x_cdw2g1_all[cut]-x_cdw2g1_all[cut-1],y_cdw2g1[cut]-y_cdw2g1[cut-1],z_cdw2g1[cut]-z_cdw2g1[cut-1]],xaxis)/(np.linalg.norm([x_cdw2g1_all[cut]-x_cdw2g1_all[cut-1],y_cdw2g1[cut]-y_cdw2g1[cut-1],z_cdw2g1[cut]-z_cdw2g1[cut-1]], axis=0)* np.linalg.norm(xaxis,axis=0)))
			lat=angleyz(y_cdw2g1[cut],z_cdw2g1[cut])
			vx_cd=np.concatenate((vx_cd,projx(np.arange(nstep_width)*0.+vt2_veco[cut],slope_vech))) #cm
			vy_cd=np.concatenate((vy_cd,projy(np.arange(nstep_width)*0.+vt2_veco[cut],slope_vech,lat)))
			vz_cd=np.concatenate((vz_cd,projz(np.arange(nstep_width)*0.+vt2_veco[cut],slope_vech,lat)))

			cote_cd=np.concatenate((cote_cd,1+0.*np.arange(nstep_width)))
			cote_cd=np.concatenate((cote_cd,2+0.*np.arange(nstep_width)))
			kT_all=np.concatenate((kT_all,kT_allo[int(cut*2.*nstep_width):int((cut+1)*2.*nstep_width)]))
			rho_all=np.concatenate((rho_all,rho_allo[int(cut*2.*nstep_width):int((cut+1)*2.*nstep_width)]))
			vol_cd=np.concatenate((vol_cd,vol_cdo[int(cut*2.*nstep_width):int((cut+1)*2.*nstep_width)]))
			sigma_all=np.concatenate((sigma_all,sigma_allo[int(cut*2.*nstep_width): int((cut+1)*2.*nstep_width)]))


	x_vec=np.concatenate((x_skew,x_cd)) #cm
	y_vec=np.concatenate((y_skew,y_cd))
	z_vec=np.concatenate((z_skew,z_cd))
	x_cdw1=np.concatenate((x_cdw1g,x_cdw1))
	y_cdw1=np.concatenate((y_cdw1g,y_cdw1))
	z_cdw1=np.concatenate((z_cdw1g,z_cdw1))
	x_cdw2=np.concatenate((x_cdw2g,x_cdw2))
	y_cdw2=np.concatenate((y_cdw2g,y_cdw2))
	z_cdw2=np.concatenate((z_cdw2g,z_cdw2))

	### Plots ###

	x_cap=np.array(xx_all[:])
	y_cap=np.array(yy_all[:])
	kT_cap=np.array(kT_all[:])
	rho_cap=np.array(rho_all[:])

	nbr_layer=int(len(x_cap)/nbr_points_shock_3D)
	x_layer=np.concatenate((x_cap[0:nbr_layer], x_cap[int(round(nbr_points_shock_3D*nbr_layer*0.5)):int(round(nbr_points_shock_3D*nbr_layer*0.5)+nbr_layer)]))
	y_layer=np.concatenate((y_cap[0:nbr_layer], y_cap[int(round(nbr_points_shock_3D*nbr_layer*0.5)):int(round(nbr_points_shock_3D*nbr_layer*0.5)+nbr_layer)]))
	kT_layer=np.concatenate((kT_cap[0:nbr_layer], kT_cap[int(round(nbr_points_shock_3D*nbr_layer*0.5)):int(round(nbr_points_shock_3D*nbr_layer*0.5)+nbr_layer)]))
	rho_layer=np.concatenate((rho_cap[0:nbr_layer], rho_cap[int(round(nbr_points_shock_3D*nbr_layer*0.5)):int(round(nbr_points_shock_3D*nbr_layer*0.5)+nbr_layer)]))

	circle1 = plt.Circle((0, 0), R1/d, color='white')
	circle2 = plt.Circle((1, 0), R2/d, color='white')
	
	levels = np.arange(10.)*math.ceil(max(kT_layer))/9.
	triang2 = tri.Triangulation(x_layer/d,y_layer/d)
	apply_mask(triang2, x_layer/d,y_layer/d, alpha=0.3)
	plt.figure(1)
	fig = plt.gcf()
	ax = fig.gca()
	ax.add_artist(circle2)
	plt.tricontourf(triang2, kT_layer,levels=levels)
	ax.add_artist(circle1)
	plt.plot(xstar1,ystar1,c='black')
	plt.plot(xstar2,ystar2,c='black')
	plt.plot(xstar1,-ystar1,c='black')
	plt.plot(xstar2,-ystar2,c='black')
	cbar=plt.colorbar(format='%.2f')
	cbar.set_label(r'Temperature [keV]')
	cbar.set_clim(0,max(kT_all))
	plt.ylabel("y/D")
	plt.xlabel("x/D")
	plt.xlim(-1,2)
	plt.ylim(-2,1)
	plt.savefig(direct+"/plots/temperature_xy_par"+str(ipar)+".pdf")
	plt.close()
	ax = 0.

	factor=1.e-13
	while (max(rho_layer)/factor < 1):
		factor = factor/10.
	circle1 = plt.Circle((0, 0), R1/d, color='white')
	circle2 = plt.Circle((1, 0), R2/d, color='white')
	levels = np.sort(np.arange(20.)*(math.ceil(max(rho_layer)/factor)-math.floor(min(rho_layer)/factor))/19.+math.floor(min(rho_layer)/factor))
	triang2 = tri.Triangulation(x_layer/d,y_layer/d)
	apply_mask(triang2, x_layer/d,y_layer/d, alpha=0.3)
	plt.figure(1)
	fig = plt.gcf()
	ax = fig.gca()
	ax.add_artist(circle2)
	plt.tricontourf(triang2, rho_layer/factor,levels=levels)
	ax.add_artist(circle1)
	plt.plot(xstar1,ystar1,c='black')
	plt.plot(xstar2,ystar2,c='black')
	plt.plot(xstar1,-ystar1,c='black')
	plt.plot(xstar2,-ystar2,c='black')
	cbar=plt.colorbar(format='%.2f')
	cbar.set_label(r'Density ['+str(factor)+' g/cm^3]')
	cbar.set_clim(min(rho_all)/factor,max(rho_all)/factor)
	plt.ylabel("y/D")
	plt.xlabel("x/D")
	plt.xlim(-1,2)
	plt.ylim(-2,1)
	plt.savefig(direct+"/plots/density_par"+str(ipar)+".pdf")
	plt.close()

	table = {'x': xx_all, 'y': yy_all, 'z': zz_all, 'kT': kT_all, 'rho': rho_all}
	ascii.write(table, direct+"/plots/char_shock_par"+str(ipar)+".dat")

	nbr_pt_tot=int(len(xx_all))

	deltaZZ=0.05*lmax #cm
	ndeltaZZ=min([100000.,math.ceil(math.sqrt((max(x_vec)-min(x_vec))**2+(max(y_vec)-min(y_vec))**2+(max(z_vec)-min(z_vec))**2)/deltaZZ)]) # max 100 000 steps
	xmin1=min(x_cdw1)
	xmax2=1.1*max(x_cdw2)
	ymin1=min(y_cdw2)
	ymax2=1.1*max(y_cdw2)

	ensemble_points=np.reshape(np.concatenate((xx_all/d,yy_all/d,zz_all/d)),(-1,int(len(xx_all)))).T
	del xx_all
	del yy_all
	del zz_all
	tree = KDTree(ensemble_points,leafsize=1000)
	treecdw1=KDTree(np.reshape(np.concatenate((x_cdw1/d,y_cdw1/d,z_cdw1/d)),(-1,int(len(x_cdw1)))).T,leafsize=1000)
	treecdw2=KDTree(np.reshape(np.concatenate((x_cdw2/d,y_cdw2/d,z_cdw2/d)),(-1,int(len(x_cdw2)))).T,leafsize=1000)
	treevec=KDTree(np.reshape(np.concatenate((x_vec/d,y_vec/d,z_vec/d)),(-1,int(len(x_vec)))).T,leafsize=1000)
	del x_cdw1
	del y_cdw1
	del z_cdw1
	del x_cdw2
	del y_cdw2
	del z_cdw2

	phase=phase-M_conj
	if (phase<0):
		phase=2.*pi+phase

	deltaX=deltaZZ*math.cos(phase)*math.sin(incl) #cm
	deltaY=deltaZZ*math.cos(incl)
	deltaZ=-deltaZZ*math.sin(phase)*math.sin(incl)
	dir_obs=np.array([deltaX,deltaY,deltaZ])/np.linalg.norm(np.array([deltaX,deltaY,deltaZ]))

	temp_emiss=[]
	emiss_part1=[]
	for i in range(nbr_bin_profile):
		temp_emiss.append([])
		emiss_part1.append([])
		fRT = open(direct+"/histograms/ray_tracing_par"+str(ipar)+"_bin"+str(i)+".data", 'w')
		fRT.write("# Temperature / Factor\n")
		fRT.close()

	### Line profile ###
	print("")
	print("***********")
	print("Ray tracing")
	print("***********")
	for i in range(nbr_pt_tot):
		if (float(i+1.)%100. == 0.):
			sentence="* "+str(i+1)+"/"+str(nbr_pt_tot)+" *"
			print(sentence)
		xx=ensemble_points[i,0]*d
		yy=ensemble_points[i,1]*d
		zz=ensemble_points[i,2]*d
		kT=kT_all[i]
		rho=rho_all[i]
		vol=abs(vol_cd[i])
		vrad=vx_cd[i]*math.cos(phase)*math.sin(incl)+vy_cd[i]*math.cos(incl)-vz_cd[i]*math.sin(incl)*math.sin(phase)

		yt1=math.sin(phase)*xx+math.cos(phase)*zz #cm
		zt1=math.cos(phase)*math.sin(incl)*xx+math.cos(incl)*yy-math.sin(phase)*math.sin(incl)*zz
		pt1=-math.cos(incl)*math.cos(phase)*xx-math.sin(incl)*yy+math.cos(incl)*math.sin(phase)*zz
		yt2=math.sin(phase)*(xx-1.)+math.cos(phase)*zz #cm
		zt2=math.cos(phase)*math.sin(incl)*(xx-1.)+math.cos(incl)*yy-math.sin(phase)*math.sin(incl)*zz
		pt2=-math.cos(incl)*math.cos(phase)*(xx-1.)-math.sin(incl)*yy+math.cos(incl)*math.sin(phase)*zz

		vis1=1.
		if (zt1 < 0. and np.sqrt(pt1**2+yt1**2) < R1):
			vis1=0.
		vis2=1.
		if (zt2 < 0. and np.sqrt(pt2**2+yt2**2) < R2):
			vis2=0.

		if (vis1*vis2 > 0.):
			ind_vrad=np.where(vvv <= vrad)[0]
			
			if (len(ind_vrad) > 0 and vol != 0.):
				temp_emiss[ind_vrad[-1]].append(kT)
				emiss_part1[ind_vrad[-1]].append(1.0e-23*(1.17*rho**2/(1.34*(mp*1000.))**2)*vol/1.0e27)

				ray_tracing=RT(ndeltaZZ, xx, yy, zz, deltaZZ, deltaX, deltaY, deltaZ, x_vec, y_vec, z_vec, d, cote_cd, kT_all, sigma_all, Mdot1, Mdot2, R1, R2, Teff_1, Teff_2, tree_prim_wind, tree_sec_wind, dwind_prim, dwind_sec, dir_obs, mp, zt1, pt1, yt1, zt2, pt2, yt2, treecdw1, treecdw2, tree, treevec, xmin1, xmax2, ymin1, ymax2)

				fRT = open(direct+"/histograms/ray_tracing_par"+str(ipar)+"_bin"+str(ind_vrad[-1])+".data", 'a')
				fRT.write(' '.join(map(str, ray_tracing))+"\n")
				fRT.close()

	if (crashing == 'yes'):
		fbin = open(direct+"/histograms/bin_par"+str(ipar)+"_crashing.data", 'w')
	else:
		fbin = open(direct+"/histograms/bin_par"+str(ipar)+".data", 'w')
	fbin.write("# Tangential velocity (km/s)\n")
	fhisto = open(direct+"/histograms/emiss_part_shock_par"+str(ipar)+".data", 'w')
	fhisto.write("# Part of the emissivity which depends only on the shock.\n")
	ftemp = open(direct+"/histograms/temp_emiss_par"+str(ipar)+".data", 'w')
	ftemp.write("# Temperature of the cell where the emission is created.\n")
	ftemp.write("# Each line of a bin of the line profile.\n")
	for i in range(nbr_bin_profile):
		ftemp.write(' '.join(map(str, temp_emiss[i]))+"\n")
		fhisto.write(' '.join(map(str, emiss_part1[i]))+"\n")
		fbin.write(str(-vvv[i]/1.E5)+"\n")
	ftemp.close()
	fhisto.close()
	fbin.close()


	return crashing

