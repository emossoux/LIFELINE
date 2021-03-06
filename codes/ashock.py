#!usr/bin/python
"""
***************************************************        
***   Program for the computation of the shock  ***
***     in adiabatic colliding wind binaries    ***    
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

Run: from ashock import ashock
     rien=ashock(Mdot1,Mdot2,vinf1,vinf2, mass_ratio,ex,omega, a, Per, R1, R2, Teff1, Teff2, nbr_points_shock_2D, nbr_points_shock_3D, nbr_points_width, nbr_bin_profile, direct, cut_lim, wind_prim, wind_sec, phase, incl,ipar, mu)

Use of the Canto et al. (1996) formalism

# Parameters 
# ==========
Mdot - Mass loss rate [Msun/yr]
vinf - terminal velocity of the winds [cm/s]
ex - eccentricity
omega - argument of the periapsis [radian]
a - semi-major axis [Rsun]
R - stellar radius [Rsun]
Teff - Effective temperature of the star [K]
nbr_points_shock_2D - number of points to create the 2D shape of the shock
nbr_points_shock_3D - number of points to create the 3D shape of the shock
nstep_width - number of points to discretize the inside of the shock
nbr_bin_profile - number of points in the line profile
direct - directory where saving files
cut_lim - Limit of the tangential velocity to see the emission
wind_prim - h5 file containing the wind distribution of the primary 
wind_sec - h5 file containing the wind distribution of the secondary 
phase - phase of the system (not corrected from omega)
incl - inclination of the orbit
ipar - index of the system parameter
mu - Mean molecular weight

# Versions
# ========
v1 - 16/03/18
"""      

# Subroutine: compute eta
# =======================
def comp_eta(x, y, R1, R2, Mdot1, Mdot2, d, tree_prim_wind, tree_sec_wind, dwind_prim, dwind_sec, vinf1, vinf2):
	import numpy as np
	vv1=100.e5
	vv2=100.e5
	if (x**2+y**2 < R1**2):
		vv1=100.e5
		# closest point
		rien, ind_close = tree_sec_wind.query([[1.-R1/d,y/d]])
		ind_close=ind_close[0]
		vv2=np.sqrt(dwind_sec.iloc[ind_close]['ux']**2+dwind_sec.iloc[ind_close]['uy']**2)
		if (np.isfinite(vv2) == False):
			vv2=0.
		if (vv2 > vinf2):
			vv2=vinf2
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
		if (vv1 > vinf1):
			vv1=vinf1
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
		if (np.isfinite(vv2) == False):
			vv2=0.
		if (vv1 > vinf1):
			vv1=vinf1
		if (vv2 > vinf2):
			vv2=vinf2
	return (Mdot2*vv2)/(Mdot1*vv1)

# Subroutine: find x0
# ===================
def find_x0(x, d, R1, R2, Mdot1, Mdot2, tree_prim_wind, tree_sec_wind, dwind_prim, dwind_sec, vinf1, vinf2):
	import math
	eta=comp_eta(x, 0., R1, R2, Mdot1, Mdot2, d, tree_prim_wind, tree_sec_wind, dwind_prim, dwind_sec, vinf1, vinf2)
	if (eta == 0.):
		return 0.
	else:
		return x-d/(1.+math.sqrt(eta))

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
def ashock(liste_param):
	from scipy.spatial import KDTree
	from scipy.optimize import fsolve  
	from scipy import spatial
	import numpy as np
	import math
	from constantes import constante
	import matplotlib.pyplot as plt
	import matplotlib 
	import os
	import pandas as pd
	from mpl_toolkits.mplot3d import Axes3D
	from raytracing import RT
	from astropy.io import ascii
	import matplotlib.tri as tri

	Mdot1,Mdot2,vinf1,vinf2, mass_ratio,ex, omega, a, per,R1,R2,Teff_1, Teff_2, nbr_points_shock_2D, nbr_points_shock_3D, nstep_width, nbr_bin_profile, direct, cut_lim, wind_prim, wind_sec, phase, incl,ipar, mu, M_conj=liste_param
        
	pi=math.pi
	dwind_prim=pd.read_hdf(wind_prim).reset_index(drop=True)
	dwind_sec=pd.read_hdf(wind_sec).reset_index(drop=True)
	tree_prim_wind=KDTree(np.reshape(np.concatenate((dwind_prim.iloc[:]['x'],dwind_prim.iloc[:]['y'])),(-1,int(len(dwind_prim.iloc[:]['x'])))).T,leafsize=1000)
	tree_sec_wind=KDTree(np.reshape(np.concatenate((dwind_sec.iloc[:]['x'],dwind_sec.iloc[:]['y'])),(-1,int(len(dwind_sec.iloc[:]['x'])))).T,leafsize=1000)
                                            
	beta=(Mdot1*vinf1)/(Mdot2*vinf2)

	Rsun=constante('Rsun') # solar radius in cm
	Msun=constante('Msun') # solar mass in g
	kb=constante('kb') # Boltzmann constant in erg/K
	mp=constante('mp')*1000. # proton mass in g
	Mdot1=Mdot1*Msun/(365.24*86400.) #g/s
	Mdot2=Mdot2*Msun/(365.24*86400.) #g/s
	R1=R1*Rsun #cm
	R2=R2*Rsun #cm
	a=a*Rsun #cm

	vmax=max([vinf1,vinf2])
	vec=np.arange(nbr_bin_profile, dtype=np.float)+1.
	vvv=(vec-(math.floor(float(nbr_bin_profile)/2.)+1.))*vmax/math.floor(float(nbr_bin_profile)/2.)

	# Use of the Canto's formalism
	# ===============================
	# Equation 28 of Canto96, theta_inf
	thetainf=pi/2.
	if (beta != 1.):
		asymp=pi/(1.-beta)
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

	M=phase
	E=M
	dE=(M-E+ex*math.sin(E))/(1.-ex*math.cos(E))
	while (abs(dE) >= 1.0E-8):
		E=E+dE
		dE=(M-E+ex*math.sin(E))/(1.-ex*math.cos(E))
	coi=(math.cos(E)-ex)/(1.-ex*math.cos(E))
	sii=math.sqrt(1.-ex*ex)*math.sin(E)/(1.-ex*math.cos(E))
	# True anomaly 
	phase_true=math.atan2(sii,coi)
	if (phase_true == 2.*pi):
		phase_true = 0.
	d=a*(1.-ex**2)/(1.+ex*math.cos(phase_true)) #cm
	ouverture_min=math.atan(R1/d)

	xstag=fsolve(find_x0, (R1+d-R2)/2., args=(d, R1, R2, Mdot1, Mdot2, tree_prim_wind,tree_sec_wind, dwind_prim, dwind_sec, vinf1, vinf2)) #cm 
	xstag=xstag[0]

	fxstag = open(direct+"/histograms/xstag_par"+str(ipar), 'w')
	fxstag.write(str(xstag)+" cm\n")
	fxstag.write(str(xstag/R1)+" R1\n")
	fxstag.write(str(xstag/d)+" d\n")
	fxstag.close()
	xstag=max([xstag,R1])

	phase=phase-M_conj
	if (phase < 0):
		phase = 2.*pi+phase

	x_cdw1, y_cdw1, z_cdw1, x_cdw2, y_cdw2, z_cdw2, x_cdw1g, y_cdw1g, z_cdw1g, x_cdw2g, y_cdw2g, z_cdw2g, kT_all, rho_all, sigma_all=([] for _ in range(15))
	x_cd, y_cd, z_cd, vol_cd, vrad_cd, x_skew, y_skew, z_skew, xx_all, yy_all, zz_all, cote_cd, vt_vec=([] for _ in range(13))
	vt1_vec, vp1_vec, vt2_vec, vp2_vec, slope_vec1, slope_vec2 =([] for _ in range(6))

	phi_vec=np.arange(nbr_points_shock_3D)*360.*pi/(180.*nbr_points_shock_3D)
	dphi=360.*pi/(180.*nbr_points_shock_3D)

	xstar1=np.linspace(-R1/d,R1/d,40)
	ystar1=plot_star(0.0,0.0,R1/d,xstar1)
	xstar2=np.linspace(-R2/d,R2/d,40)
	ystar2=plot_star(0.0,0.0,R2/d,xstar2)
	xstar2=xstar2+1.

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
	ind_half=np.where((dwind_sec['x'] > 0.5))[0]
	dwind_sec_half=dwind_sec.iloc[ind_half]
	x_old=xstag
	y_old=0.
	x_vec=[]
	eta_vec=[]
	for y_vechere in y_vec:
		r1=math.sqrt(x_old**2+y_old**2)
		r2=math.sqrt((d-x_old)**2+y_old**2)
		eta=comp_eta(x_old, y_old, R1, R2, Mdot1, Mdot2, d, tree_prim_wind,tree_sec_wind, dwind_prim, dwind_sec, vinf1, vinf2)
		if (eta != 1):
			pente=y_vechere/(x_old-(d*r1**2*math.sqrt(eta))/(r1**2*math.sqrt(eta)+r2**2))
			pente=1./pente
			x_old=pente*(y_vechere-y_old)+x_old
		y_old=y_vechere
		x_vec.append(x_old)
		eta_vec.append(eta)
	eta_vec=np.array(eta_vec)
	x_vec=np.array(x_vec)/d
	y_vec=y_vec/d

	xaxis=[1.,0.,0.]
	lmax=0.
	for icut in range(nbr_points_shock_3D):
		if (float(icut+1)%100. == 0.):
			sentence="* "+str(icut+1)+"/"+str(nbr_points_shock_3D)+" *"
			print(sentence)
		xcut=x_vec*d #cm
		ycut=y_vec*d*np.cos(phi_vec[icut]) #cm
		zcut=y_vec*d*np.sin(phi_vec[icut]) #cm

		for cut in range(int(len(x_vec))):
			if (xcut[cut]<0. and math.atan(math.sqrt(ycut[cut]**2+zcut[cut]**2)/(d-xcut[cut]))<ouverture_min):
				Rprime=math.tan(ouverture_min)*(d-xcut[cut])
				angle=math.acos(ycut[cut]/math.sqrt(ycut[cut]**2+zcut[cut]**2))
				ycut[cut]=abs(Rprime*math.cos(angle))*float(np.sign(ycut[cut]))
				zcut[cut]=abs(Rprime*math.sin(angle))*float(np.sign(zcut[cut]))

			thetah=np.arccos(np.dot([xcut[cut],ycut[cut],zcut[cut]],xaxis)/(np.linalg.norm([xcut[cut],ycut[cut],zcut[cut]], axis=0)*np.linalg.norm(xaxis,axis=0)))
			lat=angleyz(ycut[cut],zcut[cut])
			the1_vech=pi-np.arccos(np.dot([xcut[cut]-d,ycut[cut],zcut[cut]],xaxis)/(np.linalg.norm([xcut[cut]-d,ycut[cut],zcut[cut]], axis=0)*np.linalg.norm(xaxis,axis=0)))

			thetahm1=0.
			the1_vechm1=0.
			slope_vech=pi/2.
			if (cut >= 1):
				thetahm1=np.arccos(np.dot([xcut[cut-1],ycut[cut-1],zcut[cut-1]],xaxis)/(np.linalg.norm([xcut[cut-1],ycut[cut-1],zcut[cut-1]], axis=0)*np.linalg.norm(xaxis,axis=0)))
				the1_vechm1=pi-np.arccos(np.dot([xcut[cut-1]-d,ycut[cut-1],zcut[cut-1]],xaxis)/(np.linalg.norm([xcut[cut-1]-d,ycut[cut-1],zcut[cut-1]], axis=0)*np.linalg.norm(xaxis,axis=0)))
				slope_vech=np.arccos(np.dot([xcut[cut]-xcut[cut-1],ycut[cut]-ycut[cut-1],zcut[cut]-zcut[cut-1]],xaxis)/(np.linalg.norm([xcut[cut]-xcut[cut-1],ycut[cut]-ycut[cut-1],zcut[cut]-zcut[cut-1]], axis=0)* np.linalg.norm(xaxis,axis=0)))
			dtheta=abs(thetah-thetahm1)
			dthe1=abs(the1_vech-the1_vechm1)

			r1_vech=np.sqrt(xcut[cut]**2+ycut[cut]**2+zcut[cut]**2) #cm
			r2_vech=np.sqrt((d-xcut[cut])**2+ycut[cut]**2+zcut[cut]**2)

			rien, ind_close = tree_prim_wind.query([[xcut[cut]/d,np.sqrt(ycut[cut]**2+zcut[cut]**2)/d]])
			ind_close=ind_close[0]
			vv1=np.sqrt(dwind_prim.iloc[ind_close]['ux']**2+dwind_prim.iloc[ind_close]['uy']**2)
			vv1=min([vinf1,vv1])
			vv1=max([100.e5,vv1])
			rien, ind_close = tree_sec_wind.query([[1.-xcut[cut]/d,np.sqrt(ycut[cut]**2+zcut[cut]**2)/d]])
			ind_close=ind_close[0]
			vv2=np.sqrt(dwind_sec.iloc[ind_close]['ux']**2+dwind_sec.iloc[ind_close]['uy']**2)
			vv2=min([vinf2,vv2])
			vv2=max([100.e5,vv2])

			# Equations 29 & 30
			sig0=Mdot1/(2.*pi*eta_vec[cut]*d*vv1) #g/cm^2
			alpha=vv2/vv1
			term1=eta_vec[cut]*(thetah-np.sin(thetah)*np.cos(thetah))
			term1=term1+(the1_vech-np.sin(the1_vech)*np.cos(the1_vech))
			term1=term1*term1
			term2=eta_vec[cut]*np.sin(thetah)*np.sin(thetah)
			term2=term2-np.sin(the1_vech)*np.sin(the1_vech)
			term2=term2*term2
			term3=eta_vec[cut]*(1.-np.cos(thetah))+alpha*(1.-np.cos(the1_vech))
			term4=np.sin(thetah+the1_vech)/np.sin(thetah)
			term4=(term4/np.sin(the1_vech))*term3*term3
			sigma_vec=sig0*term4/np.sqrt(term1+term2) #g/cm^2

			rho1=4.*Mdot1/(4.*pi*(xcut[cut]**2+ycut[cut]**2+zcut[cut]**2)*vv1) #g/cm^3
			rho2=4.*Mdot2/(4.*pi*((d-xcut[cut])**2+ycut[cut]**2+zcut[cut]**2)*vv2) #g/cm^3
			rho=(rho1+rho2)/2. #mean density of the mixed gas

			# Shock thickness 
			thickness1=sigma_vec/rho1 #cm
			thickness2=sigma_vec/rho2 #cm
			l0_1=thickness1
			l0_2=thickness2
			if (cut >= 1):
				if (r1_vech-l0_1<rl1):
					l0_1=l011
				if (r2_vech-l0_2<rl2):
					l0_2=l022
			lmax=max(lmax,l0_1,l0_2)

			l1_width=(np.arange(nstep_width))*((2.*nstep_width-1)*l0_1/(2.*nstep_width))/(nstep_width-1)+l0_1/(2.*nstep_width)
			l2_width=(np.arange(nstep_width))*((2.*nstep_width-1)*l0_2/(2.*nstep_width))/(nstep_width-1)+l0_2/(2.*nstep_width)

			increasevol=np.arange(nstep_width)
			vol_cd=np.concatenate((vol_cd,l0_1*(r1_vech-l0_1*increasevol/nstep_width)**2*math.sin(pi/2.-thetah)*dtheta*dphi/(abs(math.sin(slope_vech-thetah))*nstep_width)))
			vol_cd=np.concatenate((vol_cd,l0_2*(r2_vech-l0_2*increasevol/nstep_width)**2*math.sin(pi/2.-the1_vech)*dthe1*dphi/(abs(math.sin(slope_vech-the1_vech))*nstep_width)))
			x_cdw1g1=projx(r1_vech-l0_1,thetah)
			y_cdw1g1=projy(r1_vech-l0_1,thetah,lat)
			z_cdw1g1=projz(r1_vech-l0_1,thetah,lat)
			x_cdw2g1=projx(r2_vech-l0_2,pi-the1_vech)+d
			y_cdw2g1=projy(r2_vech-l0_2,pi-the1_vech,lat)
			z_cdw2g1=projz(r2_vech-l0_2,pi-the1_vech,lat)
			x_cdw1g1a=projx(r1_vech-l1_width,thetah)
			y_cdw1g1a=projy(r1_vech-l1_width,thetah,lat)
			z_cdw1g1a=projz(r1_vech-l1_width,thetah,lat)
			x_cdw2g1a=projx(r2_vech-l2_width,pi-the1_vech)+d
			y_cdw2g1a=projy(r2_vech-l2_width,pi-the1_vech,lat)
			z_cdw2g1a=projz(r2_vech-l2_width,pi-the1_vech,lat)

			vt1=vv1*math.cos(pi/2.-thetah)
			vt2=vv2*math.cos(pi/2.-the1_vech)
			if (cut >= 1):
				slope_vech=np.arccos(np.dot([x_cdw1g1-x_cdw1g11,y_cdw1g1-y_cdw1g11,z_cdw1g1-z_cdw1g11],xaxis)/(np.linalg.norm([x_cdw1g1-x_cdw1g11,y_cdw1g1-y_cdw1g11,z_cdw1g1-z_cdw1g11], axis=0)* np.linalg.norm(xaxis,axis=0)))
				vt1=vv1*math.cos(slope_vech-thetah)
				slope_vech=np.arccos(np.dot([x_cdw2g1-x_cdw2g11,y_cdw2g1-y_cdw2g11,z_cdw2g1-z_cdw2g11],xaxis)/(np.linalg.norm([x_cdw2g1-x_cdw2g11,y_cdw2g1-y_cdw2g11,z_cdw1g1-z_cdw2g11], axis=0)* np.linalg.norm(xaxis,axis=0)))
				vt2=vv2*math.cos(pi-slope_vech-the1_vech)
			x_cdw1g11=x_cdw1g1
			y_cdw1g11=y_cdw1g1
			z_cdw1g11=z_cdw1g1
			x_cdw2g11=x_cdw2g1
			y_cdw2g11=y_cdw2g1
			z_cdw2g11=z_cdw2g1
			rl1=r1_vech-l0_1
			rl2=r2_vech-l0_2
			l011=l0_1
			l022=l0_2 


			vp1=math.sqrt(vv1**2-vt1**2)
			vp2=math.sqrt(vv2**2-vt2**2)

			kT1=np.arange(nstep_width)*0.+8.61732e-8*3.*mp*vp1**2/(16.*kb) #keV
			kT2=np.arange(nstep_width)*0.+8.61732e-8*3.*mp*vp2**2/(16.*kb) #keV

			# named _skew by historical reasons by not skewing applied
			x_skew=np.concatenate((x_skew,xcut)) #cm
			y_skew=np.concatenate((y_skew,ycut))
			z_skew=np.concatenate((z_skew,zcut))

			# Shock surface
			x_cdw1g=np.concatenate((x_cdw1g,[x_cdw1g1]))
			y_cdw1g=np.concatenate((y_cdw1g,[y_cdw1g1]))
			z_cdw1g=np.concatenate((z_cdw1g,[z_cdw1g1]))
			x_cdw2g=np.concatenate((x_cdw2g,[x_cdw2g1]))
			y_cdw2g=np.concatenate((y_cdw2g,[y_cdw2g1]))
			z_cdw2g=np.concatenate((z_cdw2g,[z_cdw2g1]))

			# Overall points inside the shock
			xx_all=np.concatenate((xx_all,x_cdw1g1a)) #cm
			yy_all=np.concatenate((yy_all,y_cdw1g1a))
			zz_all=np.concatenate((zz_all,z_cdw1g1a))
			xx_all=np.concatenate((xx_all,x_cdw2g1a)) #cm
			yy_all=np.concatenate((yy_all,y_cdw2g1a))
			zz_all=np.concatenate((zz_all,z_cdw2g1a))

			slope_vech=pi/2.
			if (cut >= 1):
				slope_vech=np.arccos(np.dot([x_cdw1g[-1]-x_cdw1g[-2],y_cdw1g[-1]-y_cdw1g[-2],z_cdw1g[-1]-z_cdw1g[-2]],xaxis)/(np.linalg.norm([x_cdw1g[-1]-x_cdw1g[-2],y_cdw1g[-1]-y_cdw1g[-2],z_cdw1g[-1]-z_cdw1g[-2]], axis=0)* np.linalg.norm(xaxis,axis=0)))
			lat=angleyz(y_cdw1g[-1],z_cdw1g[-1])
			vx=projx(np.arange(nstep_width)*0.+vt1,slope_vech) #cm
			vy=projy(np.arange(nstep_width)*0.+vt1,slope_vech,lat)
			vz=projz(np.arange(nstep_width)*0.+vt1,slope_vech,lat)
			vrad_cd=np.concatenate((vrad_cd,vx*math.sin(incl)*math.cos(phase)+vy*math.cos(incl)-vz*math.sin(incl)*math.sin(phase)))

			vx=projx(np.arange(nstep_width)*0.+vt2,slope_vech) #cm
			vy=projy(np.arange(nstep_width)*0.+vt2,slope_vech,lat)
			vz=projz(np.arange(nstep_width)*0.+vt2,slope_vech,lat)
			vrad_cd=np.concatenate((vrad_cd,vx*math.sin(incl)*math.cos(phase)+vy*math.cos(incl)-vz*math.sin(incl)*math.sin(phase)))

			kT_all=np.concatenate((kT_all,kT1))
			cote_cd=np.concatenate((cote_cd,1+0.*np.arange(nstep_width)))
			rho_all=np.concatenate((rho_all,np.arange(nstep_width)*0.+rho1))
			sigma_all=np.concatenate((sigma_all,np.arange(nstep_width)*0.+sigma_vec))

			kT_all=np.concatenate((kT_all,kT2))
			cote_cd=np.concatenate((cote_cd,2+0.*np.arange(nstep_width)))
			rho_all=np.concatenate((rho_all,np.arange(nstep_width)*0.+rho2))
			sigma_all=np.concatenate((sigma_all,np.arange(nstep_width)*0.+sigma_vec))

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

	circle1 = plt.Circle((0, 0), R1/d, color='white')
	circle2 = plt.Circle((1, 0), R2/d, color='white')
	factor=1.e-13
	while (max(rho_layer)/factor < 1):
		factor = factor/10.
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
	deltaX=deltaZZ*math.cos(phase)*math.sin(incl) #cm
	deltaY=deltaZZ*math.cos(incl)
	deltaZ=deltaZZ*math.sin(phase)*math.sin(incl)
	dir_obs=np.array([deltaX,deltaY,deltaZ])/np.linalg.norm(np.array([deltaX,deltaY,deltaZ]))

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
	print("*** Ray tracing ***")
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
		vrad=vrad_cd[i]

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

	return None

