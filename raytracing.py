#!usr/bin/python
"""
****************************************        
***   Program for the ray tracing    ***         
****************************************
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

Run: from raytracing import RT
     ray_tracing=RT(ndeltaZZ, xx, yy, zz, deltaZZ, deltaX, deltaY, deltaZ, x_vec, y_vec, z_vec, d, cote_cd, kT_all, sigma_all, Mdot1, Mdot2, R1, R2, Teff_1, Teff_2, dwind_prim, dwind_sec, dir_obs, mp, zt1, pt1, yt1, zt2, pt2, yt2, treecdw1, treecdw2, tree, treevec)
"""  
# Subroutine: optical depth wind
# ==============================
def tauX(zz,p,tau0,R):
	import math
	pi=math.pi
	test=math.atan(p/zz)
	if (test < 0):
		test=pi-abs(test)
	alp=math.tan((pi-test)/2.)
	tau=0.
	if (abs(alp) != 0):
		if (abs(p) < R):
			x1=R/p+math.sqrt(R**2/p**2-1.)
			x2=R/p-math.sqrt(R**2/p**2-1.)
			tau=(tau0/math.sqrt(1.-p**2/R**2))*math.log(max([abs(alp-x2),abs(alp-x1)])/min([abs(alp-x2),abs(alp-x1)]))
		if (abs(p) == R):
			tau=2.*tau0/(alp-1.)
		if (abs(p) > R):
			terme=0.5*pi-math.atan((alp*p/R-1.)/math.sqrt(p**2/R**2-1.))
			if (terme < 0.):
				terme=2.*pi+terme
			tau=(2.*tau0/math.sqrt(p**2/R**2-1.))*terme
	return abs(tau)

def RT(ndeltaZZ, xx, yy, zz, deltaZZ, deltaX, deltaY, deltaZ, x_vec, y_vec, z_vec, d, cote_cd, kT_all, sigma_all, Mdot1, Mdot2, R1, R2, Teff_1, Teff_2, tree_prim_wind, tree_sec_wind, dwind_prim, dwind_sec, dir_obs, mp, zt1, pt1, yt1, zt2, pt2, yt2, treecdw1, treecdw2, tree, treevec, xmin1, xmax2, ymin1, ymax2):
	import math
	import numpy as np
	from scipy.spatial import KDTree

	ray_tracing=[]
	pi=math.pi
	range_points=np.arange(int(ndeltaZZ))+1.
	ray=np.array([xx+range_points*deltaX,yy+range_points*deltaY,zz+range_points*deltaZ]).T
	incr_zz=0
	for _ in xrange(int(ndeltaZZ)):
		xx=xx+deltaX #cm
		yy=yy+deltaY
		zz=zz+deltaZ
		test_dist=yy
		min1=ymin1
		max2=ymax2
		if (deltaX>deltaY):
			test_dist=xx
			min1=xmin1
			max2=xmax2
		if (test_dist > min1 and test_dist < max2):
			ind_close = tree.query([[xx/d,yy/d,zz/d]])[1][0]
			distvec, ind_closevec = treevec.query([[xx/d,yy/d,zz/d]])
			distvec=distvec[0]
			ind_closevec=ind_closevec[0]

			# Which side of the shock?
			if (cote_cd[ind_close] == 1):
				cote_choc = 1
				Mdot=Mdot1
				radius=R1
				distcdw, ind_closecdw = treecdw1.query([[x_vec[ind_closevec]/d, y_vec[ind_closevec]/d, z_vec[ind_closevec]/d]])
				distcdw=distcdw[0]
				inout='in'
				if (distcdw <= distvec):
					inout='out'
					rien, ind_close2 = tree_prim_wind.query([[xx/d,math.sqrt(yy**2+zz**2)/d]])
					ind_close2=ind_close2[0]
					vwind=np.sqrt(dwind_prim.iloc[ind_close2]['ux']**2+dwind_prim.iloc[ind_close2]['uy']**2)
					T0=Teff_1*8.61732e-8
					Twind=0.4*T0+0.6*T0*(radius/np.sqrt(xx**2+yy**2+zz**2))**1.9
			else:
				cote_choc = 2
				Mdot=Mdot2
				radius=R2
				distcdw, ind_closecdw = treecdw2.query([[x_vec[ind_closevec]/d, y_vec[ind_closevec]/d, z_vec[ind_closevec]/d]])
				distcdw=distcdw[0]
				inout='in'
				if (distcdw <= distvec):
					inout='out'
					rien, ind_close2 = tree_sec_wind.query([[(d-xx)/d,math.sqrt(yy**2+zz**2)/d]])
					ind_close2=ind_close2[0]
					vwind=np.sqrt(dwind_sec.iloc[ind_close2]['ux']**2+dwind_sec.iloc[ind_close2]['uy']**2)
					T0=Teff_2*8.61732e-8
					Twind=0.4*T0+0.6*T0*(radius/np.sqrt((d-xx)**2+yy**2+zz**2))**1.9

			if (inout == 'in'):
				#Point in the shock
				rdeltan=np.array([xx,yy,zz])/np.linalg.norm(np.array([xx,yy,zz]))
				ray_tracing.append(kT_all[ind_close]) #keV
				ray_tracing.append(abs(sigma_all[ind_close])/(abs(np.dot(rdeltan,dir_obs))*(1.3*mp*1000.)))
				incr_zz=incr_zz+1
			if (inout == 'out'):
				#Point in the wind
				if (cote_choc == 2):
					ray_tracing.append(Twind)
					append_val=Mdot*tauX(zt2+(incr_zz+1.)*deltaZZ, math.sqrt(pt2**2+yt2**2), 1., radius)/(4.*pi*vwind*radius*1.3*mp*1000.)
					if (vwind<100.):
						append_val=0.
					ray_tracing.append(append_val)
					
					dist_reste,rien = treecdw2.query(ray[int(incr_zz+1.):,:])
					index_neg=np.where(((dist_reste[1:]-dist_reste[0:-1]) < 0.))[0]
					if (len(index_neg) == 0.):
						Twind_reste=0.4*T0+0.6*T0*(radius/np.sqrt((d-ray[int(incr_zz+1.):,0])**2+ray[int(incr_zz+1.):,1]**2+ray[int(incr_zz+1.):,2]**2))**1.9 #keV
						for irest in range(len(Twind_reste)/30):
							ray_tracing.append(Twind_reste[irest*30])
							rien, ind_close2 = tree_sec_wind.query([[(d-ray[int(incr_zz+1.+irest*30),0])/d,math.sqrt(ray[int(incr_zz+1.+irest*30),1]**2+ray[int(incr_zz+1.+irest*30),2]**2)/d]])
							ind_close2=ind_close2[0]
							vwind_rest=np.sqrt(dwind_sec.iloc[ind_close2]['ux']**2+dwind_sec.iloc[ind_close2]['uy']**2)
							ray_tracing.append(Mdot*tauX(zt2+(incr_zz+irest*30.+2.)*deltaZZ, math.sqrt(pt2**2+yt2**2), 1., radius)/(4.*pi*vwind_rest*radius*1.3*mp*1000.))
						break
					index_neg2=np.where(((index_neg[1:]-index_neg[0:-1]) >1))[0]
					if (len(index_neg2) == 0. and (index_neg[-1] >= len(ray[int(incr_zz):,0])-1)):
						Twind_reste=0.4*T0+0.6*T0*(radius/np.sqrt((d-ray[int(incr_zz+1.):,0])**2+ray[int(incr_zz+1.):,1]**2+ray[int(incr_zz+1.):,2]**2))**1.9 #keV
						for irest in range(len(Twind_reste)/30):
							ray_tracing.append(Twind_reste[irest*30])
							rien, ind_close2 = tree_sec_wind.query([[(d-ray[int(incr_zz+1.+irest*30),0])/d,math.sqrt(ray[int(incr_zz+1.+irest*30),1]**2+ray[int(incr_zz+1.+irest*30),2]**2)/d]])
							ind_close2=ind_close2[0]
							vwind_rest=np.sqrt(dwind_sec.iloc[ind_close2]['ux']**2+dwind_sec.iloc[ind_close2]['uy']**2)
							ray_tracing.append(Mdot*tauX(zt2+(incr_zz+irest*30.+2.)*deltaZZ, math.sqrt(pt2**2+yt2**2), 1., radius)/(4.*pi*vwind_rest*radius*1.3*mp*1000.))
						break
					if (len(index_neg2) == 0.):
						index_reprendre=index_neg[-1]
					else:
						index_reprendre=index_neg[index_neg2[0]]

					Twind_reste=0.4*T0+0.6*T0*(radius/np.sqrt((d-ray[int(incr_zz+1.):int(index_reprendre+incr_zz),0])**2 +ray[int(incr_zz+1.):int(index_reprendre+incr_zz),1]**2 +ray[int(incr_zz+1.):int(index_reprendre+incr_zz),2]**2))**1.9 #keV
					for irest in range(len(Twind_reste)/30):
						ray_tracing.append(Twind_reste[irest*30])
						rien, ind_close2 = tree_sec_wind.query([[(d-ray[int(incr_zz+1.+irest*30),0])/d,math.sqrt(ray[int(incr_zz+1.+irest*30),1]**2+ray[int(incr_zz+1.+irest*30),2]**2)/d]])
						ind_close2=ind_close2[0]
						vwind_rest=np.sqrt(dwind_sec.iloc[ind_close2]['ux']**2+dwind_sec.iloc[ind_close2]['uy']**2)
						ray_tracing.append(Mdot*tauX(zt2+(incr_zz+irest*30.+2.)*deltaZZ, math.sqrt(pt2**2+yt2**2), 1., radius)/(4.*pi*vwind_rest*radius*1.3*mp*1000.))
					incr_zz=index_reprendre+incr_zz+1

				if (cote_choc == 1):
					ray_tracing.append(Twind)
					append_val=Mdot*tauX(zt1+(incr_zz+1.)*deltaZZ, math.sqrt(pt1**2+yt1**2), 1., radius)/(4.*pi*vwind*radius*1.3*mp*1000.)
					if (vwind<100.):
						append_val=0.
					ray_tracing.append(append_val)
					dist_reste,rien = treecdw1.query(ray[int(incr_zz):,:])
					index_neg=np.where(((dist_reste[1:]-dist_reste[0:-1]) < 0.))[0]
					if (len(index_neg) == 0.):
						Twind_reste=0.4*T0+0.6*T0*(radius/np.sqrt(ray[int(incr_zz+1.):,0]**2+ray[int(incr_zz+1.):,1]**2+ray[int(incr_zz+1.):,2]**2))**1.9 #keV
						for irest in range(len(Twind_reste)/30):
							ray_tracing.append(Twind_reste[irest*30])
							rien, ind_close2 = tree_prim_wind.query([[ray[int(incr_zz+1.+irest*30),0]/d,math.sqrt(ray[int(incr_zz+1.+irest*30),1]**2+ray[int(incr_zz+1.+irest*30),2]**2)/d]])
							ind_close2=ind_close2[0]
							vwind_rest=np.sqrt(dwind_prim.iloc[ind_close2]['ux']**2+dwind_prim.iloc[ind_close2]['uy']**2)
							ray_tracing.append(Mdot*tauX(zt1+(incr_zz+irest*30.+2.)*deltaZZ, math.sqrt(pt1**2+yt1**2), 1., radius)/(4.*pi*vwind_rest*radius*1.3*mp*1000.))
						break
					index_neg2=np.where(((index_neg[1:]-index_neg[0:-1]) >1))[0]
					if (len(index_neg2) == 0. and (index_neg[-1] >= len(ray[int(incr_zz):,0])-1)):
						Twind_reste=0.4*T0+0.6*T0*(radius/np.sqrt(ray[int(incr_zz+1.):,0]**2+ray[int(incr_zz+1.):,1]**2+ray[int(incr_zz+1.):,2]**2))**1.9 #keV
						for irest in range(len(Twind_reste)/30):
							ray_tracing.append(Twind_reste[irest*30])
							rien, ind_close2 = tree_prim_wind.query([[ray[int(incr_zz+1.+irest*30),0]/d,math.sqrt(ray[int(incr_zz+1.+irest*30),1]**2+ray[int(incr_zz+1.+irest*30),2]**2)/d]])
							ind_close2=ind_close2[0]
							vwind_rest=np.sqrt(dwind_prim.iloc[ind_close2]['ux']**2+dwind_prim.iloc[ind_close2]['uy']**2)
							ray_tracing.append(Mdot*tauX(zt1+(incr_zz+irest*30.+2.)*deltaZZ, math.sqrt(pt1**2+yt1**2), 1., radius)/(4.*pi*vwind_rest*radius*1.3*mp*1000.))
						break
					if (len(index_neg2) == 0.):
						index_reprendre=index_neg[-1]
					else:
						index_reprendre=index_neg[index_neg2[0]]

					Twind_reste=0.4*T0+0.6*T0*(radius/np.sqrt(ray[int(incr_zz+1.):int (index_reprendre+incr_zz),0]**2 +ray[int(incr_zz+1.):int(index_reprendre+incr_zz),1]**2 +ray[int(incr_zz+1.):int(index_reprendre+incr_zz),2]**2))**1.9 #keV
					for irest in range(len(Twind_reste)/30):
						ray_tracing.append(Twind_reste[irest*30])
						rien, ind_close2 = tree_prim_wind.query([[ray[int(incr_zz+1.+irest*30),0]/d,math.sqrt(ray[int(incr_zz+1.+irest*30),1]**2+ray[int(incr_zz+1.+irest*30),2]**2)/d]])
						ind_close2=ind_close2[0]
						vwind_rest=np.sqrt(dwind_prim.iloc[ind_close2]['ux']**2+dwind_prim.iloc[ind_close2]['uy']**2)
						ray_tracing.append(Mdot*tauX(zt1+(incr_zz+irest*30.+2.)*deltaZZ, math.sqrt(pt1**2+yt1**2), 1., radius)/(4.*pi*vwind_rest*radius*1.3*mp*1000.))
					incr_zz=index_reprendre+incr_zz+1
		if (incr_zz >= ndeltaZZ):
			break

	return ray_tracing
