#!usr/bin/python
"""
***************************************
***   Program for the plots of the  ***
***           stellar winds         ***         
***************************************
Run: python ~/PATH/plots_winds.py

# Parameters 
# ==========

# Versions
# ========
v1 - 22/4/2020
"""        

# Subroutine: to plot stellar surface
# ===================================
def plot_star(centerx,centery,radius,x):
	import numpy as np
	return np.sqrt(radius**2-(x-centerx)**2)+centery

# Main
# ====
def plots(R1,R2,vinf1,d,nbr_star,istar,ipar,direct):
	import numpy as np
	import os
	import sys
	import matplotlib
	import matplotlib.pyplot as plt  
	import pandas as pd
	import matplotlib.mlab as mlab

	############
	### Plot ###
	############
	# Read stellar winds
	wind=direct+"/wind_star"+str(istar+1)+"_param"+str(int(ipar))+".h5"
	dwind=pd.read_hdf(wind).reset_index(drop=True)
	ux=dwind.iloc[:]['ux']
	uy=dwind.iloc[:]['uy']
	x_vec=dwind.iloc[:]['x']
	y_vec=dwind.iloc[:]['y']

	index=np.where((np.array(x_vec[1:])-np.array(x_vec[:-1]) <0))[0]
	nbr_dans_ur=index[0]+1
	nbr_dans_ut=len(np.array(y_vec))/nbr_dans_ur
	je_veux_dans_ur=12
	je_veux_dans_ut=20
	indice=[]
	for i in range(je_veux_dans_ut):
		indice=np.concatenate((indice,np.floor(np.arange(je_veux_dans_ur)*float(nbr_dans_ur)/float(je_veux_dans_ur))+float(i)*round(nbr_dans_ut/float(je_veux_dans_ut))*nbr_dans_ur))
		if i == 0:
			new_nbr_ur=len(indice)
	indice=indice[np.where((indice < len(x_vec)))[0]]
	x_vec=np.array(x_vec[indice.astype(int)])
	y_vec=np.array(y_vec[indice.astype(int)])
	ux=np.array(ux[indice.astype(int)])
	uy=np.array(uy[indice.astype(int)])

	xnew,ynew = np.mgrid[min(x_vec):max(x_vec):70j,min(y_vec):max(y_vec):70j]
	ZI = mlab.griddata(x_vec,y_vec,np.sqrt(ux**2+uy**2)/1.e5,xnew,ynew,interp='linear')

	circle1 = plt.Circle((0, 0), R1/d, color='white')
	circle2 = plt.Circle((1, 0), R2/d, color='white')
	p = [(1, 0), (1, R2/d), (2, 2.*R2/d), (2, 0)]
	xlim1=-1.6
	xlim2=1.6
	xstar1=np.linspace(-R1/d,R1/d,40)
	ystar1=plot_star(0.0,0.0,R1/d,xstar1)
	xstar2=np.linspace(-R2/d,R2/d,40)
	ystar2=plot_star(0.0,0.0,R2/d,xstar2)
	xstar2=xstar2+1.

	if (istar == 1):
		xstar1=-xstar1+1.
		xstar2=-xstar2+1.
		x_vec=-x_vec+1.
		xnew=-xnew+1.
		p = [(0, 0), (0, R2/d), (-2, 3.*R2/d), (-2, 0)]
		circle1 = plt.Circle((1, 0), R1/d, color='white')
		circle2 = plt.Circle((0, 0), R2/d, color='white')
		xlim1=-0.6
		xlim2=2.6
		ux=-ux

	plt.figure(1,figsize=(11,11.*0.41))
	fig = plt.gcf()
	ax = fig.gca()
	plt.pcolor(xnew,ynew,ZI)
	cbar=plt.colorbar()
	cbar.set_label(r'Wind velocity [km/s]')
	plt.clim(0,vinf1/1.e5)
	plt.quiver(x_vec,y_vec, ux, uy, scale=np.max(np.sqrt(ux**2+uy**2))/0.04, units='xy', width=0.002, headlength=4, headaxislength=4, headwidth=7.)
	poly=plt.Polygon(p,color='white')
	ax.add_patch(poly)
	ax.add_artist(circle1)
	ax.add_artist(circle2)
	plt.plot(xstar1,ystar1,c='black')
	plt.plot(xstar2,ystar2,c='black')
	plt.xlim(xlim1,xlim2)
	plt.ylim(0,1.8)
	plt.ylabel("x/D")
	plt.xlabel("y/D")
	plt.savefig(direct+"/wind_star"+str(istar+1)+"_param"+str(ipar)+".pdf")
	if (nbr_star==1):
		plt.savefig(direct+"/wind_star2_param"+str(ipar)+".pdf")
	plt.close()

	return None

