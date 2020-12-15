#!usr/bin/python
"""
********************************************   
***  Program for renaming the files and  ***
***      implementing the database       ***        
********************************************
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

Run: python /PATH/rename_files.py

# Parameters 
# ==========
mode - In which mode did you run the line_profile.py code [1/2/3] 1 - The overall computation; 2 - The shock characteristics and the line profile; 3 - Only the line profile.
path_param - Complete path containing the winds[_set], histograms[_set], plots[_set] and [atom]_[ion][_set] directories.
ion_dir - Name of your [atom]_[ion] directory.
file_param - Name of your stellar parameter file.
old_ipar - Index of the binary system from your stellar parameter file to include in the database (begin at 0; can be a list, ex: 0 1 3 6).

# Versions
# ========
v1 - 09/03/2020
"""      
import os
import sys
import shutil
import glob
import math
import numpy as np
from constantes import constante

grav=constante('G')*1000. # gravity  in cm^3/g/s^2
Msun=constante('Msun') #g
Rsun=constante('Rsun') #cm
pi=math.pi
set_dir_name="_set"

# Code parameters
# ===============
mode = int(raw_input("In which mode did you run the line_profile.py code? [1/2/3] 1 - The overall computation; 2 - The shock characteristics and the line profile; 3 - Only the line profile.\n"))
path_param=raw_input("Complete path containing the winds[_set], histograms[_set], plots[_set] and [atom]_[ion][_set] directories:\n")
ion_dir=raw_input("Name of your [atom]_[ion] directory: ")
file_param=raw_input("Name of your stellar parameter file: ")
old_ipar=raw_input("Which is the index of the binary system from your stellar parameter file to include in the database (begin at 0; can be a list, ex: 0 1 3 6)? ").split()
if (mode != 3):
	same_par=raw_input("If you increased the index of the results files to not overwrite old files, write 'STOP'\n")
	if (same_par == 'STOP' or same_par == 'stop'):
		print("Add false lines in your stellar parameter file so that the indexes correspond to the name of the results files.")
		sys.exit()
old_ipar = [int(a) for a in old_ipar]
nbr_par=len(old_ipar)

# Read the stellar parameter file                                           
# ================================
fpar = open(path_param+"/"+file_param, 'r')
k=0
ipar_here=0
incl=[]
phase=[]
ligne=[]
ex=[]
omega=[]
a=[]
type1=[]
type2=[]
beta=[]
while 1:
    line=fpar.readline()
    if not line: break
    vec=line.split(' ')
    if (vec[0] != '#'):
	if (k == old_ipar[ipar_here]):
	    	vec=line.split(' ')
		incl.append(float(vec[-1]))
		phase.append(float(vec[-2]))
		omega.append(float(vec[-3]))
		ex.append(float(vec[-4]))
		a.append(float(vec[-5]))
		type1.append(vec[0])
		type2.append(vec[1])

		mass_sum=(float(vec[4])+float(vec[5]))*Msun #g
		Per=2.*pi*(a[-1]*Rsun)**1.5/(math.sqrt(grav*mass_sum)*86400.)
		mass_ratio=float(vec[4])/float(vec[5])
		vinf1=2.6*math.sqrt(2.*grav*float(vec[4])*Msun/(float(vec[6])*Rsun)) #vinf=2.6*vesc
		vinf2=2.6*math.sqrt(2.*grav*float(vec[5])*Msun/(float(vec[7])*Rsun))
		beta.append((float(vec[3])*vinf2)/(float(vec[2])*vinf1))

		ligne.append(' '.join(vec[:-5]))
		ipar_here=ipar_here+1
		if (ipar_here>=nbr_par):
			break
   	k=k+1
fpar.close()

# Read the set stellar parameter file                                           
# ===================================
fpar = open(path_param+"/stellar_params_set.txt", 'r')
type1_set, type2_set, d_set =([] for _ in range(3))
line=fpar.readline()	
nbr_already=0
while 1:
    line=fpar.readline()
    if not line: break
    nbr_already=nbr_already+1
    vec=line.split(' ')
    if (vec[0] != '#'):
	type1_set.append(vec[0])
	type2_set.append(vec[1])
	d_set.append(float(vec[10]))
type1_set=np.array(type1_set)
type2_set=np.array(type2_set)
d_set=np.array(d_set)
fpar.close()

phase_old=0.+np.array(phase)
# If mode=3, the new system index will correspond to an existing one
if (mode == 3):
	ipar2_vec=[]
	for k in range(nbr_par):
		M=phase[k]*2.*pi
		E=M
		dE=(M-E+ex[k]*math.sin(E))/(1.-ex[k]*math.cos(E))
		while (abs(dE) >= 1.0E-8):
			E=E+dE
			dE=(M-E+ex[k]*math.sin(E))/(1.-ex[k]*math.cos(E))
		coi=(math.cos(E)-ex[k])/(1.-ex[k]*math.cos(E))
		sii=math.sqrt(1.-ex[k]**2)*math.sin(E)/(1.-ex[k]*math.cos(E))
		phase_true=math.atan2(sii,coi)
		d=a[k]*(1.-ex[k]**2)/(1.+ex[k]*math.cos(phase_true)) #Rsun

		# Conjunction phase. If ex==0, phase_conj=0
		if (ex[k] > 0):
			E_conj=2.*math.atan(math.sqrt((1.-ex[k])/(1.+ex[k]))*math.tan(0.5*(0.5*pi-omega[k])))
			M_conj=E_conj-ex[k]*math.sin(E_conj)
			phase[k]=phase[k]-(M_conj+pi)
			if (beta[k] >= 1):
				phase[k]=phase[k]-M_conj
			if (phase[k]<0):
				phase[k]=2.*pi+phase[k]

		ipar2 = np.where((type2_set==type2[k]) & (type1_set==type1[k]))[0]
		if (len(ipar2) > 1):
			ipar3 = np.where(abs(d_set[ipar2]-d)==min(abs(d_set[ipar2]-d)))[0]
			ipar2 = ipar2[ipar3]
		if (len(ipar2) == 0):
			ipar2 = np.where((type1_set==type2[k]) & (type2_set==type1[k]))[0]
			if (len(ipar2) > 1):
				ipar3 = np.where(abs(d_set[ipar2]-d)==min(abs(d_set[ipar2]-d)))[0]
				ipar2 = ipar2[ipar3]
		if (isinstance(ipar2,np.ndarray)):
			ipar2=ipar2[0]
		ipar2_vec.append(ipar2)


# Rename and copy the files                                         
# =========================
if (mode == 1):
	direct=['winds','histograms','plots',ion_dir]
elif (mode == 2):
	direct=['histograms','plots',ion_dir]
elif (mode == 3):
	direct=[ion_dir]
else:
	print("Please enter 1, 2 or 3 for the mode")
	sys.exit()

if not os.path.exists(path_param+"/"+ion_dir+set_dir_name):
	os.makedirs(path_param+"/"+ion_dir+set_dir_name)

# If mode == 3, we just need to copy the files (which already have the name of the inclination and phase) to the _set directory)
if (mode == 3):
	for k in range(nbr_par):
		os.chdir(path_param+'/'+ion_dir)
		for files in glob.glob("*_par"+str(int(ipar2_vec[k]))+"*"):
			shutil.copyfile(files, '../'+ion_dir+set_dir_name+'/'+files)
else:
	for k in range(nbr_par):
		for i in direct:
			os.chdir(path_param+'/'+i)
			print ipar_new
			if (i == ion_dir or i == 'histograms'):
				for files in glob.glob("*_par"+str(int(old_ipar[k]))+"*"):
	    				vec=files.split('_par')
					if (len(vec) == 3):
						vec=[vec[0]+'_par'+vec[1],vec[2]]
	    				vec2=vec[1].split('_')
					if (len(vec2) == 1):
	    					vec2=vec[1].split('.')
						shutil.copyfile(files, '../'+i+set_dir_name+'/'+vec[0]+'_par'+str(int(ipar_new))+'_incl'+str(int(round(incl[k])))+'_phase'+str(int(round(phase[k]*100./(2.*pi))))+'.'+vec2[1])
					else:
						vec2='_'.join(vec2[1:])
						shutil.copyfile(files, '../'+i+set_dir_name+'/'+vec[0]+'_par'+str(int(ipar_new))+'_incl'+str(int(round(incl[k])))+'_phase'+str(int(round(phase[k]*100.)))+'_'+vec2)
			elif (i == 'plots'):
				for files in glob.glob("*_par"+str(int(old_ipar[k]))+"*"):
	    				vec=files.split('_par')
					shutil.copyfile(files, '../'+i+set_dir_name+'/'+vec[0]+'_par'+str(int(ipar_new))+'.pdf')
			elif (i == 'winds'):
				for files in glob.glob("*_param"+str(int(old_ipar[k]))+"*"):
	    				vec=files.split('_param')
					shutil.copyfile(files, '../'+i+set_dir_name+'/'+vec[0]+'_par'+str(int(ipar_new))+'.pdf')

		ipar_new=ipar_new+1
	

# Write the new binary in the set stellar parameter file                                           
# ======================================================
if (mode != 3):
	fpar = open(path_param+"/stellar_params_set.txt", 'a')
	for k in range(nbr_par):
		M=phase_old[k]
		E=M
		dE=(M-E+ex[k]*math.sin(E))/(1.-ex[k]*math.cos(E))
		while (abs(dE) >= 1.0E-8):
			E=E+dE
			dE=(M-E+ex[k]*math.sin(E))/(1.-ex[k]*math.cos(E))
		coi=(math.cos(E)-ex[k])/(1.-ex[k]*math.cos(E))
		sii=math.sqrt(1.-ex[k]**2)*math.sin(E)/(1.-ex[k]*math.cos(E))
		phase_true=math.atan2(sii,coi)
		d=a[k]*(1.-ex[k]**2)/(1.+ex[k]*math.cos(phase_true)) #Rsun
		fpar.write(ligne[k]+" "+str(round(d*100.)/100.)+" "+str(int(mode))+"\n")
	fpar.close()

	print("")
	if (len(old_ipar) == 1):
		print("Files renamed and saved in the *_set directories with the parameter index "+str(int(ipar_new-1)))
	else:
		print("Files renamed and saved in the *_set directories from the parameter index "+str(int(ipar_new-1-(len(old_ipar)-1)))+" to "+str(int(ipar_new-1)))
	print("Stellar parameters saved at the end of the 'stellar_params_set' file")
else:
	print("")
	print("Files saved in the "+ion_dir+set_dir_name+" directory")
