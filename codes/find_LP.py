#!usr/bin/python
"""
********************************************************  
***   Program to find the names of the line profile  ***
***        file of the closest computed system       ***         
******************************************************** 
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

Run: python script_root+find_LP.py binary_parameters_file

# Parameters 
# ==========
The binary_parameters_file is the file containing the stellar parameters. 
The user must run the code from the directory where this file and 
the list of already computed binaries (stellar_params_set.txt) are located.

# Versions
# ========
v1 - 30/03/2020
"""           
import numpy as np
import sys
import math 
import glob
import os
from constantes import constante

Rsun=constante('Rsun') #cm
Msun=constante('Msun') #g
grav=constante('G')*1000. # gravity  in cm^3/g/s^2

atom=input("Atom for which you want to retrieve the line profile (ex: Fe, Ca,...): ") 
ion=int(float(input("Ion of this atom (ex: 25, 14, ...): ")))
direct_set="/local/network/htdocs/Lifeline/"

# List of the atoms until 30
# ==========================
atom_list=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn']
roman=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII', 'XIV', 'XV', 'XVI', 'XVII', 'XVIII', 'XIX', 'XX', 'XXI', 'XXII', 'XXIII', 'XXIV', 'XXV', 'XXVI', 'XXVII','XXVIII','XXIX','XXX']
if atom not in atom_list:
	print("Please enter an atom between H and Zn (Take care of the capital letter)")
	sys.exit()
atomic_number=atom_list.index(atom)+1.
if (ion > atomic_number):
	print("Ion > Atomic number. Please enter a valid ion for this element.")
	sys.exit()
liste=glob.glob(direct_set+"/"+str(atom)+"_"+str(roman[int(ion-1)])+"_*set")
if len(liste) == 0:
	print("No line profile was computed for this ion.")
	sys.exit()
if len(liste) == 1:
	direct_lp=liste[0]
if len(liste) > 1:
	vec=[]
	for i in range(int(len(liste))):
		liste_split=liste[i].split('/')
		vec.append(liste_split[-1])
	direct_lp=direct_set+input("Please choose between "+', '.join(map(str, vec))+": ")

# Read the binary parameters file
# ===============================
fres = open(sys.argv[1], 'r')
type1_th, type2_th, Mdot1_th, Mdot2_th, M1_th, M2_th, R1_th, R2_th, a_th, ex_th, omega_th, phase_th, incl_th =([] for _ in range(13))
while 1:
	line=fres.readline()
	if not line: break
	vec=line.split(' ')
	if (vec[0] != '#'):
		if (len(vec) != 15):
			print("Your parameter file must contain 14 columns: Type 1  Type 2  Mdot1 [Msun/yr]  Mdot2 [Msun/yr]  M1 [Msun]  M2 [Msun]  R1 [Rsun]  R2 [Rsun]  Teff1 [K]  Teff2 [K]  a [Rsun]  ex  omega [degree]  phase  incl [degree]")
			sys.exit()
		else:
			type1_th.append(vec[0])
			type2_th.append(vec[1])
			Mdot1_th.append(float(vec[2]))
			Mdot2_th.append(float(vec[3]))
			M1_th.append(float(vec[4]))
			M2_th.append(float(vec[5]))
			R1_th.append(float(vec[6]))
			R2_th.append(float(vec[7]))
			a_th.append(float(vec[10]))
			ex_th.append(float(vec[11]))
			omega_th.append(float(vec[12]))
			phase_th.append(float(vec[13]))
			incl_th.append(float(vec[14]))
type1_th=np.array(type1_th)
type2_th=np.array(type2_th)
Mdot1_th=np.array(Mdot1_th)
Mdot2_th=np.array(Mdot2_th)
M1_th=np.array(M1_th)
M2_th=np.array(M2_th)
R1_th=np.array(R1_th)
R2_th=np.array(R2_th)
a_th=np.array(a_th)
ex_th=np.array(ex_th)
omega_th=np.array(omega_th)*math.pi/180.
phase_th=np.array(phase_th)
incl_th=np.array(incl_th)*math.pi/180.
fres.close()

print("")
print(("Directory: "+str(direct_lp)+"\n"))
print("Index of your binary (starts at 0)       |       Closest file")
print("-------------------------------------------------------------")
fres = open("results_find_LP.txt", 'w')

fres.write("Directory: "+str(direct_lp)+"\n")
fres.write("\n")
fres.write("Index of your binary (starts at 0)       |       Closest file\n")
fres.write("-------------------------------------------------------------\n")
for ipar in range(int(len(type1_th))):
	type1 = type1_th[ipar]
	type2 = type2_th[ipar]
	a = a_th[ipar] #Rsun
	Mdot1 = Mdot1_th[ipar] #Msun/yr
	Mdot2 = Mdot2_th[ipar] #Msun/yr
	M1 = M1_th[ipar] #Msun
	M2 = M2_th[ipar] #Msun
	R1 = R1_th[ipar] #Rsun
	R2 = R2_th[ipar] #Rsun
	omega=omega_th[ipar]
	phase=phase_th[ipar]
	ex=ex_th[ipar]
	incl=incl_th[ipar]

	# Conjunction phase. If ex==0, phase_conj=0
	M_conj=0.
	if (ex > 0):
		vinf1=2.6*math.sqrt(2.*grav*M1*Msun/(R1*Rsun)) #vinf=2.6*vesc
		vinf2=2.6*math.sqrt(2.*grav*M2*Msun/(R2*Rsun))
		beta=(Mdot2*vinf2)/(Mdot1*vinf1)
		E_conj=2.*math.atan(math.sqrt((1.-ex)/(1.+ex))*math.tan(0.5*(0.5*math.pi-omega)))
		M_conj=E_conj-ex*math.sin(E_conj)+math.pi*(beta < 1)

	phase=phase-M_conj
	if (phase<0):
		phase=2.*math.pi+phase

	M=phase
	E=M
	dE=(M-E+ex*math.sin(E))/(1.-ex*math.cos(E))
	while (abs(dE) >= 1.0E-8):
		E=E+dE
		dE=(M-E+ex*math.sin(E))/(1.-ex*math.cos(E))
	coi=(math.cos(E)-ex)/(1.-ex*math.cos(E))
	sii=math.sqrt(1.-ex**2)*math.sin(E)/(1.-ex*math.cos(E))
	phase_true=math.atan2(sii,coi)
	d=a*(1.-ex**2)/(1.+ex*math.cos(phase_true)) #Rsun

	liste=glob.glob(direct_lp+"/line_profile_par*.data")
	par_vec=[]
	for ilist in range(int(len(liste))):
		namei=liste[ilist]
		vec=namei.split('_')
		if (len(vec) >= 7):
			parname=vec[-3]
			par_vec.append(int(parname[3:]))
	par_vec=np.unique(par_vec)

	# Read set of parameters
	fparam = open("stellar_params_set.txt", 'r')
	type1_set, type2_set, d_set =([] for _ in range(3))
	line_count=0
	ipar_vec=0
	while 1:
		line=fparam.readline()
		if not line: break
		vec=line.split(' ')
		if (vec[0] != '#'):
			line_count=line_count+1
			if ipar_vec > len(par_vec)-1: break
			if (line_count != par_vec[ipar_vec]):
				continue
			ipar_vec=ipar_vec+1
			type1_set.append(vec[0])
			type2_set.append(vec[1])
			d_set.append(float(vec[10]))
	type1_set=np.array(type1_set)
	type2_set=np.array(type2_set)
	d_set=np.array(d_set)
	fparam.close()

	ipar2 = np.where((type2_set==type2) & (type1_set==type1))[0]
	if (len(ipar2) > 1):
		ipar3 = np.where(abs(d_set[ipar2]-d)==min(abs(d_set[ipar2]-d)))[0]
		ipar2 = ipar2[ipar3]
	if (len(ipar2) == 0):
		ipar2 = np.where((type1_set==type2) & (type2_set==type1))[0]
		if (len(ipar2) > 1):
			ipar3 = np.where(abs(d_set[ipar2]-d)==min(abs(d_set[ipar2]-d)))[0]
			ipar2 = ipar2[ipar3]
		if (len(ipar2) == 0):
			typeall=np.concatenate((type1_set,type2_set))
			typeall=np.unique(typeall)
			typeallstr=', '.join([x for x in typeall])
			print(("Allowed spectral type: "+typeallstr))
			print("For an other spectral type, please compute the radiative inhibition")
			sys.exit()
	if (isinstance(ipar2,np.ndarray)):
		ipar2=ipar2[0]

	liste=glob.glob(direct_lp+"/line_profile_par"+str(int(par_vec[ipar2]))+"*.data")
	phase_vec=[]
	incl_vec=[]
	for ilist in range(int(len(liste))):
		namei=liste[ilist]
		vec=namei.split('_')
		if (len(vec) >= 7):
			inclname=vec[-2]
			incl_vec.append(float(inclname[4:]))
			phasename=vec[-1]
			phasename=phasename[:-5]
			phase_vec.append(float(phasename[5:]))
	incl_vec=np.array(incl_vec)
	phase_vec=np.array(phase_vec)

	iphaseincl =  np.where((np.sqrt((phase_vec/100.-phase)**2+(incl_vec*math.pi/180.-incl)**2)==min(np.sqrt((phase_vec/100.-phase)**2+(incl_vec*math.pi/180.-incl)**2))))[0]
	if (isinstance(iphaseincl,np.ndarray)):
		iphaseincl=iphaseincl[0]

	print((str(ipar)+"   |   line_profile_par"+str(int(par_vec[ipar2]))+"_incl"+str(int(incl_vec[iphaseincl]))+"_phase"+str(int(phase_vec[iphaseincl]))+".data"))
	fres.write(str(ipar)+"   |   line_profile_par"+str(int(par_vec[ipar2]))+"_incl"+str(int(incl_vec[iphaseincl]))+"_phase"+str(int(phase_vec[iphaseincl]))+".data\n")

print("")
direct_lp=direct_lp.split('/')
print("The needed files can be downloaded with: wget http://lifeline.astro.uliege.be/Lifeline/"+str(direct_lp[-1])+"/[file_names]")
fres.write("")
fres.write("The needed files can be downloaded with: wget http://lifeline.astro.uliege.be/Lifeline/"+str(direct_lp[-1])+"/[file_names]")
fres.close()

