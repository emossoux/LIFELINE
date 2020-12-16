#!usr/bin/python
"""
******************************************************    
***   Program for the computation of line profile  ***
***           in colliding wind binaries           ***         
****************************************************** 
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

Run: python script_root+line_profile.py [code_parameters_file]

# Parameters 
# ==========
If the file code_parameters is given, the parameters are read in this file else, the program asks the parameters one by one.

List of paramters:
atom - Atom for which you want to compute the line (ex: Fe, Ca,...)
ion - Ion of this element (ex: 25, 14, ...)
energy - Energy of the line [keV]
directory - PATH where to create the directory to save the results in it. Ex: /home/mossoux/Documents/Post-doc
This directory must contain the cross section file "cross_section.tab", the "cooling_functions" directory, the param file, the apec_linelist.fits file, and optionally the ionisation fraction file "ionisation_fraction.tab" if it is already computed (crea_ion_file='no') or the filemap if not.
param - File containing the stellar parameters.
binstart - Index of the first binary to compute (begins at 0). Default value: 0.
binnum - Number of binaries to compute. Default value: the length of the stellar parameters file.
mode - In which mode do you want to compute? [1/2/3] 1 - Perform the overall computation, i.e., the velocity distribution, the shock characteristics and the line profile; 2 - Compute shock characteristics and the line profile using velocity distribution already computed; 3 - Compute only the line profile using the shock characteristics and the velocity distribution already computed.
crea_ion_file - compute the ionisation fraction file for radiative shocks? [yes/no]
convolve - Convolve the theoretical light curve with the instrumental response? [yes/no]
direct_rmf_arf - Only if convolve=yes. Complete path to the response matrix file (RMF) and ancillary response file (ARF)
RMF - Only if convolve=yes. Response matrix file
ARF - Only if convolve=yes. Ancillary response file
mu [val] - Mean molecular weight. Default value: 0.62 (totally ionized gas)
H [val] - Fractional mass abundance of H. Default value: 0.7381 (Sun)
betai [val] - Parameter of the beta law of the wind velocity for each star. Default value: 1.


# Versions
# ========
v1 - 16/05/2018 - Modify the FeK_radiative_paral.py script to choose beteen adiabatic/radiative and with/out Coriolis
v2 - 01/03/2019 - Add the radiative inhibition
v3 - 15/07/2019 - Loop on a grid of stellar parameters
v4 - 24/10/2019 - Let the user choose the stellar parameters
v5 - 05/03/2020 - Add the convolution option
"""      
from scipy.interpolate import interp1d            
import numpy as np
import os,glob
import sys
import math
from constantes import constante
import matplotlib.pyplot as plt  
import pyfits
import time


# Code parameters
# ===============
if (len(sys.argv) == 1):
        atom=raw_input("Atom for which you want to compute the line (ex: Fe, Ca,...): ") 
	atom=atom.capitalize()
	ion=int(float(raw_input("Ion of this atom (ex: 25, 14, ...): ")))
	energy=float(raw_input("Energy of the line [keV]: "))
	directory =raw_input("Complete path where to create the directory to save the results in it. \n This directory must contain the cross section file \"cross_section.tab\", the \"cooling_functions\" directory, the param file, the apec_linelist.fits file, and optionally the ionisation fraction file \"ionisation_fraction.tab\" if it is already computed (crea_ion_file='no') or the filemap if not: ")
	param =raw_input("File containing the stellar parameters: ")
	binstart = raw_input("Index of the first binary to compute (begins at 0 which is the default value): ")
	if binstart == "":
		binstart="0"
	binstart=float(binstart)
	binnum = raw_input("Number of binaries to compute (Default value: the length of the stellar parameters file): ")
	mode = int(float(raw_input("What do you want to compute? [1/2/3]\n1 - The overall computation, i.e., the velocity distribution, the shock characteristics and the associated line profile;\n2 - The shock characteristics for a different binary conbination and the associated line profile using velocity distribution already computed;\n3 - Only the line profile for an ion which is not already computed using the shock characteristics computed for the set of systems: ")))
	if mode == 1:
		inhibition='yes'
		comp_shock='yes'
	if mode == 2:
		inhibition='no'
		comp_shock='yes'
	if mode == 3:
		inhibition='no'
		comp_shock='no'
	crea_ion_file = raw_input("Compute the ionisation fraction file for radiative shocks? [yes/no] ")
	convolve = raw_input("Convolve the theoretical light curve with the instrumental response? [yes/no] ")
	if (convolve == 'yes'):
		direct_rmf_arf = raw_input("Complete path to the response matrix file (RMF) and ancillary response file (ARF): ")
		RMF=raw_input("Enter the RMF [h for help on this file]: ") 
		if (RMF == 'h'):
			print("The response matrix file (RMF) contains the propability that a photon of energy E is recorded in a channel I.")
			RMF=raw_input("Enter the RMF: ") 
		ARF=raw_input("Enter the ARF [h for help on this file]. If included in the RMF, leave blank: ") 
		if (ARF == 'h'):
			print("The ancillary response file (ARF) contains the effective area of the instrument.")
			ARF=raw_input("Enter the ARF: ") 
		distance=raw_input("Enter the distance to the source in kpc to obtain a convolved profile in erg/s. If left blank, a typical distance of 1.5kpc is assumed: ")
		if (distance == ''):
			distance=1.5
		else:
			distance=float(distance)
	print("Sun abundance. Default value: wilm") 
	print("  - wilm: Wilms, Allen & McCray (2000, ApJ 542, 914)")
	print("  - angr: Anders E. & Grevesse N. (1989, Geochimica et Cosmochimica Acta 53, 197)")
	print("  - aspl: Asplund M., Grevesse N., Sauval A.J. & Scott P. (2009, ARAA, 47, 481)")
	print("  - feld: Feldman U.(1992, Physica Scripta 46, 202)")
	print("  - aneb: Anders E. & Ebihara (1982, Geochimica et Cosmochimica Acta 46, 2363)")
	print("  - grsa: Grevesse, N. & Sauval, A.J. (1998, Space Science Reviews 85, 161)")
	print("  - lodd: Lodders, K (2003, ApJ 591, 1220)")
	sunabund = raw_input("Chosen sun abundance: ")
	if sunabund == "":
		sunabund="wilm"
	while (sunabund != "" and sunabund != "wilm" and sunabund != "angr" and sunabund != "aspl" and sunabund != "feld" and sunabund != "aneb" and sunabund !=  "grsa" and sunabund !=  "lodd"):
		sunabund = raw_input("Must be wilm, angr, aspl, feld, aneb, grsa or lodd")
	mu=0.62
	H_mass_frac=0.7381
	beta1=1.
	beta2=1.
	muuser = raw_input("Mean molecular weight (Default value: 0.62 for a totally ionized gas): ")
	if muuser != "":
		mu=float(muuser)
	Huser = raw_input("Fractional mass abundance of H (Default value: 0.7381 for the Sun): ")
	if Huser != "":
		H_mass_frac=float(Huser)
	betauser = raw_input("Parameter of the beta law of the wind velocity for star 1 (Default value: 1): ")
	if betauser != "":
		beta1=float(betauser)
	betauser = raw_input("Parameter of the beta law of the wind velocity for star 2 (Default value: 1): ")
	if betauser != "":
		beta2=float(betauser)

	fres = open(directory+"/results_line_profile.txt", 'w')
	fres.write("*** Computation of the line profile from line_profile.py ***\n")
	fres.write("\n")
	fres.write("Here are the code parameters you entered (may be copied in a new file to use it as a intput):\n")
	fres.write(atom+" #Atom for which you want to compute the line (ex: Fe, Ca,...)\n")
	fres.write(str(ion)+" #Ion of this atom (ex: 25, 14, ...)\n")
	fres.write(str(energy)+" #Energy of the line [keV]\n")
	fres.write(directory+" #Path where to create the directory to save the results in it. This directory must contain the cross section file 'cross_section.tab', the \"cooling_functions\" directory, the param file, the apec_linelist.fits file, and optionally the ionisation fraction file \"ionisation_fraction.tab\" if it is already computed (crea_ion_file='no') or the filemap if not.\n")
	fres.write(param+" #File containing the stellar parameters\n")
	fres.write(str(binstart)+" #Index of the first binary to compute (begins at 0). Default value: 0.\n")
	fres.write(str(binnum)+" #Number of binaries to compute. Default value: the length of the stellar parameters file.\n")
	fres.write(str(mode)+" #What do you want to compute? [1/2/3] 1 - The overall computation, i.e., the velocity distribution, the shock characteristics and the associated line profile; 2 - The shock characteristics for a different binary conbination and the associated line profile using velocity distribution already computed; 3 - Only the line profile for an ion which is not already computed using the shock characteristics computed for the set of systems.\n")
	fres.write(crea_ion_file+" #Compute the ionisation fraction file for radiative shocks? [yes/no]\n")
	fres.write(convolve+" #Convolve the theoretical light curve with the instrumental response? [yes/no]\n")
	if (convolve == 'yes'):
		fres.write(direct_rmf_arf+" #Path to the response matrix file (RMF) and ancillary response file (ARF)\n")
		fres.write(RMF+" #Response matrix file\n")
		fres.write(ARF+" #Ancillary response file\n")
		fres.write(str(distance)+" #Distance to the binary [kpc]\n")
	fres.write(sunabund+" #Sun abundance: wilm - Wilms, Allen & McCray (2000, ApJ 542, 914), angr - Anders E. & Grevesse N. (1989, Geochimica et Cosmochimica Acta 53, 197), aspl - Asplund M., Grevesse N., Sauval A.J. & Scott P. (2009, ARAA, 47, 481), feld - Feldman U.(1992, Physica Scripta 46, 202), aneb - Anders E. & Ebihara (1982, Geochimica et Cosmochimica Acta 46, 2363), grsa - Grevesse, N. & Sauval, A.J. (1998, Space Science Reviews 85, 161), lodd - Lodders, K (2003, ApJ 591, 1220)\n")
	fres.write(str(mu)+" #Mean molecular weight. 0.62 for totally ionized gas\n")
	fres.write(str(H_mass_frac)+" #Fractional mass abundance of H. 0.7381 for the Sun\n")
	fres.write(str(beta1)+" #Parameter of the beta law of the wind velocity for star 1. Default value: 1.\n")
	fres.write(str(beta2)+" #Parameter of the beta law of the wind velocity for star 2. Default value: 1.\n")
	fres.write("\n")


elif (len(sys.argv) == 2):
	# Read code parameter file                                           
	# ========================
	fparam = open(sys.argv[1], 'r')
	par_vec=[]
	while 1:
	    line=fparam.readline()
	    if not line: break
	    vec=line.split(' ')
	    par_vec.append(vec[0])
	fparam.close()

	if (len(par_vec) == 15):
        	atom, ion, energy, directory, param, binstart, binnum, mode, crea_ion_file, convolve, sunabund, mu, H_mass_frac, beta1, beta2=par_vec
		if (convolve == 'yes'):
			print("")
			print("Please enter the path and the names of the RMF and ARF as well as the distance to the binary to convolve the line profile.")
			sys.exit()
	elif (len(par_vec) == 19):
        	atom, ion, energy, directory, param, binstart, binnum, mode, crea_ion_file, convolve, direct_rmf_arf, RMF, ARF, distance, sunabund, mu, H_mass_frac, beta1, beta2=par_vec
	else: 
		print("")
		print("Your input file does not contain the right number of lines. Here are the lines it must contain:")
		print("Atom for which you want to compute the line (ex: Fe, Ca,...)")
		print("Ion of this atom (ex: 25, 14, ...)")
		print("Energy of the line [keV]")
		print("Complete path where to create the directory to save the results in it. \nThis directory must contain the cross section file 'cross_section.tab', the \"cooling_functions\" directory, the param file and optionally the ionisation fraction file 'ionisation_fraction.tab' if it is already computed (crea_ion_file='no').")
		print("File containing the stellar parameters")
		print("Index of the first binary to compute (begins at 0). Default value: 0.")
		print("Number of binary to compute. Default value: the length of the stellar parameters file.")
		print("What do you want to compute? [1/2/3]\n 1 - The overall computation, i.e., the velocity distribution, the shock characteristics and the associated line profile;\n 2 - The shock characteristics for a different binary conbination and the associated line profile using velocity distribution already computed;\n 3 - Only the line profile for an ion which is not already computed using the shock characteristics computed for the set of systems.")
		print("Compute the ionisation fraction file for radiative shocks? [yes/no]")
		print("Convolve the theoretical light curve with the instrumental response? [yes/no]")
		print("Complete path to the response matrix file (RMF) and ancillary response file (ARF). Only if convolve=yes.\nFor Athena, the files are given on the main directory of this program.\nRMF: athena_xifu_1190_onaxis_pitch249um_v20160401.rsp; ARF: athena_xifu_1190_onaxis_pitch249um_v20160401.arf.")
		print("Response matrix file. Only if convolve=yes.")
		print("Ancillary response file. Only if convolve=yes.")
		print("Distance to the binary [kpc]. Only if convolve=yes.")
		print("Sun abundance: wilm - Wilms, Allen & McCray (2000, ApJ 542, 914), angr - Anders E. & Grevesse N. (1989, Geochimica et Cosmochimica Acta 53, 197), aspl - Asplund M., Grevesse N., Sauval A.J. & Scott P. (2009, ARAA, 47, 481), feld - Feldman U.(1992, Physica Scripta 46, 202), aneb - Anders E. & Ebihara (1982, Geochimica et Cosmochimica Acta 46, 2363), grsa - Grevesse, N. & Sauval, A.J. (1998, Space Science Reviews 85, 161), lodd - Lodders, K (2003, ApJ 591, 1220). Default value: wilm")
		print("Mean molecular weight. 0.62 for totally ionized gas")
		print("Fractional mass abundance of H. 0.7381 for the Sun")
		print("Parameter of the beta law of the wind velocity for star 1. Default value: 1.")
		print("Parameter of the beta law of the wind velocity for star 2. Default value: 1.")
		sys.exit()
	atom=atom.capitalize()
	ion=int(float(ion))
	energy=float(energy)
	mode=int(mode)
	try:
		binstart=float(binstart)
	except ValueError:
		binstart=0
	try:
		mu=float(mu)
	except ValueError:
		mu=0.62
	try:
		H_mass_frac=float(H_mass_frac)
	except ValueError:
		H_mass_frac=0.7381
	try:
		beta1=float(beta1)
	except ValueError:
		beta1=1.
	try:
		beta2=float(beta2)
	except ValueError:
		beta2=1.
	
	if mode == 1:
		inhibition='yes'
		comp_shock='yes'
	if mode == 2:
		inhibition='no'
		comp_shock='yes'
	if mode == 3:
		inhibition='no'
		comp_shock='no'

	fres = open(directory+"/results_line_profile.txt", 'w')
	fres.write("*** Computation of the line profile from line_profile.py ***\n")
	fres.write("\n")
	fres.write("File used to define the code parameter: "+sys.argv[1]+"\n")
	fres.write("\n")
else: 
	print("")
        print("Run: python script_root+line_profile.py [code_parameter_file]")
        sys.exit()

Rsun=constante('Rsun') #cm
Msun=constante('Msun') #g
grav=constante('G')*1000. # gravity  in cm^3/g/s^2
light_vel=constante('c')*100. # proton mass in cm/s             
pi=math.pi

# Read parameter file                                           
# ===================
fparam = open(directory+"/"+param, 'r')        # System parameters
type1_th, type2_th, Mdot1_th, Mdot2_th, M1_th, M2_th, R1_th, R2_th, Teff1_th, Teff2_th, a_th, ex_th, omega_th, phase_th, incl_th =([] for _ in range(15))
while 1:
    line=fparam.readline()
    if not line: break
    vec=line.split(' ')
    if (vec[0] != '#'):
	if (len(vec) != 15):
		print("")
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
		Teff1_th.append(float(vec[8]))
		Teff2_th.append(float(vec[9]))
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
Teff1_th=np.array(Teff1_th)
Teff2_th=np.array(Teff2_th)
a_th=np.array(a_th)
ex_th=np.array(ex_th)
omega_th=np.array(omega_th)*pi/180.
phase_th=np.array(phase_th)*2.*pi
incl_th=np.array(incl_th)*pi/180.
fparam.close()
if (binnum==""):
	binnum=len(type1_th)-binstart
binnum=float(binnum)
if (binnum+binstart>len(type1_th)):
	print("The number of binaries you want to compute exceeds the number of lines of your parameter file.")
	sys.exit()

# Initialisation                                             
# ==============    
nbr_bin_profile=100
nbr_points_shock_2D=30  # from theta=0 to theta=theta_inf
nbr_points_shock_3D=50  # from 0 degree to 360 degree
nbr_points_width=20     # number of points inside the shock
cut_lim=0.85            # the velocity reach cut_lim% of the terminal velocity
fres.write("Number of bins in the line_profile: "+str(nbr_bin_profile)+"\n")
fres.write("Number of angular steps used to create the shock in 2D: "+str(nbr_points_shock_2D)+"\n")
fres.write("Number of angular steps used to expend the shock in 3D: "+str(nbr_points_shock_3D)+"\n")
fres.write("Number of radial steps used to discretize the interior of the shock: "+str(int(2.*nbr_points_width))+"\n")
fres.write("Velocity cutoff for the definition of the shock cap (percent of the terminal velocity): "+str(cut_lim*100.)+"\n")
fres.write("\n")

# List of the atoms until 30
# ==========================
atom_list=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn']
roman=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII', 'XIV', 'XV', 'XVI', 'XVII', 'XVIII', 'XIX', 'XX', 'XXI', 'XXII', 'XXIII', 'XXIV', 'XXV', 'XXVI', 'XXVII','XXVIII','XXIX','XXX']
if atom not in atom_list:
	print("")
	print("Please enter an atom between H and Zn.")
	sys.exit()
atomic_number=atom_list.index(atom)+1.
if (ion > atomic_number):
	print("")
	print("Ion > Atomic number. Please enter a valid ion for this element.")
	sys.exit()

# Emissivities
# ============
aa = pyfits.open(directory+'/apec_linelist.fits')
index = np.where((aa[1].data.field('Element') == atomic_number))[0]
good_ion=aa[1].data[index]
del aa
index = np.where((good_ion.field('Ion') == ion))[0]
good_ion=good_ion[index]
index = np.where((abs(good_ion.field('Lambda')-light_vel/(2.41838e9*energy)) == min(abs(good_ion.field('Lambda')-light_vel/(2.41838e9*energy)))))[0]
T=good_ion[index].field('Temperature')[0]
emis=good_ion[index].field('Emissivity')[0]
T=T[np.where(emis != 0)]
emis=emis[np.where(emis != 0)]

plt.figure(1)
plt.plot(T*8.61732e-8,emis/1.e-16)
plt.yscale('log')
plt.xlabel(r"$\mathrm{Temperature\ (\mathrm{keV})}$")
plt.ylabel(r"$\mathrm{Emissivity\ (10^{-16}\ photons\ cm^3\ s^{-1})}$")
plt.savefig(directory+"/Emissivity_"+atom+"_"+str(roman[int(ion-1)])+"_"+str(energy)+".pdf")
plt.close()
fres.write("Plot of the emissivity: Emissivity_"+atom+"_"+str(roman[int(ion-1)])+"_"+str(energy)+".pdf\n")
fres.write("\n")

emis=-np.log10(emis*1.986e-8/(10**-23*light_vel/(2.41838e9*energy)))
qinterp = interp1d(np.array(np.log10(T).astype('f8')), np.array(emis), kind='cubic')

if not os.path.exists(directory+"/histograms"):
	os.makedirs(directory+"/histograms")
histo_dir=directory+"/histograms"
if not os.path.exists(directory+"/plots"):
	os.makedirs(directory+"/plots")

compteur=0.
vec_par=np.arange(int(binnum))+binstart
vec_par=vec_par.astype(int)
for ipar in vec_par:
	if (ipar == 0):
		ipar_new=ipar
	else:
		if (ipar_new == ipar-1):
			ipar_new=ipar
		else:
			ipar_new=ipar_new+1
	if (os.path.exists(histo_dir+"/ray_tracing_par"+str(int(ipar))+"_bin1.data") and mode != 3):
		choix=raw_input("Files already exist for the binary system "+str(int(ipar))+". Do you want to replace them? [yes/no] ")
		if choix == "no":
			choix=raw_input("Increase the identifying number of the binary system? (the files will be numbered as nbr_systems_in_histograms+1) [yes/no] ")
			if choix == "no":
				print("")
				print("Please remove the files of the system parameter "+str(int(ipar)))
				sys.exit()
			else:
				ipar_new=ipar
				while os.path.exists(histo_dir+"/ray_tracing_par"+str(int(ipar_new))+"_bin1.data"):
					ipar_new=ipar_new+1
	sentence="*** Binary "+str(ipar+1)+" ***"
	print("")
	print("*"*len(sentence))
	print(sentence)
	print("*"*len(sentence))
	fres.write("*"*len(sentence)+"\n")
	fres.write(sentence+"\n")
	fres.write("*"*len(sentence)+"\n")
	fres.write("\n")
	type1 = type1_th[ipar]
	type2 = type2_th[ipar]
	Mdot1 = Mdot1_th[ipar] #Msun/yr
	Mdot2 = Mdot2_th[ipar] #Msun/yr
	M1 = M1_th[ipar] #Msun
	M2 = M2_th[ipar] #Msun
	R1 = R1_th[ipar] #Rsun
	R2 = R2_th[ipar] #Rsun
	Teff1 = Teff1_th[ipar]
	Teff2 = Teff2_th[ipar]
	a = a_th[ipar] #Rsun
	phase=phase_th[ipar]
	ex=ex_th[ipar]
	omega=omega_th[ipar]
	incl=incl_th[ipar]

	time1=time.time()

	mass_sum=(M1+M2)*Msun #g
	Per=2.*pi*(a*Rsun)**1.5/(math.sqrt(grav*mass_sum)*86400.)
	mass_ratio=M1/M2
	vinf1=2.6*math.sqrt(2.*grav*M1*Msun/(R1*Rsun)) #vinf=2.6*vesc
	vinf2=2.6*math.sqrt(2.*grav*M2*Msun/(R2*Rsun))
	beta=(Mdot2*vinf2)/(Mdot1*vinf1)

	M=phase
	E=M
	dE=(M-E+ex*math.sin(E))/(1.-ex*math.cos(E))
	while (abs(dE) >= 1.0E-8):
		E=E+dE
		dE=(M-E+ex*math.sin(E))/(1.-ex*math.cos(E))
	coi=(math.cos(E)-ex)/(1.-ex*math.cos(E))
	sii=math.sqrt(1.-ex**2)*math.sin(E)/(1.-ex*math.cos(E))
	phase_true=math.atan2(sii,coi)
	d=a*Rsun*(1.-ex**2)/(1.+ex*math.cos(phase_true)) #cm

	beta_old=beta
	if (beta < 1):
		print("")
		print("Stronger wind: star 1 of your file")
		fres.write("Stronger wind: star 1 of your file\n")
		fres.write("\n")
		Mdot1, Mdot2 = Mdot2, Mdot1
		vinf1, vinf2 = vinf2, vinf1
		beta=(Mdot2*vinf2)/(Mdot1*vinf1)
		R1, R2 = R2, R1
		Teff1, Teff2 = Teff2, Teff1
		mass_ratio = 1./mass_ratio
		beta1, beta2 = beta2, beta1
		type1, type2 = type2, type1
		xcm1, xcm2 = xcm2, xcm1


	# Compute the wind distribution
	if (inhibition == 'yes'):
		if not os.path.exists(directory+"/winds"):
			os.makedirs(directory+"/winds")
		print("")
		print("*** Computation of the wind distribution including radiative inhibition ***")
		fres.write("***  Wind distribution including with radiative inhibition computed ***\n")
		fres.write("\n")
		from radiative_inhibition import radiative_inhibition_comp
		radiative_inhibition_comp(type1, type2, Mdot1, Mdot2, Per, mass_ratio, R1, R2, Teff1, Teff2, d, directory, mu, H_mass_frac, beta1, beta2, ipar_new)
		wind_prim=directory+"/winds/wind_star1_param"+str(ipar_new)+".h5"
		wind_sec=directory+"/winds/wind_star2_param"+str(ipar_new)+".h5"
		fres.write("Files for the wind of each star: wind_star1_par"+str(ipar_new)+".h5 and wind_star2_par"+str(ipar_new)+".h5\n")
		fres.write("\n")
		time2=time.time()
	time2=time.time()
        compute_shock_radiative='no'
        compute_shock_adiabatic='no'

        if (compteur==0):
                # Read the cross sections
                # =======================
		fcs = directory+"/cross_section.tab"  # Cross section file
                lire=0
                with open(fcs) as f:
                        for line in f:
                                if (lire != 0):
                                        line2=line.split()
                                        ener=float(line2[0]) #keV
                                        cs_new=np.array(map(float,line2[1:]))
                                        if (ener == energy):
                                                break        
                                else:
                                        lire=1
                                        line2=line.split()
                                        T_cs=np.array(map(float, line2)) #keV
                cs_interp = interp1d(T_cs, cs_new, kind='cubic') # cross sections [cm^2] for different temperatures at the energy ENERGY
		compteur=1

	crashing = 'no'
	# Conjunction phase. If ex==0, phase_conj=0
	M_conj=0.
	if (ex > 0):
		E_conj=2.*math.atan(math.sqrt((1.-ex)/(1.+ex))*math.tan(0.5*(0.5*pi-omega)))
		M_conj=E_conj-ex*math.sin(E_conj)+pi*(beta_old < 1)
        # Compute the shock characterisics
        # ================================
	if (comp_shock == 'yes'):
		xstag=math.sqrt(1./beta)*d/(1+math.sqrt(1./beta)) #cm
		if (xstag < R1*Rsun):
			xstag=R1*Rsun
		vlocal=[vinf1*(1.-0.99*R1*Rsun/xstag)**beta1,vinf2*(1.-0.99*R2*Rsun/(d-xstag))**beta2] #cm/s
		print("")
		if ((min(vlocal*np.sign(vlocal))/1.0E8)**4*(xstag/1.0e12)/(max([Mdot1/(1.0e-7),Mdot2/(1.0e-7)]))>1):
			print("Adiabatic wind collision")
			fres.write("Adiabatic wind collision\n")
			fres.write("\n")
			compute_shock_adiabatic='yes'
		else:
			print("Radiative wind collision")
			fres.write("Radiative wind collision\n")
			fres.write("\n")
			#os.system("cp "+directory+"/cooling_functions/cooling_function_"+str(atom)+str(ion)+".tab "+directory+"/cooling_function.tab")
			compute_shock_radiative='yes'

		compute_shock_coriolis='no'
		vorb=math.sqrt(grav*mass_sum*(2./d-1./(a*Rsun)))
		if (vorb > 0.1*min(vlocal)):
			print("Coriolis effects included")
			fres.write("Coriolis effects included\n")
			fres.write("\n")
			compute_shock_coriolis = 'yes'

		if (inhibition == 'no'):
			# Read set of already computed stellar parameters
			fparam = open(directory+"/stellar_params_set.txt", 'r')
			type1_set, type2_set, d_set, mode_set =([] for _ in range(4))
			while 1:
			    line=fparam.readline()
			    if not line: break
			    vec=line.split(' ')
			    if (vec[0] != '#'):
				type1_set.append(vec[0])
				type2_set.append(vec[1])
				d_set.append(float(vec[10]))
				mode_set.append(int(vec[11]))

			type1_set=np.array(type1_set)
			type2_set=np.array(type2_set)
			mode_set=np.array(mode_set)
			d_set=np.array(d_set)
			fparam.close()

			ipar2 = np.where((mode_set == 1) & (type2_set==type2) & (type1_set==type1))[0]
			if (len(ipar2) > 1):
				ipar3 = np.where(abs(d_set[ipar2]-d/Rsun)==min(abs(d_set[ipar2]-d/Rsun)))[0]
				ipar2 = ipar2[ipar3]
			if (len(ipar2) == 0):
				ipar2 = np.where((mode_set == 1) & (type1_set==type2) & (type2_set==type1))[0]
				if (len(ipar2) > 1):
					ipar3 = np.where(abs(d_set[ipar2]-d/Rsun)==min(abs(d_set[ipar2]-d/Rsun)))[0]
					ipar2 = ipar2[ipar3]
				if (len(ipar2) == 0):
					print("")
					print("Allowed spectral type: O3I, O5I, O7I, O9I, O3III, O5III, O7III, O9III, O3V, O5V, O7V and O9")
					print("For an other spectral type, please compute the radiative inhibition")
					sys.exit()	
			if (isinstance(ipar2,np.ndarray)):
				ipar2=ipar2[0]
			wind_prim=directory+"/winds_set/wind_star1_param"+str(int(ipar2))+".h5"
			fres.write("Radiative inhibition not computed.\n")
			fres.write("Wind files used from the wind_set directory:\n")
			fres.write("  wind_star1_param"+str(int(ipar2))+".h5\n")
			wind_sec=directory+"/winds_set/wind_star2_param"+str(int(ipar2))+".h5"
			fres.write("  wind_star2_param"+str(int(ipar2))+".h5\n")
			fres.write("\n")

		# Compute the ionisation fraction
		if (compute_shock_radiative == 'yes' and crea_ion_file == "yes"):
			import pyatomdb
			os.environ["ATOMDB"] = 'https://hea-www.cfa.harvard.edu/AtomDB/'
			
			# write user_data (where filemap and APED are located)
			fuser = open(directory+"/userdata", 'w')
			fuser.write("filemap="+directory+"/filemap\n")
			fuser.write("atomdbroot=https://hea-www.cfa.harvard.edu/AtomDB/\n")
			fuser.close()
			# read user_data
			setting=pyatomdb.util.load_user_prefs(adbroot=directory)

			print("")
			print("*** Computation of the ionisation fraction ***")
			fres.write("*** Ionisation fraction computed ***\n")
			fres.write("\n")
			nz=25
			kT_ion=np.linspace(0.1,10.,200)
			test=np.arange(nz)+1
			fion = open(directory+"/ionisation_fraction.tab", 'w')
			fion.write("# Temperature [keV]  Z  Ionisation fraction\n")
			for ikt in kT_ion:
				ionpop=pyatomdb.apec.calc_full_ionbal(ikt, tau=1.0e14, init_pop=False, Te_init=False, Zlist=test.tolist(), teunit='keV', extrap=False, cie=True, settings=setting)
				for inz in range(nz):
				 	fion.write(str(ikt)+" "+str(inz+1)+" "+" ".join(map(str,ionpop[inz+1]))+"\n")
			fion.close()

		print("")
		print("*** Computation of the shock characteristics ***")
		fres.write("*** Shock computed ***\n")
		fres.write("\n")
		add=""
		add2=""
		if (compute_shock_coriolis == 'yes'):
			if (compute_shock_radiative == 'yes'):
				liste_param=Mdot1,Mdot2,vinf1,vinf2, mass_ratio,ex, omega,a, Per, R1, R2,beta1,beta2, Teff1, Teff2, nbr_points_shock_2D, nbr_points_shock_3D, nbr_points_width, nbr_bin_profile, directory, cut_lim, wind_prim, wind_sec, phase, incl, ipar_new, sunabund, mu, T, M_conj, atom, ion
				from rshock_coriolis import rshock
				rien=rshock(liste_param)
				print("skew angle (degree): "+str(rien[0]*180./pi))
				fres.write("skew angle (degree): "+str(rien[0]*180./pi)+"\n")
				if (rien[1] == 'yes'):
					add="_crashing"
					print("Crashing wind: the line profile computed from this region will not be correct since other physical processes must be taken into account.")
					fres.write("Crashing wind: the line profile computed from this region will not be correct since other physical processes must be taken into account.\n")
				fres.write("\n")
			if (compute_shock_adiabatic == 'yes'):
				liste_param=Mdot1,Mdot2,vinf1,vinf2, mass_ratio,ex, omega,a, Per, R1, R2,Teff1, Teff2, nbr_points_shock_2D, nbr_points_shock_3D, nbr_points_width, nbr_bin_profile, directory, cut_lim, wind_prim, wind_sec, phase, incl, ipar_new, mu, T, M_conj
				from ashock_coriolis import ashock
				skew=ashock(liste_param)
				print("skew angle (degree): "+str(skew*180./pi))
				fres.write("skew angle (degree): "+str(skew*180./pi)+"\n")
				fres.write("\n")
		else:
			if (compute_shock_radiative == 'yes'):
				liste_param=Mdot1,Mdot2,vinf1,vinf2, mass_ratio,ex, omega,a, Per, R1, R2,beta1,beta2, Teff1, Teff2, nbr_points_shock_2D, nbr_points_shock_3D, nbr_points_width, nbr_bin_profile, directory, cut_lim, wind_prim, wind_sec, phase, incl, ipar_new, sunabund, mu, M_conj, atom, ion
				from rshock import rshock
				crashing=rshock(liste_param)
				if (crashing == 'yes'):
					add="_crashing"
					print("Crashing wind: the line profile computed from this region will not be correct since other physical processes must be taken into account.")
					fres.write("Crashing wind: the line profile computed from this region will not be correct since other physical processes must be taken into account.\n")
				fres.write("\n")
			if (compute_shock_adiabatic == 'yes'):
				liste_param=Mdot1,Mdot2,vinf1,vinf2, mass_ratio,ex, omega,a, Per, R1, R2, Teff1, Teff2, nbr_points_shock_2D, nbr_points_shock_3D, nbr_points_width, nbr_bin_profile, directory, cut_lim, wind_prim, wind_sec, phase, incl, ipar_new, mu, M_conj
				from ashock import ashock
				rien=ashock(liste_param)
		time3=time.time()
		print("I have worked during "+str((time3-time2)/60.)+" min or "+str((time3-time2)/3600.)+" hours.")
		direct=directory+"/histograms"
		ipar2=ipar_new
	else:
		phase=phase-M_conj
		if (phase<0):
			phase=2.*pi+phase

		# Read set of already computed stellar parameters
		fparam = open(directory+"/stellar_params_set.txt", 'r')
		type1_set, type2_set, d_set =([] for _ in range(3))
		while 1:
		    line=fparam.readline()
		    if not line: break
		    vec=line.split(' ')
		    if (vec[0] != '#'):
			type1_set.append(vec[0])
			type2_set.append(vec[1])
			d_set.append(float(vec[10]))
		type1_set=np.array(type1_set)
		type2_set=np.array(type2_set)
		d_set=np.array(d_set)
		fparam.close()

		ipar2 = np.where((type2_set==type2) & (type1_set==type1))[0]
		if (len(ipar2) > 1):
			ipar3 = np.where(abs(d_set[ipar2]-d/Rsun)==min(abs(d_set[ipar2]-d/Rsun)))[0]
			ipar2 = ipar2[ipar3]
		if (len(ipar2) == 0):
			ipar2 = np.where((type1_set==type2) & (type2_set==type1))[0]
			if (len(ipar2) > 1):
				ipar3 = np.where(abs(d_set[ipar2]-d/Rsun)==min(abs(d_set[ipar2]-d/Rsun)))[0]
				ipar2 = ipar2[ipar3]
			if (len(ipar2) == 0):
				print("")
				print("Allowed spectral type: O3I, O5I, O7I, O9I, O3III, O5III, O7III, O9III, O3V, O5V, O7V and O9")
				print("For an other spectral type, please compute the radiative inhibition")
				sys.exit()
		if (isinstance(ipar2,np.ndarray)):
			ipar2=ipar2[0]

		liste=glob.glob(directory+"/histograms_set/emiss_part_shock_par"+str(ipar2)+"_*.data")
		phase_vec=[]
		incl_vec=[]
		for ilist in range(len(liste)):
			namei=liste[ilist]
			vec=namei.split('_')
			inclname=vec[-2]
			inclname=float(inclname[4:])
			incl_vec.append(inclname)
			phasename=vec[-1]
			phasename=phasename[:-5]
			phasename=float(phasename[5:])
			phase_vec.append(phasename)
		incl_vec=np.array(incl_vec)
		phase_vec=np.array(phase_vec)
		iphaseincl =  np.where((np.sqrt((phase_vec/100.-phase)**2+(incl_vec*pi/180.-incl)**2)==min(np.sqrt((phase_vec/100.-phase)**2+(incl_vec*pi/180.-incl)**2))))[0]
		if (isinstance(iphaseincl,np.ndarray)):
			iphaseincl=iphaseincl[0]
		add="_incl"+str(int(incl_vec[iphaseincl]))+"_phase"+str(int(phase_vec[iphaseincl]))
		fres.write("Characteristics of the shock not computed.\n")
		fres.write("Closest binary system (from stellar_params_set.txt): "+str(ipar2)+"\n")
		fres.write("Closest orbital parameters: phase="+str(phase_vec[iphaseincl]/100.)+" and incl="+str(incl_vec[iphaseincl])+" degree\n")
		direct=directory+"/histograms_set"
		os.chdir(direct)
		add2=add
		for files in glob.glob("bin_par"+str(ipar2)+add+"*.data"):
			if (files[-13:-5]== 'crashing'):
				print("")
				print("Crashing wind: the line profile computed from this region will not be correct since other physical processes must be taken into account.")
				fres.write("Crashing wind: the line profile computed from this region will not be correct since other physical processes must be taken into account.\n")
				crashing == 'yes'
				add=add+"_crashing"
	time3=time.time()


	# Compute the line profile
	# ========================
	os.chdir(directory)
	print("")
        print("*** Computation of the line profile of "+str(atom)+" "+str(roman[int(ion-1)])+" ***")
	from compute_profile import profile
	if not os.path.exists(directory+"/"+str(atom)+"_"+str(roman[int(ion-1)])):
		os.makedirs(directory+"/"+str(atom)+"_"+str(roman[int(ion-1)]))
	direct2=directory+"/"+str(atom)+"_"+str(roman[int(ion-1)])

	M=M-M_conj
	if (M<0):
		M=2.*pi+M

	p=profile(direct, direct2, ipar2, nbr_bin_profile, qinterp, T, cs_interp, T_cs, energy, M, add, add2)
	fres.write("*** Line profile computed ***\n")
	fres.write("\n")

	# Convolve the line profile
	# =========================
	if (convolve == 'yes'):
		from observed_line_profile_main import convolve
		print("*** Convolution of the line profile ***")
		pconv=convolve(direct2, "line_profile_par"+str(ipar2)+add+".data", direct_rmf_arf, RMF, ARF, float(distance))
		fres.write("*** Line profile convolved ***\n")
		fres.write("\n")

	time4=time.time()
	print("")
	print("I have worked during "+str(time4-time3)+" sec or "+str((time4-time3)/60.)+" min.")
print("******************************")
fres.close()

