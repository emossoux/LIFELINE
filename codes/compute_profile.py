#!usr/bin/python
"""
**********************************************    
***   Program for the computation of the   ***
***              line profile              ***         
**********************************************
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

Run: from compute_profile import profile
     p=profile(direct, direct2, ipar, nbr_bin_profile, q1, T, cs_interp, T_cs, energy, M, add, add2, ipar_set)

Read the files created by [r/a]shock[_coriolis].py in direct

# Parameters 
# ==========
direct - directory where the shock characteritics are saved
direct2 - directory where to save files
ipar - index on the set of stellar parameter
nbr_bin_profile - number of bins in the line profile
q1 - interpolation function of the emissivity
T - tabulated temperatures used to define q1 
cs_interp - interpolation function of the cross-section
T_cs - tabulated temperatures used to define cs_interp
energy - rest energy of the line
M - orbital phase
add[2] - additional name when reading the set of binary systems
ipar_set - index of the binary system for which we read the ray-tracing files

# Output
# ======
p - numpy array containing the emission in each bin

# Versions
# ========
v1 - 21/10/19
"""      

# Main
# ====
def profile(direct, direct2, ipar, nbr_bin_profile, q1, T, cs_interp, T_cs, energy, M, add, add2, ipar_set):
	import numpy as np
	import math
	from constantes import constante
	import matplotlib.pyplot as plt

	fbin = open(direct+"/bin_par"+str(ipar_set)+add+".data", 'r')
	while 1:
		line=fbin.readline().split()
		bin1=line
    		if not line: break
		try:
			float(line[0])
			break
		except ValueError:
			pass
	fhisto = open(direct+"/emiss_part_shock_par"+str(ipar_set)+add2+".data", 'r')
	while 1:
		line=fhisto.readline().split() 
		histo1=line
    		if not line: break
		try:
			float(line[0])
			break
		except ValueError:
			pass
	ftemp = open(direct+"/temp_emiss_par"+str(ipar_set)+add2+".data", 'r')
	while 1:
		line=ftemp.readline().split()
		temp1=line
    		if not line: break
		try:
			float(line[0])
			break
		except ValueError:
			pass

	vtang=[]
	nom_debut="/line_profile_par"+str(ipar)+add
	fprofile = open(direct2+nom_debut+".data", 'w')	# Where to write the line profile
	fprofile.write("# Energy (keV) = "+str(energy)+"\n")
	fprofile.write("# Phases (mean anomaly, radian, M=0 -> conjunction): "+str(M)+"\n")		
	fprofile.write("# tangential velocity (km/s) | emissivity at each phase (10^27 erg/s)\n")

	emiss=[]
	for ibin in range(int(nbr_bin_profile)):
		emiss.append(0)
		if (ibin == 0):
			vtang.append(float(bin1[0])) #(km/s)
			facteur=histo1
			temper=temp1
		else:
			vtang.append(float(fbin.readline())) #(km/s)
			facteur=fhisto.readline().split()
			temper=ftemp.readline().split()
		facteur = [float(x) for x in facteur]
		temper = [float(x) for x in temper]
		fRT = open(direct+"/ray_tracing_par"+str(ipar_set)+add2+"_bin"+str(ibin)+".data", 'r')
		while 1:
			ligne=fRT.readline().split()
	    		if not ligne: break
			try:
				float(ligne[0])
				break
			except ValueError:
				pass
		for itemp in range(len(temper)):
              		if (facteur[itemp]>0):
				#emissivity
				if (temper[itemp]/8.61732e-8 < T[0] or temper[itemp]/8.61732e-8 > T[-1]):
					continue
				else:
					q = q1(math.log10(temper[itemp]/8.61732e-8)) # erg cm^3/s

				if itemp > 0:
					ligne=fRT.readline().split()
				temperRT=ligne[0::2]
				factRT=ligne[1::2]
				temperRT = [float(x) for x in temperRT]
				factRT = [float(x) for x in factRT]

				tau = 0.
				for itemp2 in range(len(temperRT)):
					# Photoionisation cross section
					if (temperRT[itemp2]<min(T_cs)):
						css = cs_interp(min(T_cs))
					elif (temperRT[itemp2]>max(T_cs)):
						css = cs_interp(max(T_cs))
					else:
						css = cs_interp(temperRT[itemp2]) #temperRT [keV]
					
					tau = tau + css*factRT[itemp2]
				if (np.isfinite(facteur[itemp]*np.exp(-tau)*10.**(-q))):
					emiss[ibin] = emiss[ibin] + facteur[itemp]*np.exp(-tau)*10.**(-q)
		fRT.close()

		fprofile.write(str(vtang[-1])+" "+str(emiss[ibin])+"\n")
	ftemp.close()
	fhisto.close()
	fbin.close()
	fprofile.close()

	plt.figure(1)
	wavelength=12.3984193/energy
	wave=wavelength*np.array(vtang)*1000./constante('c')+wavelength # vtang and c in m/s

	if (max(emiss) <= 0):
		print("No X-ray emission")
		y_to_plot=wave*0.
		y_limit=1.
	else:
		y_to_plot=emiss/np.nansum(emiss)
		y_limit=max(emiss/np.nansum(emiss))*1.1

	plt.plot(12.3984193/wave,y_to_plot)
	plt.plot([energy,energy],[0.,1.],c='r')
	plt.ylim(0.0,y_limit)
	plt.xlabel("Energy (keV)")
	plt.ylabel("Normalized emission")
	plt.savefig(direct2+nom_debut+".pdf")
	plt.close()


	return emiss

