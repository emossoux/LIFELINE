#!usr/bin/python
"""
******************************************************    
***   Program for the convolution of line profile  ***
***         with the instrumental response         ***         
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

Run: python /PATH/observed_line_profile.py

# Parameters 
# ==========
direct_LP - Directory where the line profile file is located
file_LP - file containing the theoretical line profile
direct_rmf_arf - Directory where the response matrix file (RMF) and ancillary response file (ARF) are located
RMF - Response matrix file
ARF - Ancillary response file
distance - Distance to the observed binary [kpc]. If not given, a typical distance of 1.5kpc is assumed

For Athena, you can use the files computed for a XIFU mirror module radius Rmax=1190mm, a 2.3mm rib spacing and a on-axis case given in the main directory of the code.
    RMF: XIFU_CC_BASELINECONF_THICKFILTER_2018_10_10.rmf.
    ARF: XIFU_CC_BASELINECONF_THICKFILTER_2018_10_10.arf.
# Versions
# ========
v1 - 03/03/2020
"""      
from scipy.interpolate import interp1d            
import numpy as np
import os
import sys
import math
from constantes import constante
import matplotlib.pyplot as plt  
import pandas as pd
import pyfits

print("")
print("Note: The ARF and RMF files of Athena are given in the main directory of the line profile program.")
print("Taken for a XIFU mirror module radius Rmax=1190mm, a 2.3mm rib spacing and a on-axis case.")
print("RMF: XIFU_CC_BASELINECONF_THICKFILTER_2018_10_10.rmf.")
print("ARF: XIFU_CC_BASELINECONF_THICKFILTER_2018_10_10.arf.")
print("")

# Code parameters
# ===============
direct_LP=raw_input("Enter the directory where the line profile file is located: ") 
file_LP=raw_input("File containing the theoretical line profile: ") 
direct_rmf_arf=raw_input("Enter the directory where the response matrix file (RMF) and ancillary response file (ARF) are located: ")
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
expo_time=raw_input("Enter the observing time in seconds to obtain the number of counts observed in the line profile. If left blank, a typical observation time of 10ks is assumed: ")
if (expo_time == ''):
	expo_time=10000.
else:
	expo_time=float(expo_time)


# Read line profile file                                           
# ======================
fLP = open(direct_LP+"/"+file_LP, 'r')
vtang=[]
emiss_th=[]
phrases=[]
while 1:
    line=fLP.readline()
    if not line: break
    vec=line.split(' ')
    if (vec[1] == 'Energy'):
	energy=float(vec[4])
    if (vec[0] != '#'):
	vtang.append(float(vec[0]))
	emiss_th.append(float(vec[1]))
    else:
	phrases.append(line)
fLP.close()
emiss_th=np.array(emiss_th)[np.argsort(-np.array(vtang))] #10^27 erg/s
vtang=-np.sort(-np.array(vtang))

wavelength=12.3984193/energy
wave=wavelength*vtang*1.0e3/constante('c')+wavelength
energy_th=12.3984193/wave
emiss_th=1.e27*emiss_th/(4.*math.pi*(distance*1.e5*constante('pc_m'))**2) #10^27 erg/s/cm^2
bins=(energy_th[2:]-energy_th[:-2])/2.
bin_length=np.concatenate(([energy_th[1]-energy_th[0]],bins,[energy_th[-1]-energy_th[-2]])) #keV
bin_length=np.median(bin_length)

# Discretize the RMF                                             
# ================== 
rmf = pyfits.open(direct_rmf_arf+'/'+RMF)
rmf_header = rmf[1].header
ext_1=1
ext_2=2
if (not 'ENERG_LO' in rmf_header):
	ext_1=2
	ext_2=1
energy_lo=rmf[ext_1].data.field('ENERG_LO')
energy_hi=rmf[ext_1].data.field('ENERG_HI')
matrix=rmf[ext_1].data.field('MATRIX') # variable length array
f_chan=rmf[ext_1].data.field('F_CHAN') # chanel number where matrix start
if (len(f_chan.shape) > 1):
	f_chan=f_chan[:][0]
n_chan=rmf[ext_1].data.field('N_CHAN') # number of chanels in the matrix
if (len(n_chan.shape) > 1):
	n_chan=n_chan[:][0]
channel=rmf[ext_2].data.field('CHANNEL')
e_min=rmf[ext_2].data.field('E_MIN')
e_max=rmf[ext_2].data.field('E_MAX')
energy_lo_hi=(e_max+e_min)/2.
matrix_total=[]
if (min(energy_lo)>max(energy_th+bin_length/2.) or max(energy_hi)<min(energy_th-bin_length/2.)):
	print("At least part of the line profile is not covered by the RMF.")
	sys.exit()
	
for ibin in range(len(energy_th)):
	energy_bin=energy_th[ibin]
	ind1=np.where((energy_lo<energy_bin-bin_length/2.))[0] #begin at ind1[-1]
	ind2=np.where((energy_hi>energy_bin+bin_length/2.))[0] #end at ind2[0]

	matrix_new=np.zeros((len(channel),int(ind2[0]-ind1[-1]+1)))

	for iind in range(int(ind2[0]-ind1[-1]+1)):
		f_chan_here=f_chan[int(iind+ind1[-1])]
		n_chan_here=n_chan[int(iind+ind1[-1])]
		cum_sum=np.cumsum(n_chan_here)
		if isinstance(n_chan_here,list):
			ind_cumsum=0
			if (len(n_chan_here)>1):
				ind=np.argmax(matrix[int(iind+ind1[-1])][:])
				ind_cumsum=np.where((cum_sum<ind))[0]
				ind_cumsum=ind_cumsum[0]-1
			f_chan_here=f_chan_here[ind_cumsum]
			n_chan_here=n_chan_here[ind_cumsum]
			matrix_here=matrix[int(iind+ind1[-1])][ind_cumsum:ind_cumsum+n_chan_here+1]
		else:
			matrix_here=matrix[int(iind+ind1[-1])][:]
		
		matrix_new[f_chan_here:f_chan_here+n_chan_here,iind]=[x for x in matrix_here]
	mean_matrix=np.mean(matrix_new,axis=1)

	finterp = interp1d(energy_lo_hi.astype('f8'), mean_matrix.astype('f8'), kind='cubic')
	matrix_interp=finterp(energy_th)
	matrix_total.append(matrix_interp)

# Read the ARF                                             
# ============ 
arf_bin=energy_th*0.+1.
if (ARF != ''):
	arf = pyfits.open(direct_rmf_arf+'/'+ARF)
	energy_lo=arf[1].data.field('ENERG_LO')
	energy_hi=arf[1].data.field('ENERG_HI')
	specresp=arf[1].data.field('SPECRESP') #cm^2
	energy_lo_hi=(energy_lo+energy_hi)/2.
	finterp = interp1d(energy_lo_hi.astype('f8'), specresp.astype('f8'), kind='cubic')
	arf_bin=finterp(energy_th)

# Convolve                                             
# ========
LP_conv=[]
for ibin in range(len(energy_th)) :
	LP_conv.append(abs(np.sum(emiss_th*matrix_total[:][ibin]*arf_bin))) #erg/s

# Number of photons received
# ==========================
LP_conv=expo_time*np.array(LP_conv)*6.242e8/np.array(energy_th)
nbr_photon=np.nansum(LP_conv)
print("The number of photons observed during "+str(expo_time/1000.)+"ks  is: "+str(int(nbr_photon)))
LP_conv=LP_conv/expo_time

# Save                                             
# ====
fprofile = open(direct_LP+"/"+file_LP[:-5]+"_convolved.data", 'w')
for iphrase in range(len(phrases)):
	if (iphrase == len(phrases)-1):
		fprofile.write("# The number of photons observed during "+str(expo_time/1000.)+"ks  is "+str(int(nbr_photon))+"\n")
		fprofile.write("# tangential velocity (km/s) | spectrum (photon/s)\n")
		break
	fprofile.write(phrases[iphrase])
for i in range(len(vtang)):
	fprofile.write(str(vtang[-1])+" "+str(LP_conv[i])+"\n")
fprofile.close()

LP_conv=LP_conv/np.nansum(LP_conv)
y_limit=max(LP_conv)*1.1
plt.plot(energy_th,LP_conv)
plt.plot([energy,energy],[0.,1.],c='r')
plt.ylim(0.0,y_limit)
plt.xlabel("Energy (keV)")
plt.ylabel("Normalized convolved emission")
plt.savefig(direct_LP+"/"+file_LP[:-5]+"_convolved.pdf")
plt.close()

print("")
print("The output file is "+file_LP[:-5]+"_convolved.pdf")
