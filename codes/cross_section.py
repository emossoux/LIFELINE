#!usr/bin/python
"""
***************************************************************         
***   Program for the computation of the cross section      ***         
*************************************************************** 
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

Run: python script_root+cross_section.py
Run from the directory where the ionization_fraction.tab file is.
Internet connection required.

# Parameters 
# ==========

# Versions
# ========
v1 - 9/8/2017
"""

# Subroutine: compute the cross section
# =====================================
def cs_sub(energy,kT,direct):
	import numpy as np
	from sun_abund import abund
	import pyatomdb
	# read user_data
	setting=pyatomdb.util.load_user_prefs(adbroot=direct)
	number_fraction=abund(ref="wilm")
	number_fraction=np.array(number_fraction)/sum(number_fraction)
	ionbalfile = direct+"/ionisation_fraction.tab"		# Ionisation balance file

	nz=27
	cs=0.
	for zz in range(nz):
		Z=zz+1
		Xz=pyatomdb.atomdb.get_ionfrac(ionbalfile, Z, kT/8.61732e-8, z1=-1) # Ionization fraction; assumes ionization equilibrium; 0=not ionized, -1=totally ionized
		for zzi in range(Z):
			Zi=zzi
			if (Xz[Zi] > 0):
				# load level data
				lvdata = pyatomdb.get_data(Z, Zi, 'LV')
				if (lvdata != False):
					# load XSTAR PI data if it exists
					if (lvdata[1].data['phot_type'][0] == 3):
						pidata = pyatomdb.get_data(Z, Zi, 'PI')
						cs=cs+number_fraction[Z]*Xz[Zi]*pyatomdb.atomdb.sigma_photoion(energy, Z, Zi, lvdata[1].data['phot_type'][0], lvdata[1].data['phot_par'][0], xstardata=pidata, xstarfinallev=1)
					else:
						if (np.mean(lvdata[1].data['phot_par'][0]) != 0.):
							cs=cs+number_fraction[Z]*Xz[Zi]*pyatomdb.atomdb.sigma_photoion(energy, Z, Zi, 2, lvdata[1].data['phot_par'][0], xstardata=False, xstarfinallev=1)
	return cs # Photoionisation cross section [cm^2]

# Main
# ====
import numpy as np
import os

# write user_data (where filemap and APED)
directory=os.getcwd()
fuser = open(directory+"/userdata", 'w')
fuser.write("filemap="+directory+"/filemap\n")
fuser.write("atomdbroot=http://sao-ftp.harvard.edu/AtomDB/\n")
fuser.close()

enerq=np.arange(3798)*(10.-0.05)/3797.+0.05
T_vec=np.arange(30)*(10.-0.5)/29.+0.5 #keV
fion = open("cross_section.tab", 'w')
fion.write(" "+" ".join(str(x) for x in T_vec)+"\n")
print("*** Energy ***")
for i in range(int(len(enerq))):
	print(("* "+str(enerq[i])+" keV *"))
	cssv=[]
	for j in range(int(len(T_vec))):
		cssv.append(cs_sub(enerq[i], T_vec[j],directory))
	fion.write(str(enerq[i])+" "+" ".join(str(x) for x in cssv)+"\n")
fion.close()

