#!usr/bin/python
"""
************************************ 
*      Solar number fraction       *
************************************
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

Run: 
  from sun_abund import abund
  number_fraction=abund(ref="name")

Abundances taken from: https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node117.html#abund

# Parametres du modele 
# ====================
ref - String. The reference for the abundances among:
	wilm - Wilms, Allen & McCray (2000, ApJ 542, 914)
	angr - Anders E. & Grevesse N. (1989, Geochimica et Cosmochimica Acta 53, 197)
	aspl - Asplund M., Grevesse N., Sauval A.J. & Scott P. (2009, ARAA, 47, 481)
	feld - Feldman U.(1992, Physica Scripta 46, 202)
	aneb - Anders E. & Ebihara (1982, Geochimica et Cosmochimica Acta 46, 2363)
	grsa - Grevesse, N. & Sauval, A.J. (1998, Space Science Reviews 85, 161)
	lodd - Lodders, K (2003, ApJ 591, 1220)

"""
def abund(ref=None):
	if ref is None: ref="wilm"

	if (ref=="wilm"):
		nh=1.00E+00
		nhe=9.77E-02
		nli=0
		nbe=0
		nb=0
		nc=2.40E-04
		nn=7.59E-05
		no=4.90E-04
		nf=0
		nne=8.71E-05
		nna=1.45E-06
		nmg=2.51E-05
		nal=2.14E-06
		nsi=1.86E-05
		np=2.63E-07
		ns=1.23E-05
		ncl=1.32E-07
		nar=2.57E-06
		nk=0
		nca=1.58E-06
		nsc=0
		nti=6.46E-08
		nv=0
		ncr=3.24E-07
		nmn=2.19E-07
		nfe=2.69E-05
		nco=8.32E-08
		nni=1.12E-06
		ncu=0
		nzn=0

	elif (ref=="angr"):
		nh=1.00E+00
		nhe=9.77E-02
		nli=1.45E-11
		nbe=1.41E-11
		nb=3.98E-10
		nc=3.63E-04
		nn=1.12E-04
		no=8.51E-04
		nf=3.63E-08
		nne=1.23E-04
		nna=2.14E-06
		nmg=3.80E-05
		nal=2.95E-06
		nsi=3.55E-05
		np=2.82E-07
		ns=1.62E-05
		ncl=3.16E-07
		nar=3.63E-06
		nk=1.32E-07
		nca=2.29E-06
		nsc=1.26E-09
		nti=9.77E-08
		nv=1.00E-08
		ncr=4.68E-07
		nmn=2.45E-07
		nfe=4.68E-05
		nco=8.32E-08
		nni=1.78E-06
		ncu=1.62E-08
		nzn=3.98E-08

	elif (ref=="aspl"):
		nh=1.00E+00
		nhe=8.51E-02
		nli=1.12E-11
		nbe=2.40E-11
		nb=5.01E-10
		nc=2.69E-04
		nn=6.76E-05
		no=4.90E-04
		nf=3.63E-08
		nne=8.51E-05
		nna=1.74E-06
		nmg=3.98E-05
		nal=2.82E-06
		nsi=3.24E-05
		np=2.57E-07
		ns=1.32E-05
		ncl=3.16E-07
		nar=2.51E-06
		nk=1.07E-07
		nca=2.19E-06
		nsc=1.41E-09
		nti=8.91E-08
		nv=8.51E-09
		ncr=4.37E-07
		nmn=2.69E-07
		nfe=3.16E-05
		nco=9.77E-08
		nni=1.66E-06
		ncu=1.55E-08
		nzn=3.63E-08

	elif (ref=="feld"):
		nh=1.00E+00
		nhe=9.77E-02
		nli=1.26E-11
		nbe=2.51E-11
		nb=3.55E-10
		nc=3.98E-04
		nn=1.00E-04
		no=8.51E-04
		nf=3.63E-08
		nne=1.29E-04
		nna=2.14E-06
		nmg=3.80E-05
		nal=2.95E-06
		nsi=3.55E-05
		np=2.82E-07
		ns=1.62E-05
		ncl=3.16E-07
		nar=4.47E-06
		nk=1.32E-07
		nca=2.29E-06
		nsc=1.48E-09
		nti=1.05E-07
		nv=1.00E-08
		ncr=4.68E-07
		nmn=2.45E-07
		nfe=3.24E-05
		nco=8.32E-08
		nni=1.78E-06
		ncu=1.62E-08
		nzn=3.98E-08

	elif (ref=="aneb"):
		nh=1.00E+00
		nhe=8.01E-02
		nli=2.19E-09
		nbe=2.87E-11
		nb=8.82E-10
		nc=4.45E-04
		nn=9.12E-05
		no=7.39E-04
		nf=3.10E-08
		nne=1.38E-04
		nna=2.10E-06
		nmg=3.95E-05
		nal=3.12E-06
		nsi=3.68E-05
		np=3.82E-07
		ns=1.89E-05
		ncl=1.93E-07
		nar=3.82E-06
		nk=1.39E-07
		nca=2.25E-06
		nsc=1.24E-09
		nti=8.82E-08
		nv=1.08E-08
		ncr=4.93E-07
		nmn=3.50E-07
		nfe=3.31E-05
		nco=8.27E-08
		nni=1.81E-06
		ncu=1.89E-08
		nzn=4.63E-08

	elif (ref=="grsa"):
		nh=1.00E+00
		nhe=8.51E-02
		nli=1.26E-11
		nbe=2.51E-11
		nb=3.55E-10
		nc=3.31E-04
		nn=8.32E-05
		no=6.76E-04
		nf=3.63E-08
		nne=1.20E-04
		nna=2.14E-06
		nmg=3.80E-05
		nal=2.95E-06
		nsi=3.55E-05
		np=2.82E-07
		ns=2.14E-05
		ncl=3.16E-07
		nar=2.51E-06
		nk=1.32E-07
		nca=2.29E-06
		nsc=1.48E-09
		nti=1.05E-07
		nv=1.00E-08
		ncr=4.68E-07
		nmn=2.45E-07
		nfe=3.16E-05
		nco=8.32E-08
		nni=1.78E-06
		ncu=1.62E-08
		nzn=3.98E-08

	elif (ref=="lodd"):
		nh=1.00E+00
		nhe=7.92E-02
		nli=1.90E-09
		nbe=2.57E-11
		nb=6.03E-10
		nc=2.45E-04
		nn=6.76E-05
		no=4.90E-04
		nf=2.88E-08
		nne=7.41E-05
		nna=1.99E-06
		nmg=3.55E-05
		nal=2.88E-06
		nsi=3.47E-05
		np=2.88E-07
		ns=1.55E-05
		ncl=1.82E-07
		nar=3.55E-06
		nk=1.29E-07
		nca=2.19E-06
		nsc=1.17E-09
		nti=8.32E-08
		nv=1.00E-08
		ncr=4.47E-07
		nmn=3.16E-07
		nfe=2.95E-05
		nco=8.13E-08
		nni=1.66E-06
		ncu=1.82E-08
		nzn=4.27E-08

	else:
		print("This reference is undefined")
		return []

	return [nh, nhe, nli, nbe, nb, nc, nn, no, nf, nne, nna, nmg, nal, nsi, np, ns, ncl, nar, nk, nca, nsc, nti, nv, ncr, nmn, nfe, nco, nni, ncu, nzn]
