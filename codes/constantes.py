#!usr/bin/python
"""
*********************************************** 
*      Definition of physical constants       *
***********************************************
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
  from constantes import constante
  C=constante(CName)

# Parameters 
# ==========
CName - The name of the constante you want among:
	Rsun [cm]
	Msun [g]
	Lsun [W]
	sigma_k (Stefan-Boltzmann in W/m^2/K^4)
	G (gravity  in m^3/kg/s^2)
	AU_Rsun (1 au in Rsun)
	AU_m (1 au in m)
	c (light velociy in m/s)
	pc_m (1 pc in m)
	K_kev (1 K in keV)
	sigma_t (Thomson cross section in m^2)
	me (electron mass in kg)
	mp (proton mass in kg)
	mu0 (Permeability of free space in N/A/m)
	e (electric charge in C)
	alpha (fine structure constant)
	kb (Boltzmann constant in erg/K)
	h (planck constant in eV s)

"""
def constante(CName):
	if (CName == 'Rsun'):
		return 6.957E10
	elif (CName == 'Msun'):
		return 1.991E33
	elif (CName == 'Lsun'):
		return 3.828E26
	elif (CName == 'sigma_k'):
		return 5.670373E-8
	elif (CName == 'G'):
		return 6.67408E-11
	elif (CName == 'AU_Rsun'):
		return 215.032155
	elif (CName == 'AU_m'):
		return 149597870700.
	elif (CName == 'c'):
		return 299792458.
	elif (CName == 'pc_m'):
		return 3.0857E16
	elif (CName == 'K_kev'):
		return 1.16E7
	elif (CName == 'sigma_t'):
		return 6.6524E-29
	elif (CName == 'me'):
		return 9.1094E-31
	elif (CName == 'mp'):
		return 1.6726E-27
	elif (CName == 'mu0'):
		return 1.2566E-6
	elif (CName == 'e'):
		return 1.6021766208e-19 
	elif (CName == 'alpha'):
		return 0.0072973525664 
	elif (CName == 'kb'):
		return 1.38064852e-16
	elif (CName == 'h'):
		return 4.135667662e-15
	else:
		print("This constant is not defined. Please add it in constantes.py or choose in Rsun, Msun, Lsun, sigma_k, G, AU_Rsun, AU_m, c, pc_m, K_kev, sigma_t,me, mp, mu0, e, alpha, kb or h")
		return None

