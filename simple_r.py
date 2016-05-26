# PARAMTERS FOR THE RAMACHANDRAN NUMBER
# The ranges for phi and psi are [-180,180]. 
# Any other value will be garbled (so, remember to 
# convert your angles so that it fits this range.
import math
bound = 360.0 # This does not chang
rho_scaling = 10.0 # This is sigma in the manuscript (Mannige, Kundu, Whitelam, 2016)
multiplier = int(round(rho_scaling*(bound*(2.0**0.5)),0)) # For internal reference

# To get the raw Ramachandran number from phi and psi:
def raw_ramachandran_number_collapse(phi,psi):
	phi = float(phi)
	psi = float(psi)
	a = round(rho_scaling*(phi-psi + bound)/math.sqrt(2.0),0)
	b = round(rho_scaling*(phi+psi + bound)/math.sqrt(2.0),0)
	return round(a,0) + round(b,0)*multiplier
#
# Original phi and psi values are easily obtained from the *raw* Ramachandran number
def raw_ramachandran_number_expand(z):
	z = float(z)
	phi = (math.sqrt(2.0)*np.mod(z,multiplier)/rho_scaling     - bound + math.sqrt(2.0)*np.floor(z / multiplier)/rho_scaling - bound )/2.0
	psi = (math.sqrt(2.0)*np.floor(z / multiplier)/rho_scaling - bound - math.sqrt(2.0)*np.mod(z, multiplier)/rho_scaling + bound)/2.0
	return phi,psi

# First getting the lowest and highest possible unnormalized R numbers
raw_R_min = raw_ramachandran_number_collapse(-180,-180)
raw_R_max = raw_ramachandran_number_collapse(180,180)
# To get the normalized Ramachandran number from the raw Ramachandran number 
def normalized_ramachandran_number(phi,psi):
	raw_R = raw_ramachandran_number_collapse(phi,psi)
	R = float(raw_R - raw_R_min)/(raw_R_max-raw_R_min)
	return R

current_phi = -60
current_psi = -60
# To get the normalized Ramachandran number:
R = normalized_ramachandran_number(current_phi,current_psi)
print "R(phi="+str(current_phi)+",psi="+str(current_psi)+") = "+str(R)