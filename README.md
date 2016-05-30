--------------------------------------------------
<pre>        _       _     __  __    _    ____  
  _ __ | | ___ | |_  |  \/  |  / \  |  _ \ 
 | '_ \| |/ _ \| __| | |\/| | / _ \ | |_) |
 | |_) | | (_) | |_  | |  | |/ ___ \|  __/ 
 | .__/|_|\___/ \__| |_|  |_/_/   \_\_|      v 0.0.0.0.0.0.0...
 |_|                 (Multi-angle Picture)
</pre>
This tool provides easily readable "pictures" of protein conformations, 
ensembles, and trajectories saved as either a combined protein databank 
(PDB) structure file, or a directory of such files, and produces graphs.

This tool provides easily readable "pictures" of protein conformations, 
ensembles, and trajectories saved as either a combined protein databank 
(PDB) structure file, or a directory of such files, and produces graphs.
-----<br>
#Usage<br>
-----<br>
```python plotmap.py -pdb ProteinDatabankStructureFilename.pdb
python plotmap.py -pdb /directory/containing/pdbs/```
------<br>
#Output (the x-axis always represents the models/structures listed in the PDB)<br>
------<br>
filename.rcode.eps      (y-axis: residue #; color: R number based on "-signed" and <rcode_cmap>)<br>
filename.rcode.his.eps  (y-axis: Ramachandran number (R); color: frequency of R in model)<br>
filename.rcode.rmsf.eps (y-axis: residue #; color: RMSF in R from the previous model)<br>
---------------<br>
#Additional tags<br>
---------------<br>
-h       -     Prints this message<br>
-ss      -     Color the ramachandran number codes (R-codes) by <br>
               secondary structure (default: color by chirality and sign)<br>
-signed  -     Use the signed version of the ramachandran number<br>
-rmsd    -     Also producee "filename.rcode.rmsd.eps"<br>
               (y-axis: residue #; color: RMSD in R from first model)<br>
---------------<br>
Each graph is also accompanied by "_colorbar.eps", which are keys.<br>
---------------<br>
The Ramachandran number concept is discussed in the manuscript:<br>
Mannige, Kundu, Whitelam (2016) "The Ramachandran number: an order parameter for protein geometry" <br>
Preprint at: http://arxiv.org/abs/1511.03011<br>
--------------------------------------------------

If you just want the code to calculate a Ramachandran number, then don't worry about cloning this Github account: just use the python code at the bottom:

```python
# PARAMTERS FOR THE RAMACHANDRAN NUMBER
# The ranges for phi and psi are [-180,180]. 
# Any other value will be garbled (so, remember to 
# convert your angles so that it fits this range.
import math
bound = 360.0 # This does not chang
sigma = 10.0 # This is sigma in the manuscript (Mannige, Kundu, Whitelam, 2016)
multiplier = int(round(sigma*(bound*(2.0**0.5)),0)) # For internal reference
multiplier_by_two = round(sigma*(bound*(2.0**0.5)),0)/2.0 # For internal reference

# To get the raw Ramachandran number from phi and psi:
def raw_ramachandran_number_collapse(phi,psi):
	phi = float(phi)
	psi = float(psi)
	a = round(sigma*(phi-psi + bound)/math.sqrt(2.0),0)
	b = round(sigma*(phi+psi + bound)/math.sqrt(2.0),0)
	return a + b*multiplier
#
# Original phi and psi values are easily obtained from the *raw* Ramachandran number
def raw_ramachandran_number_expand(z):
	z = float(z)
	phi = (math.sqrt(2.0)*np.mod(z,multiplier)/sigma     - bound + math.sqrt(2.0)*np.floor(z / multiplier)/sigma - bound )/2.0
	psi = (math.sqrt(2.0)*np.floor(z / multiplier)/sigma - bound - math.sqrt(2.0)*np.mod(z, multiplier)/sigma + bound)/2.0
	return phi,psi

# First getting the lowest and highest possible unnormalized R numbers
raw_R_min = raw_ramachandran_number_collapse(-180,-180)
raw_R_max = raw_ramachandran_number_collapse(180,180)

# If signed == 0, then we have the normal R number ranging from 0 to 1.
# If signed != 0, then the R number ranges from -1 to 1. Here, those R
#                 numbers associated with regions to the RIGHT of the 
#                 positive-sloped diagonal are multipled by -1.
# Unsigned version of this function returns values identical to:
# (raw_ramachandran_number_collapse(phi,psi)-raw_R_min)/(raw_R_max-raw_R_min)
def normalized_ramachandran_number(phi,psi, signed=0):
	phi = float(phi)
	psi = float(psi)
	a = round(sigma*(phi-psi + bound)/math.sqrt(2.0),0)
	b = round(sigma*(phi+psi + bound)/math.sqrt(2.0),0)
	raw_r = a + b*multiplier # identical to raw_ramachandran_number_collapse(psi,phi)
	final_r = float(raw_r - raw_R_min)/float(raw_R_max - raw_R_min)
	if signed:
		if a >= multiplier_by_two:
			final_r = final_r * -1.0
	return final_r

# To see how the Ramachandran number function works, set "if 0:" to "if 1:"
# 
current_phi = -55
current_psi = -60
# unsigned R
R = normalized_ramachandran_number(current_phi,current_psi,signed=0)
print "       R(phi="+str(current_phi)+",psi="+str(current_psi)+") = "+str(R)
# signed R
R = normalized_ramachandran_number(current_phi,current_psi,signed=1)
print "signed R(phi="+str(current_phi)+",psi="+str(current_psi)+") = "+str(R)
print
# 
current_phi = -65
current_psi = -60
# unsigned R
R = normalized_ramachandran_number(current_phi,current_psi,signed=0)
print "       R(phi="+str(current_phi)+",psi="+str(current_psi)+") = "+str(R)
# signed R
R = normalized_ramachandran_number(current_phi,current_psi,signed=1)
print "signed R(phi="+str(current_phi)+",psi="+str(current_psi)+") = "+str(R)
exit()
```


