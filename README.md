# plotmap
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

Usage:
python plotmap.py -pdb ProteinDatabankStructureFilename.pdb
python plotmap.py -pdb /directory/containing/pdbs/

Output (the x-axis always represents the models/structures listed in the PDB):<br>
filename.rcode.his.eps  (y-axis: Ramachandran number (R); color: frequency of R in model)<br>
filename.rcode.ss.eps   (y-axis: residue #; color: by secondary structure HELIX: red, SHEET: blue, PPII: cyan)<br>
filename.rcode.raw.eps  (y-axis: residue #; color: by chirality L: Blue, D: Red: Extended: White)<br>
filename.rcode.rmsd.eps (y-axis: residue #; color: RMSD in R from first model)<br>
filename.rcode.rmsf.eps (y-axis: residue #; color: RMSF in R from the previous model)

Additionally, each graph is accompanied by "_colorbar.eps", which are keys.

The Ramachandran number concept is based on the manuscript:<br>
Mannige, Kundu, Whitelam (2016) "The Ramachandran number: an order parameter for protein geometry" <br>
Preprint at: http://arxiv.org/abs/1511.03011

--------------------------------------------------

If you just want the code to calculate a Ramachandran number, then don't worry about cloning this Github account: just use the python code at the bottom:

```
# PARAMTERS FOR THE RAMACHANDRAN NUMBER
bound = 360.0 # This does not change
rho_scaling = 10.0 # This is sigma in the manuscript (Mannige, Kundu, Whitelam, 2016)
multiplier = int(round(rho_scaling*(bound*(2.0**0.5)),0)) # For internal reference

def ramachandran_number_collapse(x,y):
	a = round(rho_scaling*(x-y + bound)/numpy.sqrt(2.0),0)
	b = round(rho_scaling*(x+y + bound)/numpy.sqrt(2.0),0)
	return round(a,0) + round(b,0)*multiplier
#
def ramachandran_number_collapse_zigzag(x,y):
	a = round(rho_scaling*(x-y + bound)/numpy.sqrt(2.0),0)
	b = round(rho_scaling*(x+y + bound)/numpy.sqrt(2.0),0)
	if b % 2.0 == 0:
		# "b" is even
		a = multiplier - a
	return round(a,0) + round(b,0)*multiplier
#
def ramachandran_number_expand(z):
	z = float(z)
	x = (numpy.sqrt(2.0)*np.mod(z,multiplier)/rho_scaling     - bound + numpy.sqrt(2.0)*np.floor(z / multiplier)/rho_scaling - bound )/2.0
	y = (numpy.sqrt(2.0)*np.floor(z / multiplier)/rho_scaling - bound - numpy.sqrt(2.0)*np.mod(z, multiplier)/rho_scaling + bound)/2.0
	return x,y
```


