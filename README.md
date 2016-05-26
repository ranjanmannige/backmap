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
```
python plotmap.py -pdb ProteinDatabankStructureFilename.pdb
python plotmap.py -pdb /directory/containing/pdbs/
```

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

```python
# PARAMTERS FOR THE RAMACHANDRAN NUMBER
# The ranges for phi and psi are [-180,180]. 
# Any other value will be garbled (so, remember to 
# convert your angles so that it fits this range.
bound = 360.0 # This does not chang
rho_scaling = 10.0 # This is sigma in the manuscript (Mannige, Kundu, Whitelam, 2016)
multiplier = int(round(rho_scaling*(bound*(2.0**0.5)),0)) # For internal reference

def raw_ramachandran_number_collapse(x,y):
	x = float(x)
	y = float(y)
	a = round(rho_scaling*(x-y + bound)/numpy.sqrt(2.0),0)
	b = round(rho_scaling*(x+y + bound)/numpy.sqrt(2.0),0)
	return round(a,0) + round(b,0)*multiplier
#
def raw_ramachandran_number_expand(z):
	z = float(z)
	x = (numpy.sqrt(2.0)*np.mod(z,multiplier)/rho_scaling     - bound + numpy.sqrt(2.0)*np.floor(z / multiplier)/rho_scaling - bound )/2.0
	y = (numpy.sqrt(2.0)*np.floor(z / multiplier)/rho_scaling - bound - numpy.sqrt(2.0)*np.mod(z, multiplier)/rho_scaling + bound)/2.0
	return x,y

phi = -60
psi = -60
# To get the raw Ramachandran number:
raw_R = raw_ramachandran_number_collapse(phi,psi)
# To get the (much more useful) normalized Ramachandran number:
# First getting the lowest and highest possible unnormalized R numbers
raw_R_min = raw_ramachandran_number_collapse(-180,-180)
raw_R_max = raw_ramachandran_number_collapse(180,180)
# Finally, the normalized R number is ...
R = float(raw_R - raw_R_min)/(raw_R_max-raw_R_min)
```


