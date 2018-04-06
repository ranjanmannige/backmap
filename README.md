```
  ____             _    __  __          _____  
 |  _ \           | |  |  \/  |   /\   |  __ \ 
 | |_) | __ _  ___| | _| \  / |  /  \  | |__) |
 |  _ < / _` |/ __| |/ / |\/| | / /\ \ |  ___/ 
 | |_) | (_| | (__|   <| |  | |/ ____ \| |     
 |____/ \__,_|\___|_|\_\_|  |_/_/    \_\_|     
                       (Multi-angle Picture)                                             
```

This tool provides easily readable "pictures" of protein conformations, 
ensembles, and trajectories saved as either a combined protein databank 
(PDB) structure file, or a directory of such files, and produces graphs.

# Installation
```
pip install backmap
```

# RE-installation
```
pip install -I backmap
```

# Usage

## Module Usage 

```python
import backmap
print backmap.R(phi=0,psi=0)
```
See the manuscript for more information regarding uses.

## Stand-alone Usage 

```
python plotmap.py -pdb ProteinDatabankStructureFilename.pdb
python plotmap.py -pdb /directory/containing/pdbs/
```

## Output 

The x-axis always represents the models/structures listed in the PDB: 

```
filename.rcode.eps      (y-axis: residue #; color: R number based on "-signed" and <rcode_cmap>)
filename.rcode.his.eps  (y-axis: Ramachandran number (R); color: frequency of R in model)
filename.rcode.rmsf.eps (y-axis: residue #; color: RMSF in R from the previous model)
```

## Additional tags
```
-h       -     Prints this message
-ss      -     Color the ramachandran number codes (R-codes) by 
               secondary structure (default: color by chirality and sign)
-signed  -     Use the signed version of the ramachandran number
-rmsd    -     Also producee "filename.rcode.rmsd.eps"
               (y-axis: residue #; color: RMSD in R from first model)
```

# Publications

The Ramachandran number concept is discussed in the following manuscripts (this tool is discussed in the first reference):

1. Mannige (2018) "A simpler Ramachandran number can simplify the life of a protein simulator" Manuscript Prepared/Submitted
2. Mannige, Kundu, Whitelam (2016) "The Ramachandran number: an order parameter for protein geometry" PLoS ONE 11(8): e0160023. 
Full Text: https://doi.org/10.1371/journal.pone.0160023
