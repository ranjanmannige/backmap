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
python -m backmap.__init__ -pdb ./pdbs/ProteinDatabankStructureFilename.pdb
python -m backmap.__init__ -pdb /directory/containing/pdbs/
```

The `.__init__` is needed because the main file we are referencing is `backmap/__init__.py`.

### Expected output to the stand alone mode

Three graphs (in both png/raster and pdf/vector format)
```
./pdbs/reports/filename.rcode.pdf/png      (y-axis: residue #; color: R number based on "-signed" and <rcode_cmap>)
./pdbs/reports/filename.rcode.his.pdf/png  (y-axis: Ramachandran number (R); color: frequency of R in model)
./pdbs/reports/filename.rcode.rmsf.pdf/png (y-axis: residue #; color: RMSF in R from the previous model)
./pdbs/reports/filename.rcode.rmsd.pdf/png (y-axis: residue #; color: RMSD in R from the previous model)
```
(the last two files may not be created if only one model exists in the PDB file.)

## Additional tags
```
-h       -     Prints a help file
-ss      -     Color the ramachandran number codes (R-codes) by 
               secondary structure (default: color by chirality and sign)
-signed  -     Use the signed version of the ramachandran number
```

# Publications

The Ramachandran number concept is discussed in the following manuscripts (this tool is discussed in the first reference):

1. Mannige (2018) "The Backmap Python Module: How a Simpler Ramachandran Number Can Simplify the Life of a Protein Simulator" Manuscript Prepared. Preprint available 
the [manuscript/manuscript](manuscript/manuscript/plotmap.pdf) subdirectory of this repo.
2. Mannige, Kundu, Whitelam (2016) "The Ramachandran Number: An Order Parameter for Protein Geometry" PLoS ONE 11(8): e0160023. 
Full Text: https://doi.org/10.1371/journal.pone.0160023
