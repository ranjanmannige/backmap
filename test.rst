<img src="./manuscript/manuscript/figures/banner.png" width="75%">
	

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
import backmap as bm
print bm.R(phi=0,psi=0)
```
See the manuscript more information regarding module usage.

## Stand-alone Usage 

After installation, the following commands produce a variety of graphs (exampled below)

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


### Additional tags
```
-h       -     Prints a help file
-signed  -     Use the signed version of the ramachandran number
```

## Example output I: a structurally stable protein ([1XQQ](https://www.rcsb.org/structure/1XQQ))
The histogram below shows within the protein the presence of both helices (at R \~ 0.34) and sheets (at R \~ 0.52). The panel next to the histogram plot describes the various conformations within the ensemble.
 
<img src="./manuscript/manuscript/automated_figures/pdb_1xqq_his.png" width="25%">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="./manuscript/manuscript/figures/1xqq.png" width="15%">

Additionally, the per-residue plots below (colored by two different metrics or CMAPs), show that most of the protein backbone remains relatively stable (e.g., few fluctuations are evident over the frame \#)

<img src="./manuscript/manuscript/automated_figures/pdb_1xqq_r_Chirality.png" width="25%"><img src="./manuscript/manuscript/automated_figures/pdb_1xqq_r_SecondaryStructure.png" width="25%">

The following plots describe a) the extent towards which a single residue's state has deviated from the first frame (left), and 
b) the extent towards which a single residues state has deviated from the state in the previous frame (right). Both these graphs, 
as expected from the graphs above, show that this protein is relatively conformationally stable.

<img src="./manuscript/manuscript/automated_figures/pdb_1xqq_rmsd.png" width="25%"><img src="./manuscript/manuscript/automated_figures/pdb_1xqq_rmsf.png" width="25%">

## Example output II: an intrinsically disordered protein ([2FFT](https://www.rcsb.org/structure/2FFT))

As compared to the conformationally stable protein above, an intrinsically disordered protein [2FFT](https://www.rcsb.org/structure/2FFT)
is much more flexible. This is especially evident in the fact that each frame within the histogram graph displays a diverse range of R.

<img src="./manuscript/manuscript/automated_figures/pdb_2fft_his.png" width="25%"><img src="./manuscript/manuscript/figures/2fft.png" width="25%">

Interestingly, the graphs below show that while the conformational state of almost every residue dramatically fluctuates, 
except residues 15 through 25, which is the only stable region of the protein. This trend would be hard to recognize by simply looking at the structure.

<img src="./manuscript/manuscript/automated_figures/pdb_2fft_r_Chirality.png" width="25%"><img src="./manuscript/manuscript/automated_figures/pdb_2fft_r_SecondaryStructure.png" width="25%">

The stable region (residues 15 through 25) is also evident when looking at fluctuations.

<img src="./manuscript/manuscript/automated_figures/pdb_2fft_rmsd.png" width="25%"><img src="./manuscript/manuscript/automated_figures/pdb_2fft_rmsf.png" width="25%">

# Publications

The Ramachandran number concept is discussed in the following manuscripts (this tool is discussed in the first reference):

1. Mannige (2018) "The Backmap Python Module: How a Simpler Ramachandran Number Can Simplify the Life of a Protein Simulator" Manuscript Prepared. Preprint available 
the [manuscript/manuscript](manuscript/manuscript/plotmap.pdf) subdirectory of this repo.
2. Mannige, Kundu, Whitelam (2016) "The Ramachandran Number: An Order Parameter for Protein Geometry" PLoS ONE 11(8): e0160023. 
Full Text: https://doi.org/10.1371/journal.pone.0160023
