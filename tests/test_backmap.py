#!/usr/bin/env python
# -*- coding: utf-8 -*-

#import pytest
from unittest import TestCase


import matplotlib.pyplot as plt

import backmap
import backmap as bm

__author__ = "ranjanmannige"
__copyright__ = "Ranjan Mannige"
__license__ = "MIT"


def test_ramachandran_number():
	assert backmap.R( 180, 180) == 1.0
	assert backmap.R(-180,-180) == 0.0
	assert backmap.R(   0,   0) == 0.5

def test_load_pdb():
	# Set pdb name 
	pdbfn = 'tests/pdbs/1mba.pdb'
	# READ PDB in the form of a matrix with columns
	data = bm.read_pdb(pdbfn)
	assert set(data[1:,0]) == set([1])
	assert len(data[1:,0]) == 146
#

def test_plot_graph(): 
	# Set pdb name 
	pdbfn = 'tests/pdbs/1mba.pdb'
		
	# READ PDB in the form of a matrix with columns
	data = bm.read_pdb(pdbfn)
	# <data> is a 2d array with four columns
	# ['model','chain','resid','R'] ... first row is a header, with values that follow

	# setting the name of the colormap
	cmap = "Chirality"

	# Getting only those values for the particular chain 
	grouped_data = bm.group_data_by(data,group_by = 'chain', columns_to_return=['model','resid','R'])
	for chain in list(grouped_data.keys()):
		print(r'\t',chain)
		# Getting the X,Y,Z values for each entry
		models, residues, Rs = grouped_data[chain]
		
		# Finally, creating (but not showing) the graph 
		response = False
		response = bm.draw_xyz(X = models  ,      Y = residues  ,     Z = Rs
				   ,xlabel =r'Frame #', ylabel =r"Residue #",zlabel =r'$\mathcal{R}$'
				   ,  cmap = cmap    ,  title = cmap      , vmin=0, vmax=1)
		#
		assert response
		#
		# Now, we display the graph:
		#plt.show() # ... one can also use plt.savefig() to save to file
	#
