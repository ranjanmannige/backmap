from __future__ import division, print_function, absolute_import

import backmap as bm
import os
import matplotlib.pyplot as plt

print(" ---- \t---------")
print(" TEST \tTEST NAME")
print(" ---- \t---------")
print(" 1.   \tNanosheet Test")
# Set pdb name 
pdbfn = '../tests/pdbs/nanosheet_birth_U7.pdb'
# READ PDB in the form of a matrix with columns
data = bm.read_pdb(pdbfn)

# <data> is a 2d array with four columns
# ['model','chain','resid','R'] ... first row is a header, with values that follow

# setting the name of the colormap
cmap = "Chirality"
# DRAWING A SINGLE GRAPH
# Getting only those values for the particular chain 
grouped_data = bm.group_data_by(data,group_by = 'chain', columns_to_return=['model','resid','R'])
for chain in list(grouped_data.keys()):
	# Getting the X,Y,Z values for each entry
	models, residues, Rs = grouped_data[chain]
	
	# Finally, creating (but not showing) the graph 
	bm.draw_xyz(X = models  ,      Y = residues  ,     Z = Rs
	           , xlabel ='Frame #', ylabel ="Residue #",zlabel ='$\mathcal{R}$'
	           , title=r'1. Nanosheet Test (time vs residue[$\mathcal{R}$])'
		       ,  cmap = cmap    ,  vmin=0, vmax=1)
	#
	# Now, we display the graph:
	plt.show() # ... one can also use plt.savefig() to save to file
#
# DRAWING A SINGLE GRAPH
# Getting only those values for the particular chain 

print(" 2.  \tHistogram Test (PDB: 1xqq)")
import numpy as np # for handling 2d arrays and calculating histograms
pdbfn = '../tests/pdbs/1xqq.pdb'
# READ PDB in the form of a matrix with columns
data = bm.read_pdb(pdbfn)
# Getting only those values for the particular chain 
grouped_data = bm.group_data_by(data,group_by='chain', 
					columns_to_return=['model','resid','R'])
for chain in list(grouped_data.keys()):
	# Getting the X,Y,Z values for each entry
	models, residues, Rs = grouped_data[chain]
	X = []; Y=[]; Z=[]; # Will set X=model, Y=R, Z=P(R)
	# Bundling the three lists into one 2d array
	new_data =  np.array(list(zip(models,residues,Rs)))
	# Getting all R values, model by model
	for m in sorted(set(new_data[:,0])): # column 0 is the model column
		# Getting all Rs for that model #
		current_rs = new_data[np.where(new_data[:,0]==m)][:,2] # column 2 contains R
		# Getting the histogram
		a,b = np.histogram(current_rs,bins=np.arange(0,1.01,0.01))
		max_count = float(np.max(a))
		for i in range(len(a)):
			X.append(m); Y.append((b[i]+b[i+1])/2.0); Z.append(a[i]/float(np.max(a)));
	
	# Finally, creating (but not showing) the graph 
	bm.draw_xyz(X = X       ,      Y = Y  ,                Z = Z
	   ,xlabel ='Frame #', ylabel ="$\mathcal{R}$",zlabel ="$P'(\mathcal{R})$:"
		 ,cmap = 'Greys', ylim=[0,1],title=r'2. Histogram Test (time vs $\mathcal{R}$ [$P(\mathcal{R})$])')
	plt.yticks(np.arange(0,1.00001,0.2))
	# Now, we display the graph:
	#plt.savefig('manuscript/automated_figures/example2.pdf',dpi=200,bbox_inches="tight")
	#print "Saved to:",'manuscript/automated_figures/example.pdf'
	plt.show() # ... one can also use plt.savefig() to save to file

# TRYING OUT VARIOUS CMAPS
# Getting only those values for the particular chain 
grouped_data = bm.group_data_by(data,group_by = 'chain', columns_to_return=['model','resid','R'])
for chain in list(grouped_data.keys()):
	# Getting the X,Y,Z values for each entry
	models, residues, Rs = grouped_data[chain]
	subindex = 0
	for cmap in ['Greys','SecondaryStructure','Chirality']: #, 'Chirality_r', 'SecondaryStructureHard']:
		subindex+=1
		# Finally, creating (but not showing) the graph 
		print(' 3.'+str(subindex)+' \tColor Test (PDB: '+os.path.split(pdbfn)[-1][:-len('.pdb')]+'; CMAP: '+cmap+')')

		bm.draw_xyz(X = models  ,      Y = residues  ,     Z = Rs
		   ,xlabel ='Frame #', ylabel ="Residue #",zlabel ='$\mathcal{R}$'
			 ,cmap = cmap    
			 ,title = '3.'+str(subindex)+' Color Test (CMAP: '+cmap+')'
			 ,vmin=0,vmax=1)
		# Now, we display the graph:
		#plt.show() # ... one can also use plt.savefig() to save to file
		FN = '../manuscript/automated_figures/example_%s.pdf' %(cmap)
		plt.savefig(FN,dpi=200,bbox_inches="tight")
		print "Saved to:",FN
		plt.show()
	#
#



print(" ---- \t---------")
