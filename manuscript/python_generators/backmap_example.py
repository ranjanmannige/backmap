from __future__ import division, print_function, absolute_import

import backmap as bm
import os, copy
import matplotlib.pyplot as plt

import numpy as np # for handling 2d arrays and calculating histograms
pdbfn = '../../tests/pdbs/1xqq.pdb'
pdbfn = '../../tests/pdbs/2k39.pdb'
pdbfn = '../../tests/pdbs/2fft.pdb'
# READ PDB in the form of a matrix with columns
data = bm.read_pdb(pdbfn)
# Getting only those values for the particular chain 
grouped_data = bm.group_data_by(data,group_by='chain', 
					columns_to_return=['model','resid','R'])

pdbfn = os.path.split(pdbfn)[-1][:len('.pdb')]

print(" ---- \t---------")
print(" TEST \tTEST NAME")
print(" ---- \t---------")

print(" 3  \tR Test (PDB: 1xqq)")
# setting the name of the colormap
cmap = "Chirality"
# DRAWING A SINGLE GRAPH
for chain in list(grouped_data.keys()):
	# Getting the X,Y,Z values for each entry
	models, residues, Rs = grouped_data[chain]
	
	# Finally, creating (but not showing) the graph 
	plt.clf()
	bm.draw_xyz(X = models  ,      Y = residues  ,     Z = Rs
	           , xlabel ='Frame #', ylabel ="Residue #",zlabel ='$\mathcal{R}$'
	           , title=r'Per-residue $\mathcal{R}$ (PDB: %s)' %(pdbfn)
		       ,  cmap = cmap    ,  vmin=0, vmax=1)
	#
	# Now, we display the graph:
	FN = '../manuscript/automated_figures/example_%s_r.pdf' %(pdbfn)
	plt.savefig(FN,dpi=200,bbox_inches="tight")
	print("\tSaved to:",FN)
	plt.show() # ... one can also use plt.savefig() to save to file

#
# DRAWING A SINGLE GRAPH
# Getting only those values for the particular chain 
print(" 2.  \tHistogram Test (PDB: 1xqq)")
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
		 ,cmap = 'Greys', ylim=[0,1],title=r'Per-model $\mathcal{R}$-histogram (PDB: %s)'%(pdbfn))
	plt.yticks(np.arange(0,1.00001,0.2))
	# Now, we display the graph:
	FN = '../manuscript/automated_figures/example_%s_his.pdf'%(pdbfn)
	plt.savefig(FN,dpi=200,bbox_inches="tight")
	print("\tSaved to:",FN)
	plt.show() # ... one can also use plt.savefig() to save to file
#
#
print(" 3.  \tRMSF Test (PDB: {})".format(pdbfn))
for chain in list(grouped_data.keys()):
	# Getting the X,Y,Z values for each entry
	models, residues, Rs = grouped_data[chain]
	X = []; Y=[]; Z=[]; # Will set X=model, Y=R, Z=P(R)
	# Bundling the three lists into one 2d array
	new_data =  np.array(list(zip(models,residues,Rs)))
	
	reference_model_number = sorted(set(models))[0]
	
	reference_data = new_data[new_data[:,0]==reference_model_number]
	
	final_data = []
	sorted_models = sorted(set(models))
	for mx in range(1,len(sorted_models)):
		m1 = sorted_models[mx-1]
		m2 = sorted_models[mx]
		
		current_model = new_data[new_data[:,0]==m2]
		current_model[:,2] = np.abs(current_model[:,2] - new_data[new_data[:,0]==m1][:,2])
		if not len(final_data):
			final_data = copy.deepcopy(current_model)
		else:
			final_data = np.append(final_data,current_model,axis=0)
		#
		
	X = final_data[:,0]; 
	Y = final_data[:,1]; 
	Z = final_data[:,2]; 
	
	
	# Finally, creating (but not showing) the graph 
	bm.draw_xyz(X = X       ,      Y = Y  ,                Z = Z
	   ,xlabel ='Frame #', ylabel ="$Residue \#$",zlabel ="$RMSF(\mathcal{R})$:"
		 ,cmap = 'Blues', title=r'Per-residue $\mathcal{R}$-RMSF (PDB: %s)'%(pdbfn))
	#plt.yticks(np.arange(0,1.00001,0.2))
	# Now, we display the graph:
	FN = '../manuscript/automated_figures/example_%s_rmsf.pdf'%(pdbfn)
	plt.savefig(FN,dpi=200,bbox_inches="tight")
	print("\tSaved to:",FN)
	plt.show() # ... one can also use plt.savefig() to save to file
#
#
print(' 4.  \tRMSD Test (PDB: {})'.format(pdbfn))
for chain in list(grouped_data.keys()):
	# Getting the X,Y,Z values for each entry
	models, residues, Rs = grouped_data[chain]
	X = []; Y=[]; Z=[]; # Will set X=model, Y=R, Z=P(R)
	# Bundling the three lists into one 2d array
	new_data =  np.array(list(zip(models,residues,Rs)))
	
	reference_model_number = sorted(set(models))[0]
	
	reference_data = new_data[new_data[:,0]==reference_model_number]
	
	final_data = []
	for m in sorted(set(models)):
		current_model = new_data[new_data[:,0]==m]
		current_model[:,2] = np.abs(current_model[:,2] - reference_data[:,2])
		if not len(final_data):
			final_data = copy.deepcopy(current_model)
		else:
			final_data = np.append(final_data,current_model,axis=0)
		#
		
	X = final_data[:,0]; 
	Y = final_data[:,1]; 
	Z = final_data[:,2]; 
	
	
	# Finally, creating (but not showing) the graph 
	bm.draw_xyz(X = X       ,      Y = Y  ,                Z = Z
	   ,xlabel ='Frame #', ylabel ="$Residue \#$",zlabel ="$RMSF(\mathcal{R})$:"
		 ,cmap = 'Reds', title=r'Per-residue $\mathcal{R}$-RMSF (PDB: %s)'%(pdbfn))
	#plt.yticks(np.arange(0,1.00001,0.2))
	# Now, we display the graph:
	FN = '../manuscript/automated_figures/example_%s_rmsd.pdf'%(pdbfn)
	plt.savefig(FN,dpi=200,bbox_inches="tight")
	print("\tSaved to:",FN)
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
			 ,title = 'CMAP: '+cmap+' (PDB: %s)'%(pdbfn)
			 ,vmin=0,vmax=1)
		# Now, we display the graph:
		#plt.show() # ... one can also use plt.savefig() to save to file
		FN = '../manuscript/automated_figures/example_%s_%s.pdf' %(pdbfn,cmap)
		plt.savefig(FN,dpi=200,bbox_inches="tight")
		print("\tSaved to:",FN)
		plt.show()
	#
#



print(" ---- \t---------")
