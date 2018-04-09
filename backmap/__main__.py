from __future__ import division, print_function, absolute_import
import sys
from backmap import *


def main():
	if not "-pdb" in sys.argv:
		if "-h" in sys.argv or "-help" in sys.argv or "--help" in sys.argv:
			pass
		else:
			print("Must provide '-pdb' parameter. Exiting.")
			exit(0)
	show = False
	target_dir = False
	for i in range(len(sys.argv)):
		if sys.argv[i] == "-rmsd":
			showrmsd = 1
		if sys.argv[i] == "-show":
			show = True
		if sys.argv[i] == "-signed":
			print("Using the R number with range [-1,1]")
			signed = 1
			rrange = [-1,1]
		if sys.argv[i] == "-ss":
			colortype = "SecondaryStructure" # default: chirality
		if sys.argv[i] == "-h" or sys.argv[i] == "-help" or sys.argv[i] == "--help":
			print(helpme)
			exit(1)
		if sys.argv[i] == "-pdb":
			if len(sys.argv) <= i+1:
				print(helpme)
				print("MUST PROVIDE PDB NAME.")
				exit(0)
			else:
				pdbfn = str(sys.argv[i+1])
				print("# pdbfn set to:",pdbfn)
		if sys.argv[i] == "-target":
			if len(sys.argv) <= i+1:
				print(helpme)
				print("MUST PROVIDE TARGET DIR.")
				exit(0)
			else:
				target_dir = str(sys.argv[i+1])
				if os.path.isdir(target_dir):
					pass
				else:
					print('SPECIFIED TARGET DIR (%s) DOES NOT EXIST' %(target_dir))
				print("# target directory set to:",target_dir)
		elif sys.argv[i] == "-bins":
			if len(sys.argv) <= i+1:
				helpme
				print("When using '-bins', you must provide bin number. Exiting.")
				exit(0)
			else:
				if not sys.argv[i+1].isdigit():
					print(helpme)
					print("The -bin parameter must be a positive integer (provided: "+str(sys.argv[i+1])+") Exiting.")
					exit(0)
				else:
					bins = int(sys.argv[i+1])
					print("# bins set to:",bins)
					if bins == 0:
						print(helpme)
						print("Must have greater than 0 bins. Exiting.")
						exit(0)
	
	colormap_name = colortype
	if signed:
		colormap_name = colortype+'FourColor'
	print("Using color map name:",colormap_name)
	rcode_cmap = plt.get_cmap(colormap_name)
	#rcode_cmap = plt.get_cmap("deletemeSigned")
	
	pdbfn = os.path.abspath(pdbfn)
	pdbdir = os.path.dirname(pdbfn)
	pdbfilenames = []
	
	if os.path.isfile(pdbfn):
		# then this pathname leads to a FILE
		# ... so keep as is
		pdbfilenames = [pdbfn]
		name = re.split('[\/\.]',pdbfn)[-2]
	elif os.path.isdir(pdbfn):
		pdbdir = pdbfn
		pdbfilenames = sorted(glob.glob(pdbdir+"/*.pdb"))
		name = re.split('[\/\.]',pdbfn)[-1]
	else:
		print(helpme)
		exit("Either filename or directory expected. Exiting.")
	#
	if not target_dir:
		target_dir = pdbdir+"/reports/"
	if not os.path.isdir(target_dir):
		os.makedirs(target_dir)
		
	#NAME = os.path.basename(pdbfilenames[0])[:-len(".pdb")]
	target_base = target_dir.rstrip("/")+"/"
	#
	
	# JUST "CLEVERLY" ARRANGING THE FILENAMES, IF WE HAVE A SET OF FILENAMES RATHER THAN ONE
	# (e.g., pdbfilenames = [something2part1,something1part2,something1part1,something10part1]
	# pdbfilenames.sort() this list to: [something1part1,something1part2,something2part1,something10part1]
	REXP = re.compile( r'\d+' )
	def key_function( s ): return list(map(int, re.findall(REXP, s )))
	pdbfilenames.sort( key=key_function)
	
	print("# Parsing the PDB (structure) data")
	structure = np.array([])
	for pdbfn in pdbfilenames:#[:10]:
		#print(pdbfn)
		# Check if the PDB has no subunit IDs, and then check if segnames exist (right most column)
		# and renaming the subunit IDs alphabetically and then numerically
		
		# READ PDB in the form of a matrix with columns ['model','chain','resid','R']
		latest_structure = read_pdb(pdbfn,signed)
		
		# Getting column indices for each column name
		rx = {} # "rx" for Row indeX
		for c in latest_structure[0,:]:
			rx[c] = list(latest_structure[0,:]).index(c)
		#
		sorted_models        = sorted(list(set(latest_structure[1:,rx['model']])))
		current_model_number = 0
		original_to_new_model_numbers = {}
		for actual_model_number in sorted_models:
			current_model_number                              += 1
			original_to_new_model_numbers[actual_model_number] = current_model_number
		#
		
		# Checking if we already have some structures loaded (then we have to offset the model numbers)
		if len(structure):
			#
			largest_model_number = max(list(structure[1:,rx['model']]))
			for m in list(original_to_new_model_numbers.keys()):
				original_to_new_model_numbers[m] = original_to_new_model_numbers[m] + largest_model_number
			#
		#
		# Resetting model numbers 
		new_model_numbers = [original_to_new_model_numbers[actual_model_number] for actual_model_number in latest_structure[1:,rx['model']]]
		latest_structure[1:,rx['model']] = copy.deepcopy(new_model_numbers)
		#	
		
		if len(structure):
			# Copying as structure
			structure = np.append(structure,copy.deepcopy(latest_structure[1:]), axis=0)
		else:
			# Adding the current model to structure
			structure =  copy.deepcopy(latest_structure)
		#print(structure[:,rx['model']])
		#print(pdbfn)
		#input()
	#
	print("\t...done")
	
	##########################################################################################################
	# READ PDB in the form of a matrix with columns
	#data = read_pdb(pdbfn)
	data = structure
	# Getting only those values for the particular chain 
	grouped_data = group_data_by(data,group_by='chain', 
						columns_to_return=['model','resid','R'])
	
	pdbfn = os.path.split(pdbfn)[-1][:-len('.pdb')]
	
	print(" ---- \t---------")
	print(" TEST \tTEST NAME")
	print(" ---- \t---------")
	
	print(" 1  \tRamachandran number (PDB: %s)"%(name))
	# setting the name of the colormap
	for cmap in ['Greys','SecondaryStructure','Chirality']: #, 'Chirality_r', 'SecondaryStructureHard']:
		# DRAWING A SINGLE GRAPH
		for chain in list(grouped_data.keys()):
			final_name = name
			if len(chain.rstrip()):
				final_name+='-'+str(chain)
			
			# Getting the X,Y,Z values for each entry
			models, residues, Rs = grouped_data[chain]
			
			# Finally, creating (but not showing) the graph 
			plt.clf()	
			draw_xyz(X = models  ,      Y = residues  ,     Z = Rs
					   , xlabel ='Frame #', ylabel ="Residue #",zlabel ='$\mathcal{R}$'
					   , title='Per-residue $\mathcal{R}$; CMAP: %s\nPDB: %s' %(cmap,final_name)
					   ,  cmap = cmap    ,  vmin=0, vmax=1)
			#
			# Now, we display the graph:
			FN = target_base+'pdb_%s_r_%s' %(final_name,cmap)
			write_image(FN)
			print("\tSaved to:",FN)
	#
	# Getting only those values for the particular chain 
	print(" 2.  \tHistogram (PDB: 1xqq)")
	for chain in list(grouped_data.keys()):
		final_name = name
		if len(chain.rstrip()):
			final_name+='-'+str(chain)
		
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
		plt.clf()
		draw_xyz(X = X       ,      Y = Y  ,                Z = Z
		   ,xlabel ='Frame #', ylabel ="$\mathcal{R}$",zlabel ="$P'(\mathcal{R})$:"
			 ,cmap = 'Greys', ylim=[0,1],title='Per-model $\mathcal{R}$-histogram\nPDB: %s'%(final_name))
		plt.yticks(np.arange(0,1.00001,0.2))
		# Now, we display the graph:
		FN = target_base+'pdb_%s_his'%(final_name)
		write_image(FN)
		print("\tSaved to:",FN)
	#
	#
	
	print(" 3.  \tRMSF Test (PDB: {})".format(pdbfn))
	for chain in list(grouped_data.keys()):
		final_name = name
		if len(chain.rstrip()):
			final_name+='-'+str(chain)
		
		# Getting the X,Y,Z values for each entry
		models, residues, Rs = grouped_data[chain]
		
		if len(set(models)) > 1:
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
			plt.clf()
			draw_xyz(X = X       ,      Y = Y  ,                Z = Z
			   ,xlabel ='Frame #', ylabel ="$Residue \#$",zlabel ="$RMSF(\mathcal{R})$:"
				 ,cmap = 'Blues', title='Per-residue RMSF($\mathcal{R}$)\nPDB: %s'%(final_name))
			
			# Now, we display the graph:
			FN = target_base+'pdb_%s_rmsf'%(final_name)
			write_image(FN)
			print("\tSaved to:",FN)
		else:
			print('\tChain "%s" has only one model. Not drawing this graph.' %(chain))
	#
	#
	print(' 4.  \tRMSD Test (PDB: {})'.format(pdbfn))
	for chain in list(grouped_data.keys()):
		final_name = name
		if len(chain.rstrip()):
			final_name+='-'+str(chain)
		
		# Getting the X,Y,Z values for each entry
		models, residues, Rs = grouped_data[chain]
		
		
		if len(set(models)) > 1:
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
			plt.clf()
			draw_xyz(X = X       ,      Y = Y  ,                Z = Z
			   ,xlabel ='Frame #', ylabel ="$Residue \#$",zlabel ="$RMSD(\mathcal{R})$:"
				 ,cmap = 'Reds', title='Per-residue RMSD($\mathcal{R}$)\nPDB: %s'%(final_name))
			#plt.yticks(np.arange(0,1.00001,0.2))
			# Now, we display the graph:
			FN = target_base+'pdb_%s_rmsd'%(final_name)
			write_image(FN)
			print("\tSaved to:",FN)
		else:
			print('\tChain "%s" has only one model. Not drawing this graph.' %(chain))
	#
	##########################################################################################################
	##########################################################################################################
	##########################################################################################################
	##########################################################################################################
	
if __name__ == "__main__":
	main()
