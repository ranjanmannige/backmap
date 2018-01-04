import urllib2
import os,sys,glob,copy,re
import matplotlib as mpl
mpl.use('agg') # this is to be able to use .get_window_extent() without an error (before drawing a figure)
import matplotlib.pyplot as plt 
from matplotlib import colors
import numpy as np
import numpy,math
import pandas as pd
import seaborn as sns
import itertools, string
# LOCAL IMPORTS
sys.path.insert(0, "./local_imports/") # for the local imports
import locallib

# The figures that we will be creating in using script
figname1 = "manuscript/automated_figures/fig_ramachandran_plots_vs_numbers.pdf"
figname2 = "manuscript/automated_figures/fig_ramachandran_numbers_are_useful1.pdf"	

dryrun   = 0 # For testing. If 1, dont draw the Ramachandran plots (plt.contour() or plt.contourf() plots are computer intensive!)
prunerun = 0 # For testing. Work on a smaller dataset (a fraction of the original size; normally 1%).

show_graphs = 0        # just uses a simple command line that 
pdf_viewer  = "evince" # This is the viewer of choice for the author (based-Linux).

sns.set_style("whitegrid")
panel_letters = itertools.cycle(list(string.ascii_lowercase))

# SUPER IMPORTANT FOR ASTHETICS
black_color = np.ones(3)*0.5  # Darker is closer to 0; lighter is closer to 1
text_color  = np.ones(3)*0.35 # 

# CHECK ALL OTHER PARAMETERS: print plt.rcParams
sns.set_style("whitegrid")
plt.rc_context({'axes.labelcolor':text_color,'axes.linewidth': 1.0,'legend.edgecolor': black_color,
                    'lines.color':black_color,'lines.linewidth':0.5,'text.color':text_color,'axes.edgecolor':black_color,
                    'xtick.color':text_color,'ytick.color':text_color})

# Shows all other parameters to muck about with:
#print plt.rcParams

# For using latex commands in labels, etc.
mpl.rc('text',usetex=True)
mpl.rc('text.latex', preamble='\usepackage{color}')

pgf_with_latex = {
    "text.usetex": True,            # use LaTeX to write all text
    "pgf.rcfonts": False,           # Ignore Matplotlibrc
    "pgf.preamble": [
        r'\usepackage{color}'     # xcolor for colours
    ]
}
mpl.rcParams.update(pgf_with_latex)


motif_database_filename = 'local_imports/data_sequence_motifs.csv'
structure_csv_file = "local_imports/data_pdb.csv"

default_fontsize = 22

aa_three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


sscode_to_latex_code = {
	"h" :r'$\alpha$ $\mathrm{helix}$',
	"h'":r"$\alpha_\mathrm{L}$ $\mathrm{helix}$",
	"g" :r"$3_{10}$ $\mathrm{helix}$",
	"e" :r"$\beta$ $\mathrm{strand}$",
	"p" :"$\mathrm{ppII}$"
	}


number_type_to_label = {
	"R"     :r'$\mathcal{R}$',
	"d'"    :r"$d$",
	"theta" :r"$\theta$"
	}

sscode_to_latex_code_short = {
	"h" :r'$\alpha$',
	"h'":r"$\alpha_\mathrm{L}$",
	"g" :r"$3_{10}$",
	"e" :r"$\beta$",
	"p" :"$\mathrm{ppII}$"
	}

ss_codes={}
#code taken from  S. M. King and W. C. Johnson (1999)
ss_codes['xtlsstr']  = {
	'H': 0, # Alpha -helix 
	'G': 1, # 3.10 - helix 
	'g': 2, # 3.10 - helix C-cap 
	'T': 3, # turns (hydrogen-bonded)
	'N': 4, # turns -'like' (unhydrogen-bonded) 
	'P': 5, # polyproline II type (3-1) helix 
	'p': 6, # polyproline II type (3-1) helix 
	'-': 7, # coils 
	'E': 8, # extended strand
	'U': 9, # not found
}
# code taken from  Cubellis et al., 2005
ss_codes['segno']  = {
	'H':0,# Alpha -helix         
	'G':1,# 3.10 -helix         
	'I':2,# pi -helix         
	'P':3,# Polyproline II  
	'O':4,# coil
	'E':5,# extended strand (in beta-sheet) 
	'b':6,# extended strand (alone)
	'U':7,# not found
}

# Normally, we use only two contour lines for 2d contour plots (such as those used in Ramchandran plots)
# The commented out lines help choose the right two colors for the plot.
'''
for color_type in ['deep','muted','pastel','bright','dark']:
	current_palette = sns.color_palette(color_type)
	sns.palplot(current_palette)
	print color_type
	plt.show()
'''
# After some exploration, the following color set was found to be sufficient enough to see differences from afar
ctype = 'deep' #or 'muted' 
two_tone_colors = [sns.color_palette('pastel')[2],sns.color_palette('dark')[2]]#sns.color_palette('pastel')[0]]
#two_tone_colors.append(black_color)
# Uncomment to see the two colors side-by-side
#sns.palplot(two_tone_colors); plt.show();exit();

chain_index = -1
if 1:
	if not os.path.isfile(motif_database_filename):
		# Checking if the main csv file has been compiled (you do this by reading in PDBs)
		if not os.path.isfile(structure_csv_file):
			columns = ['chain_index','pdf_name','chainid','resid','resno','phi','psi','omega','d','theta','R']
			directory = os.path.split(structure_csv_file)[0]	
			if not os.path.isdir(directory):
				os.makedirs(directory)
	
			print "Populating the file '%s' with the values %s"%(structure_csv_file,str(columns))
			outputf = locallib.open_file(structure_csv_file,'w')
			outputf.write(",".join(columns)+"\n")
		
			pdb_database_dir = "/home/ranjan/Desktop/old/ensembles/pdbstyle-2.06"
			files = glob.glob(pdb_database_dir+"/*/*.ent")
		
			# Going through each PDB file
			for fn in files:#[:100]:
				chain_index += 1
			
				if chain_index % 100==0:
					print chain_index,"   of",len(files), "\t("+str(round(100.0*float(chain_index)/len(files),2))+"%)"
			
				pdf_name = ".".join(os.path.split(fn)[1].split(".")[:-1])
			
				model_to_chain_to_resno_atom_to_vals = locallib.read_pdb_manual(fn)
				for model in sorted(model_to_chain_to_resno_atom_to_vals.keys()):
					for chain in sorted(model_to_chain_to_resno_atom_to_vals[model].keys()):
						# We study chains 1x1x1 (Cloud Cult reference)
					
						# Getting the relationship between residue number and various values
						for resno in sorted(model_to_chain_to_resno_atom_to_vals[model][chain].keys()):
							vals = model_to_chain_to_resno_atom_to_vals[model][chain][resno]
							if 'r' in vals:
								aa_three_letter = vals['resname'].upper()
								if aa_three_letter in aa_three_to_one:
									aa_single_letter = aa_three_to_one[aa_three_letter]
									'''
									'resname'
									'phi'
									'psi'
									'ca'
									'c'
									'n'
									'r'
									'd'
									'theta'
									'''
									resnumber = int(resno)
									s = "%d,%s,%s,%s,%d,%f,%f,%f,%f,%f,%f\n" %(chain_index,pdf_name,chain,aa_single_letter,int(resno),vals['phi'],vals['psi'],vals['omega'],vals['d'],vals['theta'],vals['r'])
									#     |  |  |  |  |  |  |  |  | |  |
									#     |  |  |  |  |  |  |  |  | |  R
									#     |  |  |  |  |  |  |  |  | theta
									#     |  |  |  |  |  |  |  |  d
									#     |  |  |  |  |  |  |  omega
									#     |  |  |  |  |  |  psi
									#     |  |  |  |  |  phi
									#     |  |  |  |  resno
									#     |  |  |  aa_single_letter
									#     |  |  chain
									#     |  pdf_name
									#     chain_index
									outputf.write(s)
	
			outputf.close()	
			print "'%s' now contains all the structural information that we need to study sequence-motifs-to-structure relationships" %(structure_csv_file)
		
		
		
		
		print "READING '%s'" %(structure_csv_file)
		df = pd.read_csv(structure_csv_file)#, sep='\s+', header=None, names=ppII_colnames, comment='>')
		print "\t...done"
			
		sequence_to_data = {}
		
		columns  = ['window','index_of_interest','motif','phi','psi','omega','d','theta','R'] 
		motif_df = pd.DataFrame(columns=columns)
	
		# ADD AN "IF NOT EXISTS" STATEMENT HERE
		outputf = locallib.open_file(motif_database_filename,'w')
		outputf.write(",".join(columns)+"\n")
		
		chains = sorted(set(df['chain_index']))
		len_chains = len(chains)
		len_chains_minus_one = len_chains - 1
		data_items = []
	
		resnumber_key = "resno"
		
	
		dataframe_column_names = list(df.columns.values)
		dataframe_column_names.remove(resnumber_key)
		column_name_to_index = {}
		for val in dataframe_column_names:
			column_name_to_index[val] = dataframe_column_names.index(val)
	
		for chain in chains:
		
			cdf = df[(df['chain_index']==chain)]		
		
			chain_index = chains.index(chain)
			if chain_index%1==0 or chain_index == len_chains_minus_one:
				sys.stdout.write("\r\tCompleted %d of %d     (%d percent)\t\t" %(chain_index,len_chains,int(100.*float(chain_index)/len_chains)))
				sys.stdout.flush()
			
		
			resnumbers = cdf[resnumber_key].tolist()
			cdf_dict   = cdf.set_index(resnumber_key).T.to_dict('list')
		
			# Now, we go through the residues in order to pick up contiguous motifs 
			# of length <window>.
			for window in [1,2]:
				indices_of_interest = range(window)
				for rindex in range(len(resnumbers)):
					#rindex = resnumbers.index(resno)
					if rindex+window-1 < len(resnumbers):
						# If we want to search for <windows> length windows of contiguous sequence, 
						# then <resnumbers[rindex]+window-1> is the last (expected) residue number within the
						# sequence and <resnumbers[rindex+window-1]> is the ACTUAL residue number within the 
						# sequence. I.e., if both the expected and the actual numbers are the same, then we
						# can assume that we can proceed to collect that sequence motif <seq>.
						seq = ""
						if resnumbers[rindex]+window-1 == resnumbers[rindex+window-1]:
							for i in range(window):
								resno = resnumbers[rindex+i]
								aa    = cdf_dict[resno][column_name_to_index['resid']]
								seq  += aa
						if len(seq) == window:
							for index_of_interest in indices_of_interest:
								resno = resnumbers[rindex+index_of_interest]
								statistics = cdf_dict[resno][column_name_to_index['phi']:]
								s = "%d,%d,%s,%f,%f,%f,%f,%f,%f\n" % tuple([window,index_of_interest,seq]+list(statistics))
								#     |  |  |  |  |  |  |  |  |
								#     |  |  |  |  |  |  |  |  R
								#     |  |  |  |  |  |  |  theta
								#     |  |  |  |  |  |  d
								#     |  |  |  |  |  omega
								#     |  |  |  |  psi
								#     |  |  |  phi
								#     |  |  motif
								#     |  index_of_interest
								#     window
								outputf.write(s)
								#motif_df.loc[motif_df.shape[0]] = [window,index_of_interest,seq]+list(statistics)
						#
					#
				#
			#
		#
		print 'DONE'
	#
	
	#".pruned.csv" is for testing only (created by the commented-out lines below)
	if dryrun or prunerun:
		motif_database_filename+=".pruned.csv"
	
	print "READING '%s'" %(motif_database_filename) 
	df = pd.read_csv(motif_database_filename)#, sep='\s+', header=None, names=ppII_colnames, comment='>')
	# Just removing any row with one column that is not valid (nan)
	df = df.dropna(axis='index', how='any', inplace=False)
	#df = df.sample(int(len(df)*0.5))
	print "\t...done"
	
	# Reframing the d values to be all positive
	'''  min_d             0             max_d
	     |---x-------------|-----------------|---x'                         
	'''
	tmp_ds = df['d'].unique(); max_d = max(tmp_ds); min_d = min(tmp_ds); 
	df.ix[df['d'] <= 0.0, 'd'] = (df['d'] + max_d - min_d) #.where(df['d'] < 0.0)
	
	
	'''
	# MAKE A PRUNED VERSION OF THE DATAFRAME
	# Prunning during the testing stage (makes testing faster)
	df = df.sample(int(len(df)*0.01))
	df.to_csv(motif_database_filename+".pruned.csv")
	'''
	
	ss_to_ave = {}
	ss_to_std = {}
	
	# These statistics are obtained from the commented-out code immediately following these dictionaries
	# SS ANNOTATION SYSTEM: 'segno'
	ss_to_ave = {'H': {'theta': 94.869608333674606, 'phi': -64.045325858765509, 'psi': -39.80643438917587, 'd': 1.6936429390055281, 'R': 0.3557666271151626},
	             'E': {'theta': 190.71906314633824, 'phi': -111.69078623558029, 'psi': 127.98347332595094, 'd': 3.1700476058602014, 'R': 0.52258078429607835},
	             'G': {'theta': 102.30010376793231, 'phi': -68.970919767340533, 'psi': -22.882810474032762, 'd': 1.9829456611133074, 'R': 0.37242525349584371},
	             'P': {'theta': 226.29783699828431, 'phi': -77.181136626842829, 'psi': 138.60182926220506, 'd': 2.8759942748528973, 'R': 0.58526739603127531}
	             }
	ss_to_std = {'H': {'theta': 7.8677226180508546, 'phi': 9.6358076878084109, 'psi': 9.3193974578006511, 'd': 0.21356227985531803, 'R': 0.014295218204378299},
	             'E': {'theta': 24.185214919987217, 'phi': 40.984539304168315, 'psi': 52.003491023970412, 'd': 1.0580618540080533, 'R': 0.077440041491631625},
	             'G': {'theta': 9.8754628623228484, 'phi': 15.815100855679685, 'psi': 16.029652384646734, 'd': 0.24446253093295564, 'R': 0.01757487747371899},
	             'P': {'theta': 27.006717042530056, 'phi': 21.678391427421019, 'psi': 45.895255639479807, 'd': 0.86053939688978853, 'R': 0.071576188296957083}
	             }
	'''
	# SS ANNOTATION SYSTEM: 'xtlsstr'
	ss_to_ave = {'H': {'theta': 100.1701356411558, 'phi': -66.475132717792803, 'psi': -30.859036062946029, 'd': 1.7329957585062672, 'R': 0.36481609158670159}, 
	             'E': {'theta': 182.97932180475496, 'phi': -110.06812923872542, 'psi': 115.91870081407123, 'd': 3.090403596661679, 'R': 0.5080811279970725}, 
	             'G': {'theta': 111.21570651919835, 'phi': -66.757063858495968, 'psi': -15.458711621444696, 'd': 1.9535556636367633, 'R': 0.38581050563673974},
	             'P': {'theta': 218.35222256754639, 'phi': -81.665655635341849, 'psi': 131.53028062976674, 'd': 2.8785413145094005, 'R': 0.56921745676971269}}

	ss_to_std = {'H': {'theta': 28.651152050754906, 'phi': 20.263899078995486, 'psi': 40.436761988031222, 'd': 0.52447988678207591, 'R': 0.055180889874498901},
	             'E': {'theta': 32.678268223307988, 'phi': 46.283306311466561, 'psi': 64.482190153771128, 'd': 1.0842447928625454, 'R': 0.083276093319786013},
	             'G': {'theta': 36.049839915820613, 'phi': 26.164045381471915, 'psi': 51.16639117566875, 'd': 0.65231838846283563, 'R': 0.069558167458037279},
	             'P': {'theta': 33.410514502854959, 'phi': 27.522062193031577, 'psi': 55.418922097118923, 'd': 1.0032479119053197, 'R': 0.08392549028152399}}
	'''
	
	if not len(ss_to_ave) or not len(ss_to_std):
		# OBTAINING AND PRINTING DICTIONARIES UP THERE
		annotation_system = 'segno' # or 'xtlsstr'
		master_ss_file = "local_imports/data_ss.csv" #
		# LOADING
		ssdf           = pd.read_csv(master_ss_file)
		# CLEANING
		ssdf           = ssdf.dropna(axis='index', how='any', inplace=False)
		ssdf           = ssdf[(ssdf['d'] != 999)]
	
		# Getting standard deviations for each ss of importance
		ss_to_ave = {}
		ss_to_std = {}
		for ss in ['H','G','P','E']:
			print ss
			ss_to_ave[ss] = {}
			ss_to_std[ss] = {}
		
			tdf = ssdf[(ssdf['segno']==ss_codes['segno'][ss])]
			for key in ['phi','psi','d','theta','R']:
				vals = tdf[key]
				ss_to_ave[ss][key] = np.mean(vals)
				ss_to_std[ss][key] = np.std(vals)
			
				print '\t',key,'\t',np.mean(vals),'\t',np.std(vals)
	
		
		# To prevent having to calculate this again and again, paste these two values above
		print "ss_to_ave = "+str(ss_to_ave)
		print
		print "ss_to_std = "+str(ss_to_std)
	
	# drawing per-aa behavior 
	'''
	Panel A: per-aa behavior in ramachandran plots (5x4)
	 __ __ __ __ __
	|__|__|__|__|__|
	|__|__|__|__|__|
	|__|__|__|__|__|
	|__|__|__|__|__|
	
	Panel B:
	The Ramachandran number version of the same
	  R
	  |
	  |_____AA
	
	'''
	
	textsize =12
	
	max_columns = 10
	max_rows    = 2
	
	
	# VERY IMPORTANT FOR GETTING THE SHAPE OF EACH PANEL CORRECT
	fig = plt.figure(0,figsize=(max_columns*2.1/1.7,max_rows*4.8/1.7))
	# Declaring the various panels
	# The mother frame (into which two different grids are going to be placed)
	whole_gs = mpl.gridspec.GridSpec(2,1,height_ratios=[1,1]) # GridSpec(rows, columns)
	#                                |
	#                                Says that we want two panels one on top of the other
	
	# The two grids (panels)
	number_gs = mpl.gridspec.GridSpecFromSubplotSpec(1       , 6          , subplot_spec=whole_gs[1])
	rama_gs   = mpl.gridspec.GridSpecFromSubplotSpec(max_rows, max_columns, subplot_spec=whole_gs[0])
	#                                                |         |                                  |
	#                                                |         columns                            Says which panel comes firsts
	#                                                rows
	
	# If we were not putting grids within grids, then we can initiae a grid like so:
	#rama_gs = mpl.gridspec.GridSpec(max_rows, max_columns)
	
	# Declaring the various panels is done so:
	# ax = plt.subplot(gs[row_no,col_no])
	
	
	# VERY IMPORTANT FOR GETTING THE SHAPE OF EACH PANEL CORRECT
	max_columns2 = 8
	max_rows2    = 2
	fig2 = plt.figure(1,figsize=(max_columns2*1.2*1.2,max_rows2*1.2*2))
	
	# If we were not putting grids within grids, then we can initiae a grid like so:
	motif_grid = mpl.gridspec.GridSpec(max_rows2, max_columns2+1,height_ratios=[2,3.5])
	
	# Declaring the various panels is done so:
	plt.figure(1)
	figure2_L_sub_plot = plt.subplot(motif_grid[0,0:3])
	figure2_P_sub_plot = plt.subplot(motif_grid[0,3:6])
	figure2_all_sub_plot0 = plt.subplot(motif_grid[1,0:4])
	figure2_all_sub_plot1 = plt.subplot(motif_grid[1,4:8])
	# To play with 
	#plt.tight_layout();figure2_L_sub_plot.scatter(0,0,s=100); figure2_P_sub_plot.scatter(1,1,s=100); figure2_all_sub_plot0.scatter(2,2,s=100); figure2_all_sub_plot1.scatter(2,2,s=100);plt.show();exit();
	
	for window in [1,2]:#sorted(df['window'].unique()):
		for index_of_interest in sorted(df[(df['window']==window)]['index_of_interest'].unique()):
			print "#window =",window
			print "#index_of_interest =",index_of_interest
			
			cdf = df[ (df['window']==window) & (df['index_of_interest']==index_of_interest) ]
			# OR:
			#cdf = df.loc[(df['window']==window) & (df['index_of_interest']==index_of_interest), ['motif','phi','psi','d','theta','R']]
			
			bins   = 80
			number_of_interest = 'R'
			
			# Filling in a local database
			print "\tPopulating data from DataFrame"
			sequence_to_data = {}
			
			counter = 0
			motifs = sorted(cdf['motif'].unique())
			for seq in motifs:
				counter += 1
				if not pd.isnull(seq):
					# This line is very time consuming... is there a better way?
					tcdf = cdf[ (cdf['motif'] == seq) ]
					
					sys.stdout.write("\r\tRetreiving %d of %d motifs    (%d percent)\t\t" %(counter,len(motifs),int(100.*float(counter)/len(motifs))))
					sys.stdout.flush()
					
					if not seq in sequence_to_data:
						sequence_to_data[seq] = {}
					
					for tag in ['phi','psi','d','theta','R']:
						sequence_to_data[seq][tag] = list(tcdf[tag])
					#
				#
			#
			print "\n\t\t...done"
			
			print "\tArranging the sequence motifs"
			# Arranging the sequence motifs
			sequences = []
			if 0: # Make '0' if you just want to sort alphabetically
				r_val_to_seq = {}
				# Getting average values of R to to sort the sequence motifs (next)
				for seq in sequence_to_data.keys():
					R = sequence_to_data[seq][number_of_interest]
					rval = np.average(R)
					if not rval in r_val_to_seq:
						r_val_to_seq[rval]=[]
					r_val_to_seq[rval].append(seq)
				
				# Sorting the sequence motifs
				for rval in sorted(r_val_to_seq.keys()):
					for s in sorted(r_val_to_seq[rval]):
						sequences.append(s)
			else:
				# Just sorting alphabetically
				sequences = sorted(sequence_to_data.keys())
			print '\t\t...done'
			
			# Making labels
			# IF window is 3 and index_of_interest is 1, then masterpattern is "X_X"
			# IF window is 3 and index_of_interest is 0, then masterpattern is "_XX"
			# IF window is 2 and index_of_interest is 1, then masterpattern is "X_"
			# IF window is 1 and index_of_interest is 0, then masterpattern is "_"
			masterpattern = ""
			for i in range(0,index_of_interest):
				masterpattern+="X"
			masterpattern+="_"
			for i in range(index_of_interest+1,window):
				masterpattern+="X"
			#
			
			#Example of patterns when window == 2 and index_of_interest == 0
			#['_A','_C','_D','_E','_F','_G','_H','_I','_K','_L','_M','_N','_P','_Q','_R','_S','_T','_V','_W','_Y']
			#Example of patterns when window == 2 and index_of_interest == 1
			#['A_','C_','D_','E_','F_','G_','H_','I_','K_','L_','M_','N_','P_','Q_','R_','S_','T_','V_','W_','Y_']
			patterns = []
			for seq in sequences:
				patterns.append(seq[0:index_of_interest]+"_"+seq[index_of_interest+1:window])
			patterns = sorted(set(patterns))
			
			# Now, digesting!
			masterx = []
			mastery = []
			masterz = []
			master_horizontal_lines = [ ]
			master_vertical_lines   = [ ]
			masterxlabels           = [ ]
			masterxticks            = [ ]   
			final_index             =  0
			
			# OK, now, we need to start figuring out which panel to draw to
			row    =  0
			column = -1
			
			for pattern in patterns:
				column += 1
				if column >= max_columns:
					column  = 0
					row    += 1
				
				# if pattern == '_A', then get_pattern looks for '.A' and get_residue looks for '(.)A'
				get_pattern = re.compile(pattern.replace("_","."))
				get_residue = re.compile(pattern.replace("_","(.)"))
				
				x = []
				y = []
				z = []
				
				y2 = []
				z2 = []
				current_sequences = []

				# OK, now, we need to start figuring out which panel to draw to (for plt.figure(0))
				rama_row    =  0
				rama_column = -1		
				index = 0
				
				for seq in sequences:#sorted(sequence_to_data.keys()):
					if get_pattern.match(seq):
						# Getting the amino acid of importance (at position <index_of_importance> of the <seq>)
						q          = get_residue.match(seq)
						central_aa = q.group(1)
						
						rama_column += 1
						if rama_column >= max_columns:
							rama_column  = 0
							rama_row    += 1
						
						
						sys.stdout.write("\r\tPATTERN: %s \t Current sequence: %s   (%d of %d)  \t\t" %(pattern,seq,sequences.index(seq)+1,len(sequences)))
						sys.stdout.flush()
						# ==================================================
						# DRAW RAMACHANDRAN PLOTS IF NEEDED
						if window == 1: # draw ramachandran plots
							plt.figure(0)
							with sns.axes_style("darkgrid"):
								# The panel to which the Ramachandran plot will be drawn
							
								ax = plt.subplot(rama_gs[rama_row,rama_column])
							
								#Draw ramachandran plots
								#fig = plt.figure()  # create a figure object
								#ax  = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
							
								if not dryrun:
									X,Y = sequence_to_data[seq]['phi'],sequence_to_data[seq]['psi']
									# Draw individual ramachandran plots
									locallib.plot_ramachandran_histogram(X,Y,
										levels=[0.3,0.6], fill=1,lines=0,label_contours=0, smooth=2, ax=ax,
										linestyles=['solid',None],linewidths=1,cmap=two_tone_colors)
										# |
										# does not matter if lines==0
								if dryrun:
									ax.axis([-180,180,-180,180]); ax.set_xticks(range(-180,181,180)); ax.set_yticks(range(-180,181,180));
								plt.sca(ax);plt.xticks(rotation=45)
								plt.sca(ax);plt.yticks(rotation=45)
								ax.set_xticklabels([])
								ax.set_xlabel("", size=textsize)
								ax.set_yticklabels([])
								ax.set_ylabel("", size=textsize)
							
								#ax.set_title(seq)
								ax.text(90, 90, seq, style='normal', bbox={'facecolor':'white', 'alpha':1, 'pad':2})
							
								if rama_row == max_rows-1:
									# Then we need to draw the X axes elements
									#ax.set_xticks(      range(-180,181,180))
									#ax.set_xticklabels(range(-180,181,180))
									if rama_column == 0:
										ax.set_xlabel('$\phi \longrightarrow$')	
								if rama_column == 0:
									# Then we need to draw the Y axes elements
									#ax.set_yticks(     range(-180,181,180))
									#ax.set_yticklabels(range(-180,181,180))
									if rama_row == max_rows-1:
										ax.set_ylabel('$\psi \longrightarrow$')
								#else:
								#	ax.set_yticks([])
								
								#ax.set_aspect(1)
								#plt.show()
								
								"""
								fn_bar = fn+"_colorbar.eps"
								print "#WRITING TO:",fn_bar
								make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=cmap, fn=fn_bar, yscaling=0.5,xtitle="Key",ytitle="R")
								#
								"""
						# ==================================================
						
						# These are just labels for the X axis
						current_sequences.append(central_aa)
						
						# Getting the distribution 
						current_range = [0.29,0.71]
						if number_of_interest == 'd':
							current_range = [1,4]
						a,b = np.histogram(sequence_to_data[seq][number_of_interest],bins=bins,range=current_range)
						maxa = np.max(a)
						for i in range(len(a)):
							currenty = float(b[i]+b[i+1])/2.0
							currentz = 0.0
							if maxa:
								currentz = float(a[i])/maxa
							#
							x.append(index)
							y.append(currenty)
							z.append(currentz)
							masterx.append(final_index)
							mastery.append(currenty)
							masterz.append(currentz)
						#
						ds   = list(sequence_to_data[seq]['d'])
						a,b  = np.histogram(ds,bins=100,range=[1,4])#3.5])
						maxa = np.max(a)
						for i in range(len(a)):
							currenty = float(b[i]+b[i+1])/2.0
							currentz = 0.0
							if maxa:
								currentz = float(a[i])/maxa
							#x.append(index)
							y2.append(currenty)
							z2.append(currentz)
						#
						
						index+=1
						final_index+=1
					#
				#
				print #Breaking the auto-refreshing line above (see '\r')
				# Deciding where to draw vertical lines to distinguish one set of motifs from another
				master_vertical_lines.append(final_index)
				sorted_x = sorted(set(x))
				
				masterxticks.append(sorted(set(masterx))[-10])
				masterxlabels.append(pattern.replace('_',''))
				
				# This is the sorted list of distinct X elements
				sortedSetX = sorted(set(x))
				# This is half the spacing between any two x-elements
				# (we assume that all elements are equally spaced from their neighbors)
				sortedSetXspacingBy2 = float(sortedSetX[1]-sortedSetX[0])/2.0
				xlines = {}
				xlines[sortedSetX[0]-sortedSetXspacingBy2] = ""
				for xnow in sortedSetX:
					xlines[xnow+sortedSetXspacingBy2] = ""
				
				xlines = []
				
				ylines = {}
				ylines_abbreviated = {}
				for ss in ['E','H',"H'",'P']:
					sscode             = sscode_to_latex_code[ss.lower()]
					sscode_abbreviated = sscode_to_latex_code_short[ss.lower()]
					val = ss_to_ave[ss[0]][number_of_interest]
					if ss[-1] == "'":
						if number_of_interest == 'R':
							ylines[1.0-val]             = sscode
							ylines_abbreviated[1.0-val] = sscode_abbreviated
							
						# else ignore this ss
					else:
						ylines[val] = sscode #.lower()
						ylines_abbreviated[val] = sscode_abbreviated
				
				master_horizontal_lines = ylines_abbreviated # lines indicating secondary structures 
				
				if window == 1:
					plt.figure(0)
					ax = plt.subplot(number_gs[0,0:3])
					current_ticks = np.arange(len(current_sequences))
					
					ytitle = number_of_interest
					if number_of_interest in number_type_to_label:
						ytitle = number_type_to_label[number_of_interest] 
					locallib.make2Dfigure(x,y,z,fn=[],ax=ax,xticks=current_ticks, xlabels=current_sequences,
						              horizontallines=ylines,verticallines=xlines,colorbar=0,cmap=plt.get_cmap('Reds'),
						              ytitle=ytitle,xtitle=None,show=0)#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
					ax.minorticks_off()
					ax.tick_params(axis='y', which='major', direction='out', colors=plt.rcParams['xtick.color'], width=plt.rcParams['axes.linewidth'], length=plt.rcParams['axes.linewidth']*3)#,labelsize=for)
					ax.yaxis.set_ticks_position('left')
					ax.set_xlabel(''); # This is because this label indicates a more complicated pattern
					
					# ---------------------------------------------------------------
					# Putting a little vertical arrow where the proline column exists
					arrowymin,arrowymax = ax.get_ylim(); 
					arrowx = current_ticks[current_sequences.index('P')]
					ax.annotate('', xy    =[arrowx,arrowymin+0.09*(arrowymax-arrowymin)],   xycoords='data', # TO   THE TOP
							xytext=[arrowx,arrowymin+0.01*(arrowymax-arrowymin)], textcoords='data', # FROM THE BOTTOM
							ha="center", va="center",
							arrowprops=dict(arrowstyle="-|>", fc=text_color, color=text_color, linewidth = plt.rcParams['axes.linewidth']*1.5)
							)
					# ---------------------------------------------------------------
					
					
				#figure2_L_sub_plot, figure2_P_sub_plot, figure2_all_sub_plot0, figure2_all_sub_plot1
				if pattern == '_L' or pattern == '_P':
					plt.figure(1)
					
					ax = figure2_P_sub_plot
					ytitle = ''
					if pattern == '_L':
						ytitle = number_of_interest
						if number_of_interest in number_type_to_label:
							ytitle = number_type_to_label[number_of_interest]
						ax = figure2_L_sub_plot
						# Also, ignore all horizontal annotations:
						#for tmpykey in ylines_abbreviated.keys():
						#	ylines_abbreviated[tmpykey] = ' ' # put a space because no space means no annotation line
					
					locallib.make2Dfigure(x,y,z,fn=[],ax=ax,xticks=np.arange(len(current_sequences)), xlabels=current_sequences,
						              horizontallines=ylines_abbreviated,verticallines=xlines,colorbar=0,cmap=plt.get_cmap('Reds'),
						              ytitle=ytitle,xtitle=None,show=0)#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
					
					
					#ax.set_xlabel('X'); # This is because this label indicates a more complicated pattern
					ax.set_title('Motif studied: '+pattern.replace('_',r'\textbf{X}'), loc='right',ha='right'); # This is because this label indicates a more complicated pattern
					# panel letter to the left
					t = ax.set_title(r'('+next(panel_letters)+r')', size=textsize*1.1,loc='left', horizontalalignment='center')
					t.set_y(1.03)                           # Shifting the title up a little
					
					
					# ----------------------------------------------------------------
					# ADDING THE "X=" annotation
					plt.draw() #messes with plt.show()... can only write to file after this
					transf = ax.transData.inverted()
					xticklabels = ax.get_xticklabels()
					bb0 = xticklabels[0].get_window_extent().transformed(transf)
					bb1 = xticklabels[1].get_window_extent().transformed(transf)
					label_x        = bb0.x0 - (bb1.x0 - bb0.x0)/4.
					label_y        = bb0.y0
					label_fontsize = xticklabels[0].get_fontsize()
					plt.text(label_x,label_y,r'\textbf{X} =',ha='right',va='bottom',fontsize=label_fontsize)
					# ----------------------------------------------------------------
					
					# Some options negate the y-ticks. This brings them back and to the left
					ax.tick_params(axis='y', which='major', direction='out', colors=plt.rcParams['xtick.color'], width=plt.rcParams['axes.linewidth'], length=plt.rcParams['axes.linewidth']*4)#,labelsize=for)
					ax.yaxis.set_ticks_position('left')
				
				#ax.title(pattern.replace("_","-"))
				
				'''
				# Drawing some ss 'bands'
				for ss in ['H','E']:
					ss_ymin = ss_to_ave[ss][number_of_interest] - ss_to_std[ss][number_of_interest]/5.0
					ss_ymax = ss_to_ave[ss][number_of_interest] + ss_to_std[ss][number_of_interest]/5.0
					
					xmin,xmax = ax.get_xlim(); ymin,ymax = ax.get_ylim();  # getting axes min max values
					ax.plot([xmin,xmax], [ss_ymin,ss_ymin], c='Gray')
					ax.plot([xmin,xmax], [ss_ymax,ss_ymax], c='Gray')
				'''
				#
				'''
				ax = plt.subplot(number_gs[0,3:5])
				locallib.make2Dfigure(x,y2,z2,fn=[],ax=ax,xticks=np.arange(len(current_sequences)), xlabels=current_sequences,
				                      colorbar=0,cmap=plt.get_cmap('Reds'), #horizontallines=ylines,verticallines=xlines,
				                      ytitle='$d$',xtitle=pattern.replace("_","-"),show=0)#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
				if window == 1: ax.set_xlabel(''); # This is because this label indicates a more complicated pattern
				'''
				
				
				
				#ax.set_aspect(1)
				#ax.scatter(row,column,s=100)
			
			
			if window == 2:
				plt.figure(1)
				sns.axes_style("darkgrid")	
				ax = figure2_all_sub_plot1
				final_master_horizontal_lines = copy.deepcopy(master_horizontal_lines)
				ytitle = ''
				if index_of_interest == 0:
					ytitle = number_of_interest
					if number_of_interest in number_type_to_label:
						ytitle = number_type_to_label[number_of_interest]
					ax = figure2_all_sub_plot0					
				
				#
				locallib.make2Dfigure(masterx,mastery,masterz,fn=[],ax=ax,ytitle=ytitle,
				        xticks=masterxticks,xlabels=masterxlabels,colorbar=0,cmap=plt.get_cmap('Reds'),
					horizontallines=final_master_horizontal_lines,verticallines=master_vertical_lines,
					xtitle=None,show=0)#,title=pattern.replace("_","-"))#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
				
				#ax.set_xlabel('X'); # This is because this label indicates a more complicated pattern
				title = r'Motif studied: X\textbf{Y}'
				if index_of_interest == 1:
					title = r'Motif studied: \textbf{Y}X'
				ax.set_title(title,loc='right',ha='right'); # This is because this label indicates a more complicated pattern
				
				# panel letter to the left
				t = ax.set_title(r'('+next(panel_letters)+r')', size=textsize*1.1,loc='left', horizontalalignment='center')
				t.set_y(1.03)                           # Shifting the title up a little
				
				# ----------------------------------------------------------------
				# ADDING THE "Y=" annotation
				plt.draw() #messes with plt.show()... can only write to file after this
				transf = ax.transData.inverted()
				xticklabels = ax.get_xticklabels()
				bb0 = xticklabels[0].get_window_extent().transformed(transf)
				bb1 = xticklabels[1].get_window_extent().transformed(transf)
				label_x        = bb0.x0 - (bb1.x0 - bb0.x0)/4.
				label_y        = bb0.y0
				label_fontsize = xticklabels[0].get_fontsize()
				plt.text(label_x,label_y,r'\textbf{Y} =',ha='right',va='bottom',fontsize=label_fontsize)
				# ----------------------------------------------------------------
				
				
				# Some options negate the y-ticks. This brings them back and to the left
				ax.tick_params(axis='y', which='major', direction='out', colors=plt.rcParams['xtick.color'], width=plt.rcParams['axes.linewidth'], length=plt.rcParams['axes.linewidth']*4)#,labelsize=for)
				ax.yaxis.set_ticks_position('left')
				'''
				for spine in ax.spines.keys():
				
					print ax.spines[spine]
					raw_input()
					#ax.spines[spine].set_linewidth(10)
					#ax.spines[spine].set_linewidth(10)
				'''
				
			#
			if window == 1:
				plt.figure(0)
				plt.tight_layout();
				plt.subplot(rama_gs[  0,0  ]).annotate('(a)', xy=(.0, 1), xycoords='figure fraction', horizontalalignment='left', verticalalignment='top', fontsize=textsize*1.2)
				plt.subplot(number_gs[0,0:3]).annotate('(b)', xy=(.0, 0.5), xycoords='figure fraction', horizontalalignment='left', verticalalignment='top', fontsize=textsize*1.2)
				print "WRITING TO:",figname1
				plt.savefig(figname1,dpi=170,bbox_inches="tight")
				if show_graphs:
					os.system(pdf_viewer+" "+figname1)
			
			# Drawing connectors and lines
			if window == 2 and index_of_interest == 0:
				plt.figure(1)
				
				# ---------------------------------------------------------------
				# Putting a little vertical arrow where the proline column exists
				arrowymin,arrowymax = ax.get_ylim();
				arrowx = masterxticks[masterxlabels.index('P')]
				ax.annotate('', xy    =[arrowx,arrowymin+0.10*(arrowymax-arrowymin)],   xycoords='data', # TO   THE TOP
						xytext=[arrowx,arrowymin+0.01*(arrowymax-arrowymin)], textcoords='data', # FROM THE BOTTOM
						ha="center", va="center",
						arrowprops=dict(arrowstyle="-|>", fc=text_color, color=text_color, linewidth = plt.rcParams['axes.linewidth']*1.5)
						)
				# ---------------------------------------------------------------
				
				
				# Drwaing the various lines
				for ax1,ax2,letter,armlengths in [[figure2_L_sub_plot,figure2_all_sub_plot0,'L',[14,14]]
				                      ,[figure2_P_sub_plot,figure2_all_sub_plot0,'P',[17,14]]]:
					xmin,xmax = ax1.get_xlim(); ymin,ymax = ax1.get_ylim();  # getting axes min max values
					line1ax1  = [xmin,ymin-(ymax-ymin)*0.16]
					line2ax1  = [xmax,ymin-(ymax-ymin)*0.16]
					
					'''   ____________ __ymax
					     |            |
					     |            |
					     |____________|__ymin                    ax1
					     |            |
					  xmin            xmax
					      \     ______|
					 line1 |   /       line2
					 ______|__|__________ __ ax2_ymax
					|      :  :          |
					|      :  :          |                       ax2
					|      :  :          |
					|______:__:__________|__ ax2_ymin
					 aa_xmin  aa_xmax
					'''
					
					ax2_ymin,ax2_ymax = ax2.get_ylim();
					aa_index          = masterxlabels.index(letter)
					aa_xmin           = master_vertical_lines[aa_index]-20
					aa_xmax           = master_vertical_lines[aa_index]
					
					line1ax2 = [aa_xmin, ax2_ymax]
					line2ax2 = [aa_xmax, ax2_ymax]
					
					#text_A = min_avex + ss_frac*(max_avex-min_avex)
					#text_A = 1.9+next(paddings)
					#
					#point_x = ss_to_avex[ss]
					#point_y = 1.05
					#
					# length of vertical leg emanating from the lower panel
					line1_armA = 17; line2_armA = line1_armA; # 
					# length of vertical leg connecting to the upper panel
					line1_armB = 7; line2_armB = line1_armB;
					
					# Overriding the defaults:
					if armlengths is not None:
						if type(armlengths) is int or type(armlengths) is float:
							line1_armA = armlengths; line2_armA = line1_armA; # 
						if type(armlengths) is tuple or type(armlengths) is list:
							line1_armA = armlengths[0]
							line2_armA = armlengths[1]
						#
					#
					current_arrow_color = np.ones(3)*0.65 # plt.rcParams['lines.color']
					
					angleA =  90
					angleB = -90
					line1_connection_style = "arc,angleA=%f,angleB=%f,armA=%f,armB=%f,rad=0"%(angleA,angleB,line1_armA,line1_armB)
					line2_connection_style = "arc,angleA=%f,angleB=%f,armA=%f,armB=%f,rad=0"%(angleA,angleB,line2_armA,line2_armB)
					
					# line 1
					ax2.annotate('',xy    =line1ax1,   xycoords=ax1.transData,
							xytext=line1ax2, textcoords=ax2.transData, 
							ha="center", va="center",
							arrowprops=dict(arrowstyle="-",connectionstyle=line1_connection_style,color=current_arrow_color)
							)
					
					ax2.annotate('',xy    =line2ax1,   xycoords=ax1.transData,
							xytext=line2ax2, textcoords=ax2.transData, 
							ha="center", va="center",
							arrowprops=dict(arrowstyle="-",connectionstyle=line2_connection_style,color=current_arrow_color)
							)
					
					'''
					motif_grid.update(wspace=1.8,hspace=0.4)
					con = mpl.patches.ConnectionPatch(xyA=line1ax1, xyB=line1ax2, coordsA='data', coordsB='data', axesA=ax1, axesB=ax2, arrowstyle="->")#, shrinkB=5)
					ax2.add_artist(con)
					#con = mpl.patches.ConnectionPatch(xyA=xy, xyB=xy, coordsA='data', coordsB='data', axesA=ax2, axesB=ax1, arrowstyle="-") #, shrinkB=5)
					#ax1.add_artist(con)
					'''
					#plt.show()	
	
	'''
	left   # the left side of the subplots of the figure
	right  # the right side of the subplots of the figure
	bottom # the bottom of the subplots of the figure
	top    # the top of the subplots of the figure
	wspace # the amount of width reserved for blank space between subplots, expressed as a fraction of the average axis width
	hspace # the amount of height reserved for white space between subplots, expressed as a fraction of the average axis height
	'''
	
	plt.figure(1)
	plt.tight_layout();
	
	motif_grid.update(wspace=1.8,hspace=0.53)
	#plt.subplot(rama_gs[  0,0  ]).annotate('(a)', xy=(.0, 1), xycoords='figure fraction', horizontalalignment='left', verticalalignment='top', fontsize=textsize*1.2)
	#plt.subplot(number_gs[0,0:3]).annotate('(b)', xy=(.0, 0.5), xycoords='figure fraction', horizontalalignment='left', verticalalignment='top', fontsize=textsize*1.2)
	print "WRITING TO:",figname2
	plt.savefig(figname2,dpi=170,bbox_inches="tight")
	if show_graphs:
		os.system(pdf_viewer+" "+figname2)
	
	exit()
	
