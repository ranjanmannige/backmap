#!/usr/bin/env python
import urllib2
import os,sys,glob,copy,re
import matplotlib.pyplot as plt 
from matplotlib import colors
import numpy as np
import numpy,math
import plotmap

show_graphs = 1

default_fontsize = 22

figure_base="./sequence_studies/"
figure_base = figure_base.rstrip("/")
if not os.path.isdir(figure_base):
	os.makedirs(figure_base)

secondary_structure_code_to_name = {
	"G":"3-turn helix ($3_{10}$-helix)",
	"H":"4-turn helix ($\\alpha$-helix)",
	"I":"5-turn helix ($\\pi$-helix)",
	"T":"Turn (3, 4 or 5 turn)",
	"E":"$\\beta$-sheet strand",
	"B":"Isolated beta-bridge",
	"S":"Bend",
	"C":"Coil",
	"-":"Undesignated",
	"P":"PP-II"
}

cm = plt.cm.bone_r
secondary_structure_code_to_cmap = {
	"G":cm,   # "3-turn helix ($3_{10}$-helix)",
	"H":cm,      #"4-turn helix ($\\alpha$-helix)",
	"I":cm,   #"5-turn helix ($\\pi$-helix)",
	"T":cm,   #"Hydrogen bonded turn (3, 4 or 5 turn)",
	"E":cm,    #"Extended $\\beta$-sheet strand",
	"B":cm,   #"Residue in isolated beta-bridge",
	"S":cm,   #"Bend",
	"C":cm,   #"Coil",
	"-":cm,   #"Undesignated"
	"P":cm    #"PP-II"
}


aa_three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

#rhomin = plotmap.ramachandran_number_collapse(-180,-180)
#rhomax = plotmap.ramachandran_number_collapse( 180, 180)
#rhomax_minus_rhomin = rhomax - rhomin

sequences = []

#alpha_phi = -63
#alpha_psi = -42


pdb_database_dir = "/home/ranjan/Desktop/ensembles/pdbstyle-2.06"
files = glob.glob(pdb_database_dir+"/*/*.ent")

# ----------------
# NEW COLOR SCHEME: color by backbone twist (expected range: R=[0,1])
# ----------------
# This lets you later on get the cmap by name 'TwoColor': cmap = plt.get_cmap('TwoColor')
# POSITION: 0              0.25             0.5           0.75            1
#    COLOR: | white - black | yellow - white | white - red | blue - white |
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
bc = 0.8
cdict = {
#                         
	'red':   ((0.00,  1, 1), (0.0000001,  1, bc), (1.0, 0, 0)), 
	'green': ((0.00,  1, 1), (0.0000001,  1, bc), (1.0, 0, 0)),
	'blue':  ((0.00,  1, 1), (0.0000001,  1, bc), (1.0, 0, 0)) 
}
cmap = LinearSegmentedColormap('Ramachandran', cdict)
plt.register_cmap(cmap=cmap)


def write_ramachandran_plot(RAWPHI,RAWPSI,fn="deme.eps",cmap=cmap,write_labels=1):
	stepsize = 5
	divide_by = 5.0
	phi_psi_to_counts = {}
	i = 0
	for rawphi,rawpsi in zip(RAWPHI,RAWPSI):
		rawphi = round(float(rawphi)/divide_by,0)*divide_by
		rawpsi = round(float(rawpsi)/divide_by,0)*divide_by
		t = (rawphi,rawpsi)
		if not t in phi_psi_to_counts:
			phi_psi_to_counts[t] = 0.0
		phi_psi_to_counts[t] += 1.0
	
	PHI = []
	PSI = []
	WEIGHT = []
	for phi,psi in phi_psi_to_counts.keys():
		PHI.append(phi)
		PSI.append(psi)
		WEIGHT.append(phi_psi_to_counts[(phi,psi)])
	
	maxWEIGHT = max(WEIGHT)
	for i in range(len(WEIGHT)):
		WEIGHT[i] = WEIGHT[i]/maxWEIGHT
	
	PHI=np.array(PHI)
	PSI=np.array(PSI)
	WEIGHT=np.array(WEIGHT)
	
	if write_labels:
		plotmap.make2Dfigure(PHI,PSI,WEIGHT,fn,xscaling=1, cmap=cmap,xtitle="psi",ytitle="psi",xticks=range(-180,181,90),xlabels=range(-180,181,90),yticks=range(-180,181,90),ylabels=range(-180,181,90),xlim=[-180,180],ylim=[-180,180],showeps=0)
	else:
		plotmap.make2Dfigure(PHI,PSI,WEIGHT,fn,xscaling=1, cmap=cmap,xticks=[],xlabels=[],yticks=[],ylabels=[],xlim=[-180,180],ylim=[-180,180],showeps=0)
	
	return 1

counter = 0
#files = glob.glob("/home/ranjan/research/projects/ramachandran/pdb/dssp_outputs/*.ent")#[:5100]
if 1:
	window = 2
	index_of_interest = 0
	write_ramachandran_plots = 0 # this will be set to "off" when window = 0 (below)
	window_and_index_of_interests = [(1,0),(2,0),(2,1)]
	#window_and_index_of_interests = [(1,0),(2,1)]
	
	cmap_for_rcode = plt.get_cmap("gray_r")
	
	signed = 0
	
	total_windows = []
	for window, index_of_interest in window_and_index_of_interests:
		total_windows.append(window)
	total_windows = list(set(total_windows))
	
	rrange = [0.29,0.71]
	stepsize = 0.004
	binrange = np.arange(rrange[0]-stepsize/2,rrange[-1]+stepsize,stepsize)
	if signed == 1:
		rrange = [-71,71]
		binrange = np.arange(rrange[0]-stepsize,rrange[-1]+stepsize*2.0,stepsize*2.0)
	
	# ---------->
	# code for creating comma separated value files
	get_scop_classification = re.compile("REMARK\s+99\s+ASTRAL\s+SCOPe\-sccs\:\s+([^\s]+)",re.M)
	tag_preceeding_SCOP_classification = "REMARK  99 ASTRAL SCOPe-sccs:"
	len_tag_preceeding_SCOP_classification = len(tag_preceeding_SCOP_classification)
	csfcount = 0
	#
	step_sizesForCSV = [0.005,0.01,0.05,0.1]
	binrangesForCSV = []
	fCSV = []
	for step_size in step_sizesForCSV:
		rrange1 = [0,1]
		
		binrange1 = np.arange(rrange1[0]-step_size/2,rrange1[-1]+step_size,step_size)
		if signed == 1:
			rrange1 = [-1,1]
			binrange1 = np.arange(rrange1[0]-step_size,rrange1[-1]+step_size*2,step_size*2)
			
		csv_fn = "rcode_to_scop_class_rstep"+str(binrange1[1]-binrange1[0])+".csv"
		print "WRITING TO:",csv_fn
		
		binrangesForCSV.append(binrange1)
		fCSV.append(open(csv_fn,"w"))
	#
	# <----------
	
	code_for_all_matches = "-"
	sequence_to_r   = {}
	sequence_to_phi = {}
	sequence_to_psi = {}

	ss_vs_rho = {}
	for fn in files:#[:100]:
		counter += 1
		if counter%50==0:
			print counter,"   of",len(files), "\t("+str(round(100.0*float(counter)/len(files),2))+"%)"
		
		model_to_chain_to_resno_atom_to_vals = plotmap.read_pdb(fn,signed=signed)
		
		# ---------->
		# code for creating comma separated value files
		f = open(fn,"r")
		block = f.read()
		f.close()
		SCOPclassification = ""
		SCOPclassification = block[block.index(tag_preceeding_SCOP_classification)+len_tag_preceeding_SCOP_classification:]
		SCOPclassification = SCOPclassification[:SCOPclassification.index("\n")].lstrip().rstrip()
		if 0: # alternative way to get the SCOP 
			q = get_scop_classification.search(block)
			if q:
				SCOPclassification = q.group(1)
				print "Found SCOP classification:",q.group(1)
			else:
				print "Did not find SCOP classification for `"+fn+"'"
				exit()
		#
		# <----------
		
		for model in sorted(model_to_chain_to_resno_atom_to_vals.keys()):
			for chain in sorted(model_to_chain_to_resno_atom_to_vals[model].keys()):
				resno_to_data = {}
				resnos = sorted(model_to_chain_to_resno_atom_to_vals[model][chain].keys())
				resnos_len = len(resnos)
				resno_counter = 0
				rnumbers = []
				for resno in resnos: #sorted(model_to_chain_to_resno_atom_to_vals[model][chain].keys()):
					if 'r' in model_to_chain_to_resno_atom_to_vals[model][chain][resno]:
						aa            = model_to_chain_to_resno_atom_to_vals[model][chain][resno]["aa"]
						rnumbers.append(model_to_chain_to_resno_atom_to_vals[model][chain][resno]["r"])
						if len(aa) == 1:
							# options: ['resname','phi','psi','ca','c','n','r']
							for window in total_windows:
								resno_plus_window_minus_one = resno+window-1
								resno_counter_plus_window   = resno_counter+window
								resno_counter_plus_window_minus_one = resno_counter_plus_window-1
								# 0 1 2 3 4 5
								#         |
								if resno_counter_plus_window_minus_one < resnos_len:
									if resno_plus_window_minus_one == resnos[resno_counter_plus_window_minus_one]:
										#print resno,"\t",resnos[resno_counter_plus_window_minus_one]
										seq  = aa
										phis =  [model_to_chain_to_resno_atom_to_vals[model][chain][resno]["phi"]]
										psis =  [model_to_chain_to_resno_atom_to_vals[model][chain][resno]["psi"]]
										ramachandran_numbers = [model_to_chain_to_resno_atom_to_vals[model][chain][resno]["r"]]
										
										for next_resno_index in range(resno_counter+1,resno_counter_plus_window):
											current_resno = resnos[next_resno_index]
											if "r" in model_to_chain_to_resno_atom_to_vals[model][chain][current_resno]:
												seq += model_to_chain_to_resno_atom_to_vals[model][chain][current_resno]["aa"]
												phis.append(model_to_chain_to_resno_atom_to_vals[model][chain][current_resno]["phi"])
												psis.append(model_to_chain_to_resno_atom_to_vals[model][chain][current_resno]["psi"])
												ramachandran_numbers.append(model_to_chain_to_resno_atom_to_vals[model][chain][current_resno]["r"])
										
										if len(seq) == window:
											if not seq in sequence_to_r:
												sequence_to_r[seq]   = []
												sequence_to_phi[seq] = []
												sequence_to_psi[seq] = []
											#
											sequence_to_phi[seq].append(phis)
											sequence_to_psi[seq].append(psis)
											sequence_to_r[seq].append(ramachandran_numbers)
					resno_counter += 1
				#
				
				# ---------->
				# code for creating comma separated value files
				if len(rnumbers):
					for step_size in step_sizesForCSV:
						csvcount = step_sizesForCSV.index(step_size)	
						binrangeForHistogram = binrangesForCSV[csvcount]
						
						a,b = np.histogram(rnumbers,bins=binrangeForHistogram)#,range=[rrange[0],rrange[1]])#range=[0.29,0.71])
						#print b
						
						totala = np.sum(a)
						maxa = np.max(a)
						
						if csfcount == 0:
							fCSV[csvcount].write(",")
							for i in range(len(a)):
								currenty = float(b[i]+b[i+1])/2.0
								#currentz = float(a[i])
								#if maxa:
								#	currentz = currentz/maxa
								fCSV[csvcount].write(str(currenty)+",")
							fCSV[csvcount].write("class label\n")
						
						fCSV[csvcount].write(str(csfcount)+",")
						for i in range(len(a)):
							#currenty = float(b[i]+b[i+1])/2.0
							currentz = 0
							if totala:
								currentz = float(a[i])/totala
							#if maxa:
							#	currentz = currentz/maxa
							fCSV[csvcount].write(str(currentz)+",")
						fCSV[csvcount].write(SCOPclassification+"\n")
						#csf.flush()
					csfcount += 1
				# <----------
	
	ylines = {}
	sscode_to_latex_code = {
		"H" :"h", #"$\\alpha$",
		"H'":"h'", #"$\\alpha'$",
		"G" :"g",
		"E" :"e",
		"P" :"p"
		}
	
	ss_to_r = {' ' : 0.558636024139,
		'B' : 0.537760481375,
		'E' : 0.515806491698,
		'G' : 0.386323496778,
		'I' : 0.360462931421,
		'H' : 0.355013406346,
		'P' : 0.599251876766,
		'S' : 0.462813874371,
		'T' : 0.452077674478
		}
	
	'''
	"G":"3-turn helix ($3_{10}$-helix)",
	"H":"4-turn helix ($\\alpha$-helix)",
	"I":"5-turn helix ($\\pi$-helix)",
	"T":"Turn (3, 4 or 5 turn)",
	"E":"$\\beta$-sheet strand",
	"B":"Isolated beta-bridge",
	"S":"Bend",
	"C":"Coil",
	"-":"Undesignated",
	"P":"PP-II"
	'''
	
	for ss in ["E","H","H'","P"]:
		sscode = ss
		if ss[-1] == "'":
			location = 1-np.average(ss_to_r[ss[:-1]])
			if signed:
				location = location * -1.0
			ylines[location] = sscode.lower()
			
		else:
			ylines[np.average(ss_to_r[ss])] = sscode.lower()
	
	#ylines[0.599251876766] = sscode_to_latex_code["P"]
	
	
	for seq in sequence_to_r.keys():
		sequence_to_r[seq] = np.array(sequence_to_r[seq])
		sequence_to_phi[seq] = np.array(sequence_to_phi[seq])
		sequence_to_psi[seq] = np.array(sequence_to_psi[seq])
	
	for window,index_of_interest in window_and_index_of_interests:
		print "#window =",window
		print "#index_of_interest =",index_of_interest
		
		x = []
		y = []
		z = []
		bins = 100
		index = 0.0
		
		offset = 0
		
		rangewindow  = range(window)
		sequences = []
		
		for seq in sorted(sequence_to_r.keys()):
			if len(seq) == window:
				sequences.append(seq)
		#
		sequences = sorted(sequences)
		#
		masterpattern = ""
		for i in range(0,index_of_interest):
			masterpattern+="X"
		masterpattern+="_"
		for i in range(index_of_interest+1,window):
			masterpattern+="X"
		
		
		patterns = []
		for seq in sequences:
			patterns.append(seq[0:index_of_interest]+"_"+seq[index_of_interest+1:window])
		
		universal_pattern = ""
		for i in range(index_of_interest):
			universal_pattern+="."
		universal_pattern+="_"
		for i in range(index_of_interest+1,window):
			universal_pattern+="."
		
		
		patterns = [universal_pattern]+patterns
		patterns = sorted(set(patterns))
		
		masterx = []
		mastery = []
		masterz = []
		master_horizontal_lines = []
		master_vertical_lines   = []
		masterxlabels = []
		final_index = 0
		
		basedir = "motif_studies/window"+str(window)+"_index"+str(index_of_interest)
		if not os.path.isdir(basedir):
				os.makedirs(basedir)
		
		colorbarXscaling = 0.08
		cbfn = basedir+"/colorbar.eps"
		color_bar_range = np.arange(0,1.005,0.01)
		plotmap.make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=cmap_for_rcode, fn=cbfn, title="P",showeps=0,xscaling=colorbarXscaling)
		
		for pattern in patterns:
			
			stripped_pattern = pattern.replace(".","").replace("_","")
			len_stripped_pattern = len(stripped_pattern)
			
			#print [pattern,pattern.replace("_","."),pattern.replace("_","(.)"),stripped_pattern]
			#exit()
			
			get_pattern = re.compile(pattern.replace("_","."))
			get_residue = re.compile(pattern.replace("_","(.)"))
			
			x = []
			y = []
			z = []
			current_sequences = []
			
			index = 0
			
			relevant_sequences = []
			residue_to_sequence = {"-":[]}
			for seq in sequences:#sorted(sequence_to_data.keys()):
				if get_pattern.match(seq):# or len_stripped_pattern == 0:
					#relevant_sequences.append(seq)
					q = get_residue.match(seq)
					if q:
						central_aa = "-"
						if seq != "-":
							central_aa = q.group(1)
						
						if not central_aa in residue_to_sequence:
							residue_to_sequence[central_aa] = []
						#
						residue_to_sequence[central_aa].append(seq)
						residue_to_sequence["-"].append(seq)
			#
			for central_aa in sorted(residue_to_sequence.keys()):
				Y_vals = []
				RAWPHI = []
				RAWPSI = []
				for seq in residue_to_sequence[central_aa]:
					Y_vals += list(sequence_to_r[seq][:,index_of_interest].T)
					RAWPHI += list(sequence_to_phi[seq][:,index_of_interest].T)
					RAWPSI += list(sequence_to_psi[seq][:,index_of_interest].T)
				
				# MAYBE HERE WE CAN DO SOMETHING LIKE ADD A BLANK LINE IF THERE ARE NO HITS
				if len(RAWPHI):
					seq
					current_sequences.append(central_aa.replace("-","*"))
					print "----------"
					print "pattern:",pattern,"\tcentral_aa:",central_aa,"\tlen(residue_to_sequence[central_aa])",len(residue_to_sequence[central_aa]),"\tseq:",seq
					#print ["\tpattern:",pattern,"\tseq:",seq]
					write_ramachandran_plots = 1
					if window > 1:
						write_ramachandran_plots = 0
					
					if write_ramachandran_plots: # Draw individual ramachandran plots
						print "deme"
						# ---------------------------------
						# MAKING THE DIRECTORIES, IF ABSENT
						basedir_ramachandran_plot1          = basedir+"/ramachandran_plots/color1"
						basedir_ramachandran_plot1labelfree = basedir+"/ramachandran_plots/color1/label_free"
						basedir_ramachandran_plot2          = basedir+"/ramachandran_plots/color2"
						basedir_ramachandran_plot2labelfree = basedir+"/ramachandran_plots/color2/label_free"
						
						if not os.path.isdir(basedir_ramachandran_plot1):
							os.makedirs(basedir_ramachandran_plot1)
						if not os.path.isdir(basedir_ramachandran_plot2):
							os.makedirs(basedir_ramachandran_plot2)
						if not os.path.isdir(basedir_ramachandran_plot1labelfree):
							os.makedirs(basedir_ramachandran_plot1labelfree)
						if not os.path.isdir(basedir_ramachandran_plot2labelfree):
							os.makedirs(basedir_ramachandran_plot2labelfree)
						
						# ---------------------------------
						# NEW COLOR SCHEME
						cmap = plt.get_cmap("gray_r")
						fn = basedir_ramachandran_plot1+"/ramachandran_"+central_aa+".eps"
						write_ramachandran_plot(RAWPHI,RAWPSI,fn=fn,cmap=cmap,write_labels=1)
						
						if index == 0:
							# draw the colorbar, but only once
							cbfn = basedir_ramachandran_plot1+"/colorbar.eps"
							color_bar_range = np.arange(0,1.005,0.01)
							plotmap.make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=cmap, fn=cbfn, title="P",showeps=0,xscaling=colorbarXscaling)
						#
						# If you want to write ONLY the data, no labels (sometimes useful), uncomment the following and comment the previous.
						fn = basedir_ramachandran_plot1labelfree+"/ramachandran_"+central_aa+".eps"
						write_ramachandran_plot(RAWPHI,RAWPSI,fn=fn,cmap=cmap,write_labels=0)
						if index == 0:
							# draw the colorbar, but only once
							cbfn = basedir_ramachandran_plot1labelfree+"/colorbar.eps"
							color_bar_range = np.arange(0,1.005,0.01)
							plotmap.make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=cmap, fn=cbfn, title="P",showeps=0,xscaling=colorbarXscaling)
						#
						# ---------------------------------
						# NEW COLOR SCHEME
						cmap = plt.get_cmap('Ramachandran')
						fn = basedir_ramachandran_plot2+"/ramachandran_"+central_aa+".eps"
						write_ramachandran_plot(RAWPHI,RAWPSI,fn=fn,cmap=cmap,write_labels=1)
						if index == 0:
							# draw the colorbar, but only once
							cbfn = basedir_ramachandran_plot2+"/colorbar.eps"
							color_bar_range = np.arange(0,1.005,0.01)
							plotmap.make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=cmap, fn=cbfn, title="P",showeps=0,xscaling=colorbarXscaling)
						#
						# If you want to write ONLY the data, no labels (sometimes useful), uncomment the following and comment the previous.
						fn = basedir_ramachandran_plot2labelfree+"/ramachandran_"+central_aa+".eps"
						write_ramachandran_plot(RAWPHI,RAWPSI,fn=fn,cmap=cmap,write_labels=0)
						if index == 0:
							# draw the colorbar, but only once
							cbfn = basedir_ramachandran_plot2labelfree+"/colorbar.eps"
							color_bar_range = np.arange(0,1.005,0.01)
							plotmap.make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=cmap, fn=cbfn, title="P",showeps=0,xscaling=colorbarXscaling)
						#raw_input()
					
					a,b = np.histogram(Y_vals,bins=binrange)#,range=[rrange[0],rrange[1]])#range=[0.29,0.71])
					#print b
					maxa = np.max(a)
					for i in range(len(a)):
						
						currenty = float(b[i]+b[i+1])/2.0
						currentz = 0.0
						if maxa:
							currentz = float(a[i])/maxa
						
						x.append(index-offset)
						y.append(currenty)
						z.append(currentz)
						masterx.append(final_index-offset)
						mastery.append(currenty)
						masterz.append(currentz)
					index+=1
					final_index+=1
			
			sortedX = sorted(set(x))
			master_vertical_lines.append(final_index-offset)
			for i in range(len(sortedX)):
				if i == len(sortedX)/2:
					masterxlabels.append(pattern.replace("_","").replace(".","*"))
				else:
					masterxlabels.append("")
			#
			fn = basedir+"/"+pattern+".eps"
			
			sortedSetX = sorted(set(x))
			sortedSetXspacingBy2 = float(sortedSetX[1]-sortedSetX[0])/2.0
			xlines = {}
			xlines[sortedSetX[0]-sortedSetXspacingBy2] = ""
			for xnow in sortedSetX:
				xlines[xnow+sortedSetXspacingBy2] = ""
			
			master_horizontal_lines = ylines
			
			xlines = []
			plotmap.make2Dfigure(x,y,z,fn,xlabels=current_sequences,horizontallines=ylines,verticallines=xlines,ytitle="R",xtitle=pattern.replace("_","-"),showeps=0, cmap=cmap_for_rcode, ylim = rrange)
		#
		if window > 1:
			fn = basedir+"/all.eps"
			plotmap.make2Dfigure(masterx,mastery,masterz,fn,xlabels=masterxlabels,
				horizontallines=master_horizontal_lines,verticallines=master_vertical_lines,
				xtitle=masterpattern,ytitle="R",showeps=0, cmap=cmap_for_rcode, ylim = rrange)#,title=pattern.replace("_","-"))#,zlim=[0.0,0.3])# ylim = [0.3,0.71])
	