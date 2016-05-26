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

rhomin = plotmap.ramachandran_number_collapse(-180,-180)
rhomax = plotmap.ramachandran_number_collapse( 180, 180)
rhomax_minus_rhomin = rhomax - rhomin

sequences = []

#alpha_phi = -63
#alpha_psi = -42

pdb_database_dir = "/home/ranjan/Desktop/ensembles/pdbstyle-2.06"
files = glob.glob(pdb_database_dir+"/*/*.ent")

counter = 0
#files = glob.glob("/home/ranjan/research/projects/ramachandran/pdb/dssp_outputs/*.ent")#[:5100]
if 1:
	window = 2
	index_of_interest = 0
	
	sequence_to_data = {}
	ss_vs_rho = {}
	for fn in files:#[:100]:
		counter += 1
		if counter%100==0:
			print counter,"   of",len(files), "\t("+str(round(100.0*float(counter)/len(files),2))+"%)"
		
		f = open(fn,"r")
		block = f.read()
		f.close()
		
		#deme = []
		model_to_chain_to_resno_atom_to_vals = plotmap.read_pdb(block)
		for model in sorted(model_to_chain_to_resno_atom_to_vals.keys()):
			for chain in sorted(model_to_chain_to_resno_atom_to_vals[model].keys()):
				resno_to_data = {}
				for resno in sorted(model_to_chain_to_resno_atom_to_vals[model][chain].keys()):
					vals = model_to_chain_to_resno_atom_to_vals[model][chain][resno]
					if 'r' in vals:
						aa_three_letter = vals['resname'].upper()
						if aa_three_letter in aa_three_to_one:
							aa_single_letter = aa_three_to_one[aa_three_letter]
							'resname'
							'phi'
							'psi'
							'ca'
							'c'
							'n'
							'r'
							#print vals.keys()
							resnumber = int(resno)
							#ss_type   = line[16]
							#chirality = line[22]
							resno_to_data[resnumber] = {}
							resno_to_data[resnumber]["aminoacid"] = aa_single_letter
							resno_to_data[resnumber]["phi"] = vals['phi']
							resno_to_data[resnumber]["psi"] = vals['psi']
							resno_to_data[resnumber]["rho"] = vals['r']
							#resno_to_data[resnumber]["ss_type"] = ss_type
							#deme.append(vals['r'])
				
				resnumbers = sorted(resno_to_data.keys())
				for rindex in range(len(resnumbers)):
					if rindex+window-1 < len(resnumbers):
						seq = ""
						if resnumbers[rindex]+window-1 == resnumbers[rindex+window-1]:
							for i in range(window):
								resno = resnumbers[rindex+i]
								aa = resno_to_data[resno]["aminoacid"]
								seq += aa
							#aa1 = resno_to_data[rh]["aminoacid"]
						if len(seq):
							ri = resnumbers[rindex+index_of_interest]
							if not seq in sequence_to_data:
								sequence_to_data[seq] = []
							sequence_to_data[seq].append(copy.deepcopy(resno_to_data[ri]))
		
		"""
		print "#deme"
		a,b = np.histogram(deme,bins=100)
		for i in range(len(a)):
			print float(b[i]+b[i+1])/2.0, a[i]
		exit()
		"""
	
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
		'T' : 0.452077674478}
	
	"""
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
	"""
	
	for ss in ["E","H","H'","P"]:
		sscode = ss
		if ss[-1] == "'":
			ylines[1-np.average(ss_to_r[ss[:-1]])] = sscode.lower()
		else:
			ylines[np.average(ss_to_r[ss])] = sscode.lower()
	
	#ylines[0.599251876766] = sscode_to_latex_code["P"]
	
	print "#window =",window
	print "#index_of_interest =",index_of_interest
	
	x = []
	y = []
	z = []
	bins = 100
	index = 0.0
	
	offset = 0
	
	r_val_to_seq = {}
	for seq in sorted(sequence_to_data.keys()):
		rhos = []
		for i in range(len(sequence_to_data[seq])):
			phi = sequence_to_data[seq][i]["phi"]
			psi = sequence_to_data[seq][i]["psi"]
			rho = sequence_to_data[seq][i]["rho"]
			rhos.append(rho)
		
		rval = np.average(rhos)#)/len(sequence_to_data[seq])
		if not rval in r_val_to_seq:
			r_val_to_seq[rval]=[]
		r_val_to_seq[rval].append(seq)
	
	sequences = []
	if 1:
		for rval in sorted(r_val_to_seq.keys()):
			for s in sorted(r_val_to_seq[rval]):
				sequences.append(s)
	
	sequences = sorted(sequences)
	
	#sequences = "RKDEQNHSTYCMWFILVAPG"
	#sequences.reverse()
	
	
	masterpattern = ""
	for i in range(0,index_of_interest):
		masterpattern+="X"
	masterpattern+="_"
	for i in range(index_of_interest+1,window):
		masterpattern+="X"
	
	patterns = []
	for seq in sequences:
		patterns.append(seq[0:index_of_interest]+"_"+seq[index_of_interest+1:window])
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
	
	for pattern in patterns:
		basedir_ramachandran_plot = basedir+"/"+pattern
		if not os.path.isdir(basedir_ramachandran_plot):
			os.makedirs(basedir_ramachandran_plot)
		
		get_pattern = re.compile(pattern.replace("_","."))
		get_residue = re.compile(pattern.replace("_","(.)"))
		
		x = []
		y = []
		z = []
		current_sequences = []
		
		index = 0
		for seq in sequences:#sorted(sequence_to_data.keys()):
			if get_pattern.match(seq):
				q = get_residue.match(seq)
				central_aa = q.group(1)
				#print [central_aa,seq,pattern]
				#raw_input()
				current_sequences.append(central_aa)
				#sequences.append(seq)
				
				Y_vals = []
				RAWPHI = []
				RAWPSI = []
				for sinfo in sequence_to_data[seq]:
					aminoacid = sinfo["aminoacid"]
					#ss_type = sinfo["ss_type"]
					phi = sinfo["phi"]
					psi = sinfo["psi"]
					rho = sinfo["rho"]
					#rho = plotmap.ramachandran_number_collapse(phi,psi)
					Y_vals.append(rho)
					
					RAWPHI.append(phi)
					RAWPSI.append(psi)
				if 1: # Draw individual ramachandran plots
					stepsize = 10
					PHI = []
					PSI = []
					WEIGHT = []
					for refphi in range(-180,181,stepsize):
						refphi_plus_step = refphi + stepsize
						for refpsi in range(-180,181,stepsize):
							refpsi_plus_step = refpsi + stepsize
							totalcounts = 0.0
							counts = 0.0
							for rawphi,rawpsi in zip(RAWPHI,RAWPSI):
								if refphi <= rawphi and rawphi <= refphi_plus_step:
									if refpsi <= rawpsi and rawpsi <= refpsi_plus_step:
										counts += 1.0
								totalcounts += 1.0
							#
							PHI.append(refphi)
							PSI.append(refpsi)
							WEIGHT.append(counts/totalcounts)
					
					maxWEIGHT = max(WEIGHT)
					for i in range(len(WEIGHT)):
						WEIGHT[i] = WEIGHT[i]/maxWEIGHT
					
					PHI=np.array(PHI)
					PSI=np.array(PSI)
					WEIGHT=np.array(WEIGHT)
					
					fn = basedir_ramachandran_plot+"/"+seq+".eps"
					"""
					fn_bar = fn+"_colorbar.eps"
					print "#WRITING TO:",fn_bar
					make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=cmap, fn=fn_bar, yscaling=0.5,xtitle="Key",ytitle="R")
					#
					"""
					#print "#WRITING TO:",fn
					plotmap.make2Dfigure(PHI,PSI,WEIGHT,fn,yscaling=1,xtitle="phi",ytitle="psi",xlabels=range(-180,181,90),ylabels=range(-180,181,90),showeps=0)
					#raw_input()
				
				a,b = np.histogram(Y_vals,bins=bins,range=[0.29,0.71])
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
				masterxlabels.append(pattern.replace("_",""))
			else:
				masterxlabels.append("")
		
		
		
		fn = basedir+"/"+pattern+".eps"
		
		sortedSetX = sorted(set(x))
		sortedSetXspacingBy2 = float(sortedSetX[1]-sortedSetX[0])/2.0
		xlines = {}
		xlines[sortedSetX[0]-sortedSetXspacingBy2] = ""
		for xnow in sortedSetX:
			xlines[xnow+sortedSetXspacingBy2] = ""
		
		master_horizontal_lines = ylines
		
		xlines = []
		plotmap.make2Dfigure(x,y,z,fn,xlabels=current_sequences,horizontallines=ylines,verticallines=xlines,ytitle="R",xtitle=pattern.replace("_","-"),showeps=0)#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
	#
	if window > 1:
		fn = basedir+"/all.eps"
		plotmap.make2Dfigure(masterx,mastery,masterz,fn,xlabels=masterxlabels,
			horizontallines=master_horizontal_lines,verticallines=master_vertical_lines,
			xtitle=masterpattern,ytitle="R",showeps=0)#,title=pattern.replace("_","-"))#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
