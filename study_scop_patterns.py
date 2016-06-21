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


counter = 0
#files = glob.glob("/home/ranjan/research/projects/ramachandran/pdb/dssp_outputs/*.ent")#[:5100]
if 1:
	window = 1
	index_of_interest = 0
	
	cmap = plt.get_cmap('Ramachandran')
	cmap_for_rcode = plt.get_cmap("gray_r")
	#cmap = plt.get_cmap('Ramachandran')
	#cmap_for_rcode = plt.get_cmap('Ramachandran')
	signed = 0
	
	rrange = [0,1]
	step_size = 0.05
	binrange = np.arange(rrange[0]-step_size/2,rrange[-1]+step_size,step_size)
	if signed == 1:
		rrange = [-1,1]
		binrange = np.arange(rrange[0]-step_size,rrange[-1]+step_size*2,step_size*2)
	
	code_for_all_matches = "-"
	
	get_scop_classification = re.compile("REMARK\s+99\s+ASTRAL\s+SCOPe\-sccs\:\s+([^\s]+)",re.M)
	tag_preceeding_SCOP_classification = "REMARK  99 ASTRAL SCOPe-sccs:"
	len_tag_preceeding_SCOP_classification = len(tag_preceeding_SCOP_classification)
	
	'''
	sequence_to_data = {code_for_all_matches:[]}
	'''
	ss_vs_rho = {}
	
	csfcount = 0
	csv_fn = "rcode_to_scop_class_rstep"+str(binrange[1]-binrange[0])+".csv"
	print "WRITING TO:",csv_fn
	csf = open(csv_fn,"w")
	
	for fn in files:#[:100]:
		counter += 1
		if counter%50==0:
			print counter,"   of",len(files), "\t("+str(round(100.0*float(counter)/len(files),2))+"%)"
		
		model_to_chain_to_resno_atom_to_vals = plotmap.read_pdb(fn,signed=signed)
		
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
		
		for model in sorted(model_to_chain_to_resno_atom_to_vals.keys()):
			for chain in sorted(model_to_chain_to_resno_atom_to_vals[model].keys()):
				resno_to_data = {}
				rnumbers = []
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
							rnumbers.append(vals['r'])
							#resno_to_data[resnumber]["ss_type"] = ss_type
							#deme.append(vals['r'])
				
				'''
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
							sequence_to_data[code_for_all_matches].append(copy.deepcopy(resno_to_data[ri]))
				'''
				
				a,b = np.histogram(rnumbers,bins=binrange)#,range=[rrange[0],rrange[1]])#range=[0.29,0.71])
				#print b
				
				totala = np.sum(a)
				#if totala < 30 and len(resnumbers):
				#	print totala,len(resnumbers),100.0*float(totala)/len(resnumbers)
				if totala:# > 30:
					if csfcount == 0:
						csf.write(",")
						for i in range(len(a)):
							currenty = float(b[i]+b[i+1])/2.0
							currentz = float(a[i])
							#if maxa:
							#	currentz = currentz/maxa
							csf.write(str(currenty)+",")
						csf.write("class label\n")
					
					
					
					csf.write(str(csfcount)+",")
					for i in range(len(a)):
						currenty = float(b[i]+b[i+1])/2.0
						currentz = float(a[i])/totala
						#if maxa:
						#	currentz = currentz/maxa
						csf.write(str(currentz)+",")
					csf.write(SCOPclassification+"\n")
					
					#csf.flush()
					csfcount += 1
				
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
	
	colorbarXscaling = 0.08
	cbfn = basedir+"/colorbar.eps"
	color_bar_range = np.arange(0,1.005,0.01)
	plotmap.make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=cmap_for_rcode, fn=cbfn, title="P",showeps=0,xscaling=colorbarXscaling)
	
	for pattern in patterns:
		basedir_ramachandran_plot = basedir+"/"+pattern
		if not os.path.isdir(basedir_ramachandran_plot):
			os.makedirs(basedir_ramachandran_plot)
		
		# colorbar
		colorbarXscaling = 0.08
		cbfn = basedir_ramachandran_plot+"/colorbar.eps"
		color_bar_range = np.arange(0,1.005,0.01)
		plotmap.make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=cmap, fn=cbfn, title="P",showeps=0,xscaling=colorbarXscaling)
		
		get_pattern = re.compile(pattern.replace("_","."))
		get_residue = re.compile(pattern.replace("_","(.)"))
		
		x = []
		y = []
		z = []
		current_sequences = []
		
		index = 0
		for seq in sequences:#sorted(sequence_to_data.keys()):
			if seq==code_for_all_matches or get_pattern.match(seq):
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
					stepsize = 5
					PHI = []
					PSI = []
					WEIGHT = []
					
					'''
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
							PHI.append(float(refphi+refphi_plus_step)/2.0)
							PSI.append(float(refpsi+refpsi_plus_step)/2.0)
							WEIGHT.append(counts/totalcounts)
					'''
					
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
					
					fn = basedir_ramachandran_plot+"/"+seq+".eps"
					'''
					fn_bar = fn+"_colorbar.eps"
					print "#WRITING TO:",fn_bar
					make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=cmap, fn=fn_bar, 
					=0.5,xtitle="Key",ytitle="R")
					#
					'''
					#print "#WRITING TO:",fn
					#plotmap.make2Dfigure(PHI,PSI,WEIGHT,fn,xscaling=1,showeps=0, cmap=cmap, xlabels=range(-180,181,90),ylabels=range(-180,181,90), xtitle="phi",ytitle="psi")
					# If you want to write ONLY the data, no labels (sometimes useful), uncomment the following and comment the previous.
					
					plotmap.make2Dfigure(PHI,PSI,WEIGHT,fn,xscaling=1,showeps=0, cmap=cmap, xlabels=[],ylabels=[])
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
		
		plotmap.make2Dfigure(x,y,z,fn,xlabels=current_sequences,horizontallines=ylines,verticallines=xlines,ytitle="R",xtitle=pattern.replace("_","-"),showeps=0, cmap=cmap_for_rcode)#,ylim=[-0.71,0.71])#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
	#
	if window > 1:
		fn = basedir+"/all.eps"
		plotmap.make2Dfigure(masterx,mastery,masterz,fn,xlabels=masterxlabels,
			horizontallines=master_horizontal_lines,verticallines=master_vertical_lines,
			xtitle=masterpattern,ytitle="R",showeps=0, cmap=cmap_for_rcode)#,title=pattern.replace("_","-"))#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
	"""