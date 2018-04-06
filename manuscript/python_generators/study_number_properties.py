#!/usr/bin/env python
# sections marked by "FIGURE_START, shouldbeone" indicate that a new figure is going to be made.
r'''
COMPANION SCRIPT #1, TESTED ONLY ON PYTHON 2.7, FOR: 
Mannige RV (2017) (article title TBD).
'''

length_dependent_csv_file = "local_imports/data_length.csv"
master_ss_file            = "local_imports/data_ss.csv"
fake_ss_file              = "local_imports/data_fake_ss.csv"

type_to_label = {'L':'$L$ (residues)','phi':'$\phi (^\circ)$',  'psi':'$\psi (^\circ)$', 'omega':'$\omega$', 'rg':'$r_{\mathrm{g}}$','re':'$r_{\mathrm{e}}$',
		 'd':'$d$', 'd2':'$|d|$' ,'theta':r'$\theta (^\circ)$',  'R':'$\mathcal{R}$'} #, (\phi+\psi+2\pi)/4\pi$'}

type_to_label = {'L':'$L$ (residues)','phi':'$\phi$',  'psi':'$\psi$', 'omega':'$\omega$', 'rg':'$r_{\mathrm{g}}$','re':'$r_{\mathrm{e}}$',
		 'd':'$d$', 'd2':'$|d|$' ,'theta':r'$\theta$',  'R':'$\mathcal{R}$'} #, (\phi+\psi+2\pi)/4\pi$'}

show_graphs = 1        # just uses a simple command line that 
pdf_viewer  = "evince" # This is the viewer of choice for the author (based-Linux).

# GLOBAL IMPORTS:
import os, sys, copy, random, string
import matplotlib.pyplot as plt        # For utilizing colormap interpolation 
from matplotlib.colors import LogNorm
import numpy as np
import Bio.PDB       # Biopython's PDB module
import pandas as pd
import matplotlib as mpl
from scipy import interpolate
import scipy.ndimage
import scipy.stats   # for KDE
# LOCAL IMPORTS
sys.path.insert(0, "./local_imports/") # for the local imports
import Geometry, PeptideBuilder, locallib
# -----------------------
# A PRETTY GRAPH IMPORT
seaborn_exists = 1
try:
	import seaborn as sns
except:
	seaborn_exists = 0
if seaborn_exists:
	sns.set_style('ticks') # *
	#sns.set_context("paper", font_scale=1.5)#, rc={"lines.linewidth": 4})
# -----------------------

# ===================================================================================
# SUPER IMPORTANT FOR ASTHETICS (call after every sns.style() function)
def set_grays():
	black_color = np.ones(3)*.5  # Darker is closer to 0; lighter is closer to 1
	text_color  = np.ones(3)*.35 # 

	# CHECK ALL OTHER PARAMETERS: print plt.rcParams
	#sns.set_style("whitegrid")
	plt.rc_context({'axes.labelcolor':text_color,'axes.linewidth': 1.0,'legend.edgecolor': black_color,
		            'lines.color':black_color,'lines.linewidth':1.0,'text.color':text_color,'axes.edgecolor':black_color,
		            'xtick.color':text_color,'ytick.color':text_color})
	
	#to retreive any of the parameter values, just do, e.g.:
	#some_value = plt.rcParams['lines.color']
	
	# Shows all other parameters to muck about with:
	#print plt.rcParams

	# For using latex commands in labels, etc.
	#mpl.rc('text',usetex=True)
	#mpl.rc('text.latex', preamble='\usepackage{color}')
#set_grays()
# ===================================================================================

import itertools
sns_palette   = sns.color_palette() #"colorblind") # another possibility: sns.color_palette()
palette       = itertools.cycle(sns_palette)
marker_shapes = itertools.cycle(['o','s', '^','p'])

# ===================================================================================
# SETTING UP A CUSTOM COLORMAP
# 
COLORSWITCH = 0.5; bc = [1,1,1] # background (white)
import colorsys
r  = [colorsys.rgb_to_hsv(1,0,0)[0],0.75,0.5]
y  = [colorsys.rgb_to_hsv(1,1,0)[0],0.5,0.75]
c3 = colorsys.hsv_to_rgb(*r) # the '*' converts [a,b,c] into a,b,c
c4 = colorsys.hsv_to_rgb(*y)
# Now the color map dictionary
cdict = {
	'red':   ((0.00,  c3[0], c3[0]), (COLORSWITCH,  bc[0], bc[0]), (1.0, c4[0], c4[0])), 
	'green': ((0.00,  c3[1], c3[1]), (COLORSWITCH,  bc[1], bc[1]), (1.0, c4[1], c4[1])),
	'blue':  ((0.00,  c3[2], c3[2]), (COLORSWITCH,  bc[2], bc[2]), (1.0, c4[2], c4[2])) 
}
ss_name_to_label      = {} #{'H':r'$\alpha$','E':r'$\beta$','G':r'$3_{10}$','P':r'$\mathrm{pp2}$'}
ss_name_to_label_long = {'H':r'$\alpha$ (H)','E':r'E ($\beta$)','G':r'G ($3_{10}$)','P':r'P ($\mathrm{pp2}$)'}

from matplotlib.colors import LinearSegmentedColormap # For making your own colormaps
cmap = LinearSegmentedColormap('chirality_r', cdict)
plt.register_cmap(cmap=cmap)


ss_codes = {}
ss_codes['dssp']  = {
	'H':0, # (0) Alpha - helix
	'G':1, # (1) 3.10 - helix        
	'I':2, # (2) Pi - helix
	'T':3, # (3) hydrogen bonded turn
	'S':4, # (4) bend
	' ':5, # (5) coil
	'B':6, # (6) residue in isolated beta-bridge
	'E':7, # (7) extended strand, participates in beta ladder
	'U':8  # (8) not found
}
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


def do_d(X):
	return X
	#return np.abs(X)
	
def do_theta(X):
	return X #-1.0*np.cos(np.radians(X)/2.0)
	#return X

def prepare_master_ss_file():
	## SS Database study
	# Mansiaux, Joseph, Gelly, de Brevern (2011) Assignment of PolyProline II Conformation and Analysis 
	# of Sequence-Structure Relationship. PLoS One 6(3): e18401.
	
	ss_database_filename = "local_imports/PPII/PPII_cullpdb_pc30_res2.5_R1.0_d090511_chains7172.db"
	#obtained from: http://www.dsimb.inserm.fr/~debrevern/DOWN/DB/PPII/PPII_cullpdb_pc30_res2.5_R1.0_d090511_chains7172.db
	#reported in: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0018401
	if not os.path.isfile(ss_database_filename):
		print "Local PPII database '%s' not available. Exiting." %(ss_database_filename)
		exit()
		#Download from:
		#http://www.dsimb.inserm.fr/~debrevern/DOWN/DB/PPII/PPII_cullpdb_pc30_res2.5_R1.0_d090511_chains7172.db
	
	# Ignore all lines starting with ">"
	ppII_colnames = ['aminoacids', 'proteinblocks', 'dssp', 'xtlsstr', 'segno', 'pross', 
	                 'relacc','phi','psi','omega','calphax','calphay','calphaz','temp', 'u1', 'u2','u3'] # the last three are 'u' for unknown
	
	ppII_df = pd.read_csv(ss_database_filename, sep='\s+', header=None, names=ppII_colnames, comment='>')
	ppII_df = ppII_df[(np.abs(ppII_df['phi']) != 999) 
	                & (np.abs(ppII_df['psi']) != 999)
	                & (np.abs(ppII_df['omega']) != 999)
	              & ( (np.abs(ppII_df['phi']) > 1) & (np.abs(ppII_df['phi']) > 1) ) # a weird space where phi=psi=0
	              ]
	
	#df = ppII_df[(ppII_df['phi'] > 0)]
	#sns.jointplot(df['phi'],df['psi'], kind="hex", color="k"); plt.show();
	#locallib.plot_ramachandran_histogram(df['phi'],df['psi'],levels=[0.9], fill=1,lines=0, show = 1,smooth=2)#!!!
	#exit()
	
	print "obtaining d,theta,R"
	ds = []
	thetas = []
	Rs = []
	
	total   = len(ppII_df['phi'])
	counter = 0 
	for phi,psi,omega in zip(ppII_df['phi'],ppII_df['psi'],ppII_df['omega']):
		d     = 999
		theta = 999
		R     = 999
		
		if phi and psi and omega:
			if -180.0 < phi and phi < 180.0 and -180.0 < psi and psi < 180.0:
				#h, theta, d = locallib.calculate_handedness_from_theory(phi, psi, omega)
				d,theta,rs  = locallib.calculate_d_theta_r(phi, psi, omega)
				R  = locallib.R(phi=phi, psi=psi)
		
		ds.append(d)
		thetas.append(theta)
		Rs.append(R)
	
		counter += 1
		if counter % 1000 == 0:
			sys.stdout.write("\r\tCompleted %d percent\t\t" %(int(100.*float(counter)/total)))
			sys.stdout.flush()
	print
	print "populating the new column entries"	
	ppII_df = ppII_df.assign(R=Rs)
	ppII_df = ppII_df.assign(d=ds)
	ppII_df = ppII_df.assign(theta=thetas)
	
	print "Saving to fn='%s' [open using 'pd.read_csv(fn)]'" %(master_ss_file)
	ppII_df.to_csv(master_ss_file)
	
	print "DONE"

def prepare_main_csv_file():
	# Generates structures for cis and trans backbones, calculates for each structure
	# the radius of gyration (rg) and the end-to-end distance (re), and saves them to a csv file
	# Relevant values saved to the file:
	columns = ['L','phi','psi','omega','d','theta','rg','re','R','R2']
	
	step_size = 1 # angle step size
	
	directory = os.path.split(length_dependent_csv_file)[0]	
	if not os.path.isdir(directory):
		os.makedirs(directory)
	
	print "Populating the file '%s' with the values %s"%(length_dependent_csv_file,str(columns))
	outputf = locallib.open_file(length_dependent_csv_file,'w')
	outputf.write(",".join(columns)+"\n")

	counter = 0
	phipsi_range = np.arange(-180,180.1,step_size)
	omegas       = [180,0]  
	lengths      = [16*2,16]   # length of peptide
	
	total_points = len(omegas)*len(lengths)*(len(phipsi_range)**2) # just for counting
	
	for L in lengths:
		for omega in [180,0]:
			print "\n--------\nCALCULATING VALUES FOR OMEGA = "+str(omega)+" PEPTIDES"
			for phi in phipsi_range:
				for psi in phipsi_range:
					counter += 1
					if counter % 20 == 0.0:
						str_phi = "{0: >4}".format(int(phi))
						str_psi = "{0: >4}".format(int(psi))
						sys.stdout.write("\r\tTotal completed %d percent\t" %(int(100.*float(counter)/total_points)))
						sys.stdout.flush()
					
					# ------------------
					# Building a peptide (Tien et al., PeptideBuilder, PeerJ, 2013)
					angles = []
					for n in range(int(L)):
						angles.append((phi,psi,omega))
					st1 = []
					st1 = locallib.build_structure(angles)
					# Only rg and re require structures to get values out
					rg = locallib.calculate_rg(st1)
					re = locallib.calculate_re(st1)
					# ------------------
				
					# d, theta, and R do not need structures to be calculated
					d, theta, rs = locallib.calculate_d_theta_r(phi, psi, omega)
					R  = locallib.R(phi=phi, psi=psi)
					R2 = (phi + psi + 360.)/(720.) # min = -360, max = 360
					str_vals = []
					for v in [L,phi,psi,omega,d,theta,rg,re,R,R2]:
						str_vals.append(str(v))
					outputf.write(",".join(str_vals)+"\n")
			sys.stdout.write("\n")
			#
	outputf.close()
	print "DONE"

def write_fake_secondary_structure_database(sstype = 'segno',omega=180.0,sample_points=100000, sigma=5):
	fake_ss_centroids = []
	angle_range = [-45-90,-45,45,45+90]
	for x in angle_range:
		for y in angle_range:
			if not (x,y) in fake_ss_centroids:
				fake_ss_centroids.append((x,y,omega))
	phi_start = -180
	psi_start = -180
	
	dictionary = {'phi':[],'psi':[],'omega':[],'d':[],'theta':[],'R':[],'R2':[],sstype:[]}
	
	ss_index = 0 # This will act as the secondary structure name
	
	for phi,psi,omega in fake_ss_centroids:
		ss_index+=1
		
		sys.stdout.write("\r\tSecondary structures generated: %d of %d\t\t" %(ss_index,len(fake_ss_centroids)))
		sys.stdout.flush()
		
		phis = []
		psis = []
		
		for trial in range(sample_points):
			rama_phi = random.gauss(phi, sigma=sigma)
			rama_psi = random.gauss(psi, sigma=sigma)
			
			phis.append(rama_phi)
			psis.append(rama_psi)			
			'''
			rama_phi2   = (rama_phi - phi_start) % 360. + phi_start   # 'phi_start'   declared at the beginning of this file
			rama_psi2   = (rama_psi - psi_start) % 360. + psi_start   # 'psi_start'   declared at the beginning of this file
			
			phis.append(rama_phi2)
			psis.append(rama_psi2)
			if rama_phi != rama_phi2 or rama_psi != rama_psi2:				
				phis.append(rama_phi)
				psis.append(rama_psi)
			'''
		
		#sns.jointplot(np.array(phis), np.array(psis)); plt.show();
		
		for phi,psi in zip(phis,psis):
			original_phi = (phi - phi_start  ) % 360. + phi_start   # 'phi_start'   declared at the beginning of this file
			original_psi = (psi - psi_start  ) % 360. + psi_start   # 'psi_start'   declared at the beginning of this file
			
			d,theta,rs = locallib.calculate_d_theta_r(original_phi,original_psi,omega) # Miyazawa values Mannige, PeerJ, 2017
			R = locallib.R(phi=original_phi,psi=original_psi) # Ramachandran number from Mannige, Kundu, Whitelam, PLoS ONE, 2016
			R2 = (original_phi+original_psi + 360.)/(720.) # min = -360, max = 360
			
			dictionary['phi'].append(phi)
			dictionary['psi'].append(psi)
			dictionary['omega'].append(omega)
			dictionary['d'].append(d)
			dictionary['theta'].append(theta)
			dictionary['R'].append(R)
			dictionary['R2'].append(R2)
			dictionary[sstype].append(ss_index)
	print 
	df = pd.DataFrame(dictionary)
	df.to_csv(fake_ss_file)
	print "Done"

#write_fake_secondary_structure_database()
#exit()

#
def calc_rg(d,L):
	#return 0.475398 * L * d
	#Delta = (0.477931 + 0.955227/2)/2
	Delta = 0.48
	return Delta * L * np.abs(d)
#
def get_rg_estimates(x,y=None,c='k',show=0,L=None):
	x2 = []
	y2 = [] 
	minval = np.min(x)
	maxval = np.max(x)
	steps = float(maxval-minval)/200.0 
	for currentx in np.arange(minval,maxval+steps/2.,steps):
		x2.append(currentx)
		y2.append(calc_rg(currentx,L))
	#
	if show and y:
		plt.plot(x2, y2, c=c)
		plt.show()
	#
	return x2,y2

#
def calc_re(d,L):
	#return 0.70833*d*(L**(2.0 - 0.90765))
	#Delta = (0.477931*2 + 0.955227)/2
	Delta = 0.48*2
	return Delta * L * np.abs(d)
#
def get_re_estimates(x,y=None,c='k',show=0,L=None):	
	x2 = []
	y2 = [] 
	minval = np.min(x)
	maxval = np.max(x)
	steps = float(maxval-minval)/200.0 
	for currentx in np.arange(minval,maxval+steps/2.,steps):
		x2.append(currentx)
		y2.append(calc_re(currentx,L))
	#
	if show and y:
		plt.plot(x2, y2,c=c)
		plt.show()
	#
	return x2,y2

def get_his(v,extent=None,bins=120,norm=1):
	extent = [v.min(),v.max()]
	a,b = np.histogram(v,bins=bins,range=extent)
	X = []; Y = [];
	for i in range(len(a)):
		X.append(float(b[i]+b[i+1])/2)
		Y.append(float(a[i])/2)
	#
	X = np.array(X)
	Y = np.array(Y)
	if norm:
		Y = Y/Y.max()
	return X,Y

if 0:
	# The radius of gyration (rg) and the end-to-end distance (re) correlated well with Miyazawa's d.
	# However, there is a scaling factor that appears to be required for perfect corrrelation, and that 
	# scaling factor will likely be correlated with $N$. Here, we generate a set of Rg, Re, and d, and 
	# try a bunch of scaling factors to see if one allows for perfect correlation.
	
	# First, we draw the relationship between the deviation in values such as Rg/d as a function of N
	step = 40
	show = 1 # if you want to know what is going on (graphically)
	
	ci     = -1
	colors = ['b','r','g','k']

	Ns = []
	rgs_over_d = []
	res_over_d = []
	
	columns = ['L','phi','psi','omega','d','theta','rg','re','R']
	df = pd.DataFrame(columns=columns)
	
	total_points = len(np.arange(-180,181,step))**2.0
	for L in np.arange(1,9)*5:
		L = float(L)
		ds  =  []
		rgs = []
		res = []
		
		geo = Geometry.geometry ("G")
		geo.phi     = 0.0
		geo.psi_im1 = 0.0
		geo.omega   = 0.0
		
		for omega in [0,180]:
			counter   = 0
			geo.omega = omega
			for phi in np.arange(-180,181,step):
				for psi in np.arange(-180,181,step):
					
					counter += 1
					if counter % 10 == 0.0 or counter == total_points:
						str_phi = "{0: >4}".format(int(phi))
						str_psi = "{0: >4}".format(int(psi))
						sys.stdout.write("\rL=%d (omega=%d) Completed %d percent\t" %(L,omega,int(100.*float(counter)/total_points)))
						sys.stdout.flush()
					
					geo.phi     = phi
					geo.psi_im1 = psi
					
					st = PeptideBuilder.initialize_res(geo)
					for i in range(1,int(L)):
						st = PeptideBuilder.add_residue(st, geo)
					
					d,theta,rs = locallib.calculate_d_theta_r(phi,psi,omega)
					
					rg = locallib.calculate_rg(st)
					re = locallib.calculate_re(st)
					
					# Ramachandran number from Mannige, Kundu, Whitelam, PLoS ONE, 2016
					R = float(phi+psi+360.0)/720.0 # locallib.R(phi=phi,psi=psi)
					
					df.loc[df.shape[0]] = [L,phi,psi,omega,d,theta,rg,re,R]
				#
			#
		#
	# 'Database' created	
	
	for xtype,ytype in [['d','rg'],['d','re']]:
		# http://stackoverflow.com/questions/19165259/python-numpy-scipy-curve-fitting
		from scipy.optimize import curve_fit
		def fitme(x,y):
			global L
			global xtype
			global ytype
			
			def fit_func(x, a):
				global L
				return a*x
			
			params = curve_fit(fit_func, x, y)
			
			p = tuple(params[0])
			
			print "For L="+str(L)+": \trg(d) = %f*x " %p
			if show:
				# Draw a scatter plot and the line (relationship)
				p = tuple(params[0])
				minval = np.min(x)
				maxval = np.max(x)
				x2 = []
				y2 = [] 
				steps = float(maxval-minval)/200.0 
				for currentx in np.arange(minval,maxval+steps/2.,steps):
					x2.append(currentx)
					y2.append(fit_func(currentx,*p))
				plt.scatter(x, y, s=20)
				plt.plot(x2, y2)
				plt.xlabel(xtype); plt.ylabel(ytype);
				plt.show()
			return p[0]
		
		Ls = []
		As = []
		for L in sorted(set(df['L'])):
			cdf = df[(df['L']==L)]
			a = fitme(np.abs(cdf[xtype]),cdf[ytype])
			Ls.append(L)
			As.append(a)
		#
		if show:
			plt.xlabel('Peptide length $L$'); plt.ylabel('a');
			plt.scatter(Ls,As)
			plt.show()
		#
		# There is a linear dependence of a on L, we find this out here:
		def fitme_line1(x,a):
			return a*x
		params = curve_fit(fitme_line1, Ls[2:], As[2:])
		print "The equation for %s as a function of L and %s is: %f L d (place that equation in 'calc_rg()' above)" %(ytype,xtype,params[0])
	

# ====================================================================================	
if 0: # FIGURE_START, shouldbeone
	# Show the relationship between (rg or re) with respect to d for various L.
	figname = "manuscript/automated_figures/fig_many_Ls.pdf"
	
	sns.axes_style("ticks")
	set_grays()
	
	# First, we draw the relationship between the deviation in values such as Rg/d as a function of N
	step = 25
	show = 1 # if you want to know what is going on (graphically)
	
	columns = ['L','phi','psi','omega','d','theta','rg','re','R']
	df = pd.DataFrame(columns=columns)
	
	xtype = 'd'
	
	textsize  = plt.rcParams['font.size']
	linewidth = plt.rcParams['lines.linewidth']
	
	
	# VERY IMPORTANT FOR GETTING THE SHAPE OF EACH PANEL CORRECT
	fig = plt.figure(figsize=(5,8))
	
	
	gs = mpl.gridspec.GridSpec(2, 1)#, height_ratios=[1,10,10])#, width_ratios =[1,1,1,1], height_ratios=[1,20,20])
	panels = ['rg','re']
	axes = { # 'legend': plt.subplot(gs[0,0]),
	         panels[0]: plt.subplot(gs[0,0]),
	         panels[1]: plt.subplot(gs[1,0])
	       }
	
	if not os.path.isfile(length_dependent_csv_file):
		prepare_main_csv_file()
	df = pd.read_csv(length_dependent_csv_file)
	# TO CREATE THE DATABASE FROM SCRATCH, COMMENT THE LINES BELOW
	'''
	# BUILDING THE 'DATABASE'
	total_points = len(np.arange(-180,181,step))**2.0
	for L in np.arange(1,3)*5:
		L = float(L)
		ds  =  []
		rgs = []
		res = []
		
		geo = Geometry.geometry ("G")
		geo.phi     = 0.0
		geo.psi_im1 = 0.0
		geo.omega   = 0.0
		
		for omega in [0,180]:
			counter   = 0
			geo.omega = omega
			
			for phi in np.arange(-180,181,step):
				for psi in np.arange(-180,181,step):
					
					counter += 1
					if counter % 10 == 0.0 or counter == total_points:
						str_phi = "{0: >4}".format(int(phi))
						str_psi = "{0: >4}".format(int(psi))
						sys.stdout.write("\rL=%d (omega=%d) Completed %d percent\t" %(L,omega,int(100.*float(counter)/total_points)))
						sys.stdout.flush()
					
					geo.phi     = phi
					geo.psi_im1 = psi
					
					st = PeptideBuilder.initialize_res(geo)
					for i in range(1,int(L)):
						st = PeptideBuilder.add_residue(st, geo)
					
					d,theta,rs = locallib.calculate_d_theta_r(phi,psi,omega)
					rg = locallib.calculate_rg(st)
					re = locallib.calculate_re(st)
					
					# Ramachandran number from Mannige, Kundu, Whitelam, PLoS ONE, 2016, modified by current paper
					R = float(phi+psi+360.0)/720.0 # locallib.R(phi=phi,psi=psi)
					df.loc[df.shape[0]] = [L,phi,psi,omega,d,theta,rg,re,R]
				#
			#
		#
	# 'Database' created	
	'''
	
	panel_letters = itertools.cycle(list(string.ascii_lowercase))
	for ytype in panels: #['rg','re']
		#palette       = itertools.cycle(color_palette())
		palette       = itertools.cycle(sns.color_palette("Paired", 20))
		
		marker_shapes = itertools.cycle(['o','s', '^','p'])
		ax = axes[ytype]
				
		for L in reversed(sorted(set(df['L']))): # Things drawn first show up at the 'top' of the graph (in layers)
			for omega in [180.0,0.0]:
				c  = next(palette) # place out of the "for omegea in [...]:" loop if the palette is not "Paired" 
				
				cdf = df[ (df['L']==L) & (df['omega']==omega) ]
				
				X = cdf[xtype]
				Y = cdf[ytype]
					
				if 1:
					# GETTING ERROR BARS AND INTERVALS
					xsteps = float(np.max(X)-np.min(X))/200
					Xbins  = np.arange(np.min(X)-xsteps/2, np.max(X)+xsteps/2, xsteps)
					xave   = []
					yave   = []
					ystd   = []
					for i in range(len(Xbins)-1):
						# temp df
						ys_in_range = cdf[(cdf[xtype] >= Xbins[i]) & (cdf[xtype] < Xbins[i+1])][ytype]
						if len(ys_in_range):
							xave.append((Xbins[i]+Xbins[i+1])/2.)
							yave.append(np.mean(ys_in_range))
							ystd.append(np.std(ys_in_range))
					
					xave = np.array(xave)
					yave = np.array(yave)
					ystd = np.array(ystd)
					
					cis_or_trans = 'trans'
					if omega==0.0:
						cis_or_trans = 'cis'
					label = '$L=%d \quad (%s)$' %(L,cis_or_trans)
					#ax.plot(xave, yave, color=c, label=label)
					#ax.errorbar(xave, yave, yerr=ystd, color=c, label=label)#'#CC4F1B')
					ax.fill_between(xave, yave - ystd, yave + ystd, interpolate=True, 
					                color=c,facecolor=c, alpha=0.7, label=label)
		
		linestyle = 'dashed'; #'solid','dashed','dashdot','dotted'
		for L in reversed(sorted(set(df['L']))):
			allXs = df[(df['L'] == L)][xtype]
			if ytype == 'rg':
				x,y = get_rg_estimates(allXs,show=0,L=L)
				ax.plot(x,y, color='k', linestyle=linestyle, linewidth=linewidth)
			elif ytype == 're':
				x,y = get_re_estimates(allXs,show=0,L=L)
				ax.plot(x,y, color='k', linestyle=linestyle, linewidth=linewidth)
		
		
		'''
		# loc keys
		best -- 0
		upper right -- 1
		upper left -- 2
		lower left -- 3
		lower right -- 4
		right -- 5
		center left -- 6
		center right -- 7
		lower center -- 8
		upper center -- 9
		center -- 10
		'''
		
		#loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=3, fancybox=True, shadow=True)
		if ytype == panels[0]:
			ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.4), ncol=3, fontsize=textsize*.9) #shadow=True, ncol=2,
		#ax.legend(loc=9,fontsize=textsize*1.1)
		
		xlabel = xtype
		ylabel = ytype
		if xlabel in type_to_label: 
			xlabel = type_to_label[xlabel]
		if ylabel in type_to_label: 
			ylabel = type_to_label[ylabel]
		
		#
		# Annotate the line for the largest polymer length L
		'''
			 \           /
			  \         /
			   \   l   p (arrow points at 'x', label starts at 'l')
			    \     / 
			     \   /
			      \ / 
			       .
		'''
		largestL = sorted(set(df['L']))[-1]
		
		largestL_Xs  = df[(df['L'] == largestL)][xtype]
		p_x = float(np.max(abs(largestL_Xs)) + np.min(abs(largestL_Xs)))/2.0
		p_y = calc_re(p_x,largestL)
		annotation = r'$%s \approx 2 \beta L |%s|$' %(ylabel.rstrip('$').lstrip('$'),xlabel.rstrip('$').lstrip('$'))
		if ytype == 'rg':
			p_y = calc_rg(p_x,largestL)
			annotation = r'$%s \approx \beta L |%s|$' %(ylabel.rstrip('$').lstrip('$'),xlabel.rstrip('$').lstrip('$'))
		#
		fig_xmin,fig_xmax = ax.get_xlim(); fig_ymin,fig_ymax = ax.get_ylim()  # getting axes min max values
		l_x = fig_xmin + 0.5*(fig_xmax-fig_xmin)
		l_y = fig_ymin + 0.8*(fig_ymax-fig_ymin)
		
		#weight = 'normal'
		weight = 'bold'
		c = 'k'
		
		c         = plt.rcParams['text.color']
		ax.annotate(annotation, xy=(p_x, p_y), xytext=(l_x, l_y), weight=weight, color=c, bbox=dict(boxstyle="round4", ec=c, fc="w"),
		            xycoords='data', arrowprops=dict(arrowstyle="->",color=c, lw=linewidth), fontsize=textsize*1.4,
		            horizontalalignment='center', verticalalignment='center')
		
		# AXES LABELS
		ax.set_xlabel(xlabel,fontsize=textsize*1.5); ax.set_ylabel(ylabel,fontsize=textsize*1.5);
		
		# PANEL LETTERING
		
		t = ax.set_title(r'('+next(panel_letters)+r')', size=textsize*1.3,loc='left', horizontalalignment='right')
		t.set_y(1.03)                           # Shifting the title up a little
		
	
	plt.tight_layout()#pad=0.0, w_pad=0.0, h_pad=0.0)
	sns.despine(offset=10)#, trim=True);
	
	plt.savefig(figname, dpi=180, bbox_inches='tight',transparent=True) #facecolor='w', edgecolor='w',
	if show_graphs:
        	os.system(pdf_viewer+" "+figname)

# ====================================================================================	
if 0: # FIGURE_START, shouldbeone
	# Well, this is not really a figure... but a nifty display by Symbolic Python (sympy)
	# which shows that the Ramachandran number, in the limit of infinite precision (sigma)
	# is actually SO MUCH SIMPLER!!!
	import sympy as sp
	phi,psi,sigma,L = sp.symbols('phi,psi,sigma,lambda')
	
	numerator = ((phi-psi+L)*sigma/sp.sqrt(2)) + (2*sigma*L/sp.sqrt(2)) * ((phi+psi+L)*sigma/sp.sqrt(2)) - (L*sigma/sp.sqrt(2))
	denominator = (2*L*sigma/sp.sqrt(2))**2
	R = numerator/denominator
	
	print "===================================================================================="
	print "The Ramachandran number presented by (Mannige, Kundu, Whitela, PLoS ONE, 2017) takes the (unrounded) form:"
	sp.pprint(R)
	print "===================================================================================="
	print "The simplified Ramachandran number presented here takes the form:"
	sp.pprint(sp.limit(R,sigma,sp.oo))
	print "===================================================================================="

if 1: # FIGURE_START, shouldbeone
	# Shows the Miyazawa equations for handedness (the hope is to find the extreme limits for d)
	import sympy as sp
	
	# All the symbols used in the equations below
	phi,psi,omega,sigma_n,sigma_a,sigma_c,v_na,v_ac,v_cn = sp.symbols('phi,psi,omega,sigma_n,sigma_a,sigma_c,v_na,v_ac,v_cn')
	
	# All angles in radians; all distances in Angstroms
	phi     = -sp.pi
	psi     = -sp.pi
	omega   =  sp.pi
	v_na    = 1.459
	v_ac    = 1.525
	v_cn    = 1.336
	sigma_a = np.radians(111.0)
	sigma_c = np.radians(117.2)
	sigma_n = np.radians(121.7)
	
	cos_theta_by_two = sp.cos((+phi+psi+omega)/2) * sp.sin(sigma_n/2)*sp.sin(sigma_a/2)*sp.sin(sigma_c/2) \
                          -sp.cos((+phi-psi+omega)/2) * sp.sin(sigma_n/2)*sp.cos(sigma_a/2)*sp.cos(sigma_c/2) \
                          -sp.cos((+phi+psi-omega)/2) * sp.cos(sigma_n/2)*sp.sin(sigma_a/2)*sp.cos(sigma_c/2) \
                          -sp.cos((-phi+psi+omega)/2) * sp.cos(sigma_n/2)*sp.cos(sigma_a/2)*sp.sin(sigma_c/2)
	
	theta = 2 * sp.acos(cos_theta_by_two)
	
	d_sin_theta_by_two  = (+v_na+v_ac+v_cn)*sp.sin((+phi+psi+omega)/2)*sp.sin(sigma_n/2)*sp.sin(sigma_a/2)*sp.sin(sigma_c/2) \
                            - (+v_na-v_ac+v_cn)*sp.sin((+phi-psi+omega)/2)*sp.sin(sigma_n/2)*sp.cos(sigma_a/2)*sp.cos(sigma_c/2) \
                            - (+v_na+v_ac-v_cn)*sp.sin((+phi+psi-omega)/2)*sp.cos(sigma_n/2)*sp.sin(sigma_a/2)*sp.cos(sigma_c/2) \
                            - (-v_na+v_ac+v_cn)*sp.sin((-phi+psi+omega)/2)*sp.cos(sigma_n/2)*sp.cos(sigma_a/2)*sp.sin(sigma_c/2)
	
	d = d_sin_theta_by_two/sp.sin(theta/2)
	sp.pprint(d)
	print '---------------------------'
	# The 'simplest' version of this value
	sp.pprint(sp.simplify(d))
	
	exit()
	



# ====================================================================================	
if 0: # FIGURE_START, shouldbeone
	# Study how the standard deviation between R(phi,psi,sigma) and R2(phi,psi) as we increase sigma:
	figname = "manuscript/automated_figures/fig_R_vs_R2.pdf"
	min_v = -360.0
	max_v =  360.0
	
	sns.axes_style("ticks")
	set_grays()
	
	fig = plt.figure(figsize=(10*.8,6*.8))
	ax  = plt.subplot(111)
	
	textsize  = 17.0
	linewidth = 2.0
	
	sigmas   = []
	averages = [] # Difference between R2 and R
	stds     = []
	
	astep = 5
	R2s = []
	for phi in np.arange(-180,180+astep/2,astep):
		for psi in np.arange(-180,180+astep/2,astep):
			R2 = (phi + psi - min_v)/(max_v-min_v)
			R2s.append(R2)
	R2s = np.array(R2s)
	
	sigma_vals = [1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000,10000000000,100000000000,1000000000000,10000000000000]
	for sigma in sigma_vals:
		Rs  = []
		for phi in np.arange(-180,180+astep/2,astep):
			for psi in np.arange(-180,180+astep/2,astep):
				Rs.append(locallib.R(phi,psi,sigma=sigma))
		diff =  np.abs(R2s - np.array(Rs))
		average, std = np.mean(diff),np.std(diff)
		sigmas.append(sigma)
		averages.append(average)
		stds.append(std)
	
	color = sns.color_palette()[0]	
	ax.plot(sigmas, averages, '--', linewidth=linewidth, c=color)
	ax.errorbar(sigmas[1:-1], averages[1:-1], yerr=stds[1:-1],fmt='o',
	            capthick=linewidth, elinewidth=linewidth, capsize=linewidth*6,markersize=10, c=color)
	
	sns.despine(offset=10)#, trim=True);
	
	ax.set_aspect(0.5)
	ax.set_yscale("log", nonposy='clip')
	ax.set_xscale("log", nonposy='clip')
	ax.set_xlabel('$\sigma$',fontsize=textsize*2)
	ax.set_ylabel(r'$|\mathcal{R}_1(\sigma)-\mathcal{R}_2|$',fontsize=textsize*1.5)
	plt.setp(ax.get_xticklabels(),     size=textsize)
	plt.setp(ax.get_yticklabels(), fontsize=textsize)#, rotation='vertical')
	
	# Dropping alternate tick labels because the text size is too large (but necessarily large)
	
	plt.setp(ax.get_xticklabels(), visible=False)
	plt.setp(ax.get_yticklabels(), visible=False)
	
	plt.setp(ax.get_xticklabels()[::2], visible=True)
	plt.setp(ax.get_yticklabels()[::3], visible=True)
	
	
	
	
	#
	ax.minorticks_off()
	ax.tick_params(axis=u'both', direction='out', width=linewidth, length=linewidth*4)#,labelsize=for)
	for spine in ax.spines.keys():
		ax.spines[spine].set_linewidth(linewidth)
	
	#print "R - (phi+psi-2pi)/4pi shows the following distribution: %f +/- %f" %()
	plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
	#gs.update(wspace=0.0, hspace=0)
	#plt.subplots_adjust(wspace=.201)
	#plt.subplots_adjust(hspace=0.001)
	
	plt.savefig(figname, dpi=180, bbox_inches='tight',transparent=True) #facecolor='w', edgecolor='w',
	if show_graphs:
        	os.system(pdf_viewer+" "+figname)

# ====================================================================================	
if 0: # FIGURE_START, shouldbeone
	figname = "manuscript/automated_figures/fig_various_rama_plots.pdf"
	
	sns.axes_style("ticks")#
	# VERY IMPORTANT FOR GETTING THE SHAPE OF EACH PANEL CORRECT
	plt.figure(figsize=(12*1.2,6*1.2))
	set_grays()
	
	
	textsize  = plt.rcParams['font.size']
	linewidth = plt.rcParams['lines.linewidth']*.5
	
	# Needs the length_dependent_csv_file (to be read in as the pandas dataframe 'df')
	# WARNING: VERY COSTLY TO REDO:
	if not os.path.isfile(length_dependent_csv_file):
		prepare_main_csv_file()
	
	df = pd.read_csv(length_dependent_csv_file)
	df = df[(df['L']==16.0)]
	
	'''              _______________
	                |_______________| colorbar (0,1)
	                                                
	      |--------|--------|--------|--------|
	trans |  rg    |  R     | theta  |   d    |
	      | (0,0)  | (0,1)  | (0,2)  | (0,3)  |
	      |--------|--------|--------|--------|
	cis   |  rg    |  R     | theta  |   d    |
	      | (1,0)  | (1,1)  | (1,2)  | (1,3)  |
	      |--------|--------|--------|--------|
	'''
	gs = mpl.gridspec.GridSpec(3, 4, width_ratios =[1,1,1,1], height_ratios=[1,20,20])
	gs.update(wspace=0.4, hspace=0.5)
	
	ax_cmap    = plt.subplot(gs[0,1:3])
	
	ax_trans_0 = plt.subplot(gs[1,0])
	ax_trans_1 = plt.subplot(gs[1,1])
	ax_trans_2 = plt.subplot(gs[1,2])
	ax_trans_3 = plt.subplot(gs[1,3])
	
	ax_cis_0 = plt.subplot(gs[2,0])
	ax_cis_1 = plt.subplot(gs[2,1])
	ax_cis_2 = plt.subplot(gs[2,2])
	ax_cis_3 = plt.subplot(gs[2,3])
	
	
	panels_a_through_c = ['rg','theta','d','R']
	axes = {'cmap' : ax_cmap,
	        'trans':{   panels_a_through_c[0]: ax_trans_0,
	                    panels_a_through_c[1]: ax_trans_1,
	                    panels_a_through_c[2]: ax_trans_2,
	                    panels_a_through_c[3]: ax_trans_3
	                },
	        'cis':  {   panels_a_through_c[0]: ax_cis_0,
	                    panels_a_through_c[1]: ax_cis_1,
	                    panels_a_through_c[2]: ax_cis_2,
	                    panels_a_through_c[3]: ax_cis_3
	                }
	       }
	
	cmap = sns.cubehelix_palette(100, start=.5, rot=-.75, dark=0.2, light=1, as_cmap=True)
	#cmap = sns.diverging_palette(240, 10, n=2, as_cmap=True)
	
	panel_letters = itertools.cycle(list(string.ascii_lowercase))
	
	amin = -180.0; amax=180.0;
	for omega in [180.0,0.0]:
		cis_or_trans = "trans"
		if omega == 0.0:
			cis_or_trans = "cis"
		
		xtype = "phi"; ytype = "psi";
		for ztype in panels_a_through_c: #axes['trans'].keys():
			ax = axes[cis_or_trans][ztype]
			
			tdf = df[(df['omega'] == omega)] # & (df['psi']>0)]
			x = tdf[xtype]
			y = tdf[ytype]
			
			levels = 10
			trunkated_levels = levels/2
			
			real_z_type = ztype
			if real_z_type == 'd':
				z = do_d(tdf[real_z_type])
			elif  real_z_type == 'd2':
				real_z_type = 'd'
				z = np.abs(tdf[real_z_type])
			elif real_z_type == 'theta':
				z = do_theta(tdf[real_z_type])
			else:
				z = tdf[real_z_type]
			
			if real_z_type == 're' or real_z_type == 'rg' or ztype == 'd2':
				levels = 7
			
			if ztype == 'd2':
				levels = 5
				
			if np.min(z) < 0:
				trunkated_levels += 1 
			
			'''
			# This proves that 'R(psi,psi)' ~ '(phi+psi+2pi)/4pi'
			if ztype == 'R':
				min_v = -360.0
				max_v =  360.0
				z = (tdf['phi'] + tdf['psi'] - min_v)/(max_v-min_v)
			'''
			X,Y,Z = locallib.xyz_to_ndarrays(x,y,z)
			
			if 0: # smooth functions if desired (it messes with contour lines sometimes at plot edges, though)
				Z = scipy.ndimage.filters.gaussian_filter(Z, 2)#, order=0)
			
			if 0:
				ax.imshow(Z, cmap='hot', interpolation='nearest')
			else:
				# DRAWING THE FILLED CONTOUR PLOT!
				ax.contourf(X, Y, Z, levels, cmap = cmap)
				
				# DRAWING LINES (negative contours will be dashed by default)
				CS = ax.contour(X, Y, Z, levels, colors=[plt.rcParams['lines.color']],linewidth=linewidth) 
				
				# ===============================================================================
				# ALL THIS STUFF IS FOR LABELING CONOTOURS (drop alternating labels, etc)
				# Recast level labels
				fmt   = {}
				for l in CS.levels:
					s = "%1.1f" %(round(l,1))
					if float(l) == float(int(l)):
						s = "%d" %(l)
					fmt[l] = s
				
				# Dropping every <skip_number> terms (ignored if skip_number == 0)
				counter = 0
				skip_number = 2
				if skip_number:
					for collection_i in range(len(CS.collections)):
						counter+=1
						if counter % skip_number == 0:
							plt.setp(CS.collections[collection_i], linewidth=0)
							fmt[CS.levels[collection_i]] = ''
				
				xmin,xmax = ax.get_xlim(); ymin,ymax = ax.get_ylim(); # getting axes min max values
				label_pos = []
				midpoint = np.array([float(xmax+xmin)/2.0, float(ymax+ymin)/2.0])
				for line in CS.collections:
					for path in line.get_paths():
						# find closest point
						X = []
						Y = []
						
						label_x = 0.0
						label_y = 0.0
						
						distance_from_center = []
						for point in path.vertices:
							distance_from_center.append(np.linalg.norm(point-midpoint))
						
						label_position_index = distance_from_center.index(np.min(distance_from_center))
						
						if 0:# ztype == 'd':
							plt.scatter(path.vertices[:,0],path.vertices[:,1])
							plt.scatter(*path.vertices[label_position_index],s=100,c='r')
							plt.ylim((-190,190))
							plt.xlim((-190,190))
							plt.show()
						
						label_pos.append(path.vertices[label_position_index])
				
				CLS = ax.clabel(CS, fontsize=textsize*.65, inline=1, fmt=fmt,manual=label_pos, colors='k') 
				
				if 1:
					# --------------------------------------------------------------------------
					# SOMETIMES, LABELS ENCROACH ON THE AXES AND GET CLIPPED. HERE WE FIND THOSE 
					# THAT EXCEED A THRESHOLD DISTANCE AND DELETE THEM
					# Swiped from: http://stackoverflow.com/questions/25873681/matplotlib-contour-plot-labels-overlap-axes
				
					thresh = 0.02  # ratio in x/y range in border to discard
					xmin,xmax = ax.get_xlim(); ymin,ymax = ax.get_ylim()  # getting axes min max values
					Dx = xmax-xmin
					Dy = ymax-ymin
				
					# check which labels are near a border
					keep_labels = []
					labels_to_delete = []
					for label in CLS:
						lx,ly = label.get_position()
						if xmin+thresh*Dx<lx<xmax-thresh*Dx and ymin+thresh*Dy<ly<ymax-thresh*Dy:
							# inlier, redraw it later
							keep_labels.append((lx,ly))
						else:
							labels_to_delete.append(label)
				
				
					# delete the original lines, redraw manually the labels we want to keep
					# this will leave unlabelled full contour lines instead of overlapping labels
					#for cline in CS.collections:
					#	cline.remove()
					for label in labels_to_delete:
						label.remove()
				
				# LABEL STUFF DONE, DONE. D.O.N.E.
				# ===============================================================================
			
			#'-' solid line style; '--' dashed line style; '-.' dash-dot line style; ':' dotted line style
			#ax.plot([-180,180],[180,-180], 'k--')  # -ve diagonal
			#ax.plot([-180,180],[0,0], 'k:')        # Horizontal line passing (0,0) 
			#ax.plot([0,0],[-180,180], 'k:')        # Vertical   line passing (0,0)
			xticks = range(int(amin),int(amax)+1,180)
			yticks = range(int(amin),int(amax)+1,180)
			# Setting ticks
			
			ax.tick_params(axis=u'both', direction='out', width=linewidth, 
			               length=linewidth*5, pad=0,labelsize=textsize)
			
			ax.set_xticks(xticks)
			plt.setp(ax.get_xticklabels(), fontsize=textsize*.7) #rotation='vertical',
			plt.setp(ax.get_yticklabels(), fontsize=textsize*.7) #rotation='vertical',
			ax.set_xlabel(type_to_label[xtype], size=textsize)
			
			# Drawing an equal-aspect-ratio graph
			#ax.set_aspect(1)
			
			# Setting title
			if 1:#cis_or_trans == 'trans':
				addto_title = r'('+next(panel_letters)+r')'
				t = ax.set_title(addto_title, size=textsize*1.05,loc='left') 
				ax.set_title(type_to_label[ztype], size=textsize*1.15,loc='center') 
				t.set_y(1.03)                           # Shifting the title up a little
			
			'''    
			      |--------|--------|--------|--------|
			trans |  re    |  R     | theta  |  |d|   |
			      | (0,0)  | (0,1)  | (0,2)  | (0,3)  |
			      |--------|--------|--------|--------|
			cis   |  re    |  R     | theta  |  |d|   |
			      | (1,0)  | (1,1)  | (1,2)  | (1,3)  |
			      |--------|--------|--------|--------|
			'''	
	
			ax.set_yticks(yticks)				
			if ztype != panels_a_through_c[0]:
				# erase the y labels
				ax.set_yticklabels([])
				ax.set_ylabel("", size=textsize)
			elif cis_or_trans == 'cis':
				ax.set_yticks(yticks)
				ax.set_ylabel(r"cis"+"\n\n"+type_to_label[ytype], size=textsize,weight='normal',style='italic')
			elif cis_or_trans == 'trans':
				ax.set_ylabel(r"trans"+"\n\n"+type_to_label[ytype], size=textsize,weight='normal',style='italic')
			
			
			if cis_or_trans != 'cis':
				# erase the x labels
				#ax.set_xticklabels([])
				ax.set_xlabel("", size=textsize)
			
			plt.sca(ax);plt.xticks(rotation=45)
			plt.sca(ax);plt.yticks(rotation=45)
			
			for spine in ax.spines.keys():
				ax.spines[spine].set_linewidth(linewidth*2.0)
	
	
	if 1:
		# addding the colorbar
		'''              _______________
			        |_______________| colorbar (0,1)
			        
		      |--------|--------|--------|--------|
		trans |  rg    |  re    |  R     |  |d|   |
		      | (1,0)  | (1,1)  | (1,2)  | (1,3)  |
		      |--------|--------|--------|--------|
		cis   |  rg    |  re    |  R     |  |d|   |
		      | (2,0)  | (2,1)  | (2,2)  | (2,3)  |
		      |--------|--------|--------|--------|
		'''	

		ax = axes['cmap'] #plt.subplot2grid((3,4), (0,1), colspan=2)
	
		# Making a fake 'figure' for a colorbar (easier than controling an orphan colorbar)
		x = [] ; y = [] ; z = [] ;
		for i in range(0,101,5): # to be x and y
			for j in [0,5,10]: # to be i and j
				x.append(i)
				y.append(j)
				z.append(i)
		x = np.array(x); y = np.array(y); z = np.array(z);
		X,Y,Z    = locallib.xyz_to_ndarrays(x,y,z)	
		
		# Plot the colorbar
		ax.contourf(X,Y,Z, 20,cmap = cmap)
		
		# Just setting the various aspects of the ticks, labels, etc.
		ax.tick_params(axis=u'both', direction='in', width=linewidth, length=0, pad=2,labelsize=textsize*0.8,
			       top=True, labeltop=True, bottom=False, labelbottom=False)
		ax.set_xticks([4,50,96])
		ax.set_xticklabels(['low','medium','high'])
		ax.set_yticks([])
		# Resetting the frame size
		for spine in ax.spines.keys():
			ax.spines[spine].set_linewidth(linewidth*2.0)
		#
	#
	
	#plt.tight_layout()#pad=0.0, w_pad=0.0, h_pad=0.2)
	#plt.subplots_adjust(wspace=.201)
	#plt.subplots_adjust(hspace=0.01)
	plt.savefig(figname, dpi=180, bbox_inches='tight',transparent=True) #facecolor='w', edgecolor='w',
	if show_graphs:
        	os.system(pdf_viewer+" "+figname)
	#plt.show()

if 0: # FIGURE_START, shouldbeone
	figname = "manuscript/automated_figures/fig_various_relationships.pdf"
	
	#sns.reset_orig()
	sns.set_style("ticks")
	set_grays()
	
	# VERY IMPORTANT FOR GETTING THE SHAPE OF EACH PANEL CORRECT
	plt.figure(figsize=(12,5))
	
	# DRAWING A BUNCH OF RELATIONSHIPS!
	colors = ['b','r','g','k']
	
	# Shows all other parameters to muck about with:
	#print c

	textsize  = plt.rcParams['font.size']
	linewidth = plt.rcParams['lines.linewidth']
	ytype     = 'rg'
	# =======================================================================
	# Here, we generate a set of relationships
	
	
	'''        colorbar
	      |---|---|---|---|
	cis   |0,0|0,1|0,2|0,3|
	      |---|---|---|---|
	trans |1,0|1,1|1,2|0,3|
	      |---|---|---|---|
	       phi  R theta d
	'''
	
	#plt.plot([0, 1], [0, 2], sns.xkcd_rgb["medium green"], lw=3)
	#plt.plot([0, 1], [0, 3], sns.xkcd_rgb["denim blue"]
	
	c1 = sns.hls_palette(2, l=.3, s=.6)[1] # should be blue
	c2 = sns.hls_palette(2, l=.8, s=.9)[0] # should be red
	c3 = sns.hls_palette(3, l=.8, s=.9)[1] # should be green
	relationships_palette = [c1,c2,c3] #sns.color_palette("colorblind")
	
	# Declaring the various panels
	gs = mpl.gridspec.GridSpec(3, 4, width_ratios =[1,1,1,1], height_ratios=[2,10,10])
	
	# Declaring the various axes (panels)
	row_no = 1
	ax_0_0 = plt.subplot(gs[row_no,0])#    Y 
	ax_0_1 = plt.subplot(gs[row_no,1],sharey=ax_0_0)
	ax_0_2 = plt.subplot(gs[row_no,2],sharey=ax_0_0)
	ax_0_3 = plt.subplot(gs[row_no,3],sharey=ax_0_0)
	#ax_0_4 = plt.subplot(gs[row_no,4],sharey=ax_0_0)
	row_no = 2                        #    Y             X
	ax_1_0 = plt.subplot(gs[row_no,0]              ,sharex=ax_0_0)
	ax_1_1 = plt.subplot(gs[row_no,1],sharey=ax_1_0,sharex=ax_0_1)
	ax_1_2 = plt.subplot(gs[row_no,2],sharey=ax_1_0,sharex=ax_0_2)
	ax_1_3 = plt.subplot(gs[row_no,3],sharey=ax_1_0,sharex=ax_0_3)
	#ax_1_4 = plt.subplot(gs[row_no,4],sharey=ax_1_0,sharex=ax_0_4)
	
	# Placing the axes in an easy to access dictionary of dictionaries (axes['cis/trans']['xtype'])
	# ytype is always 'rg' or 're' (or whatever you want it to be... it should be a label in our Pandas dataframe)
	axes = {'trans':{
	                   'phi': ax_0_0,
	                 'theta': ax_0_1,
	                     'd': ax_0_2,
	                     'R': ax_0_3
	                },
	        'cis':  {
	                   'phi': ax_1_0,
	                 'theta': ax_1_1,
	                     'd': ax_1_2,
	                     'R': ax_1_3
	                }
	}
	
	#
	step = 2
	colors  = ['b','r','g','k']
	lengths = np.arange(1,3)*8
	omegas  = [0,180]
	
	# COSTLY:
	if not os.path.isfile(length_dependent_csv_file):
		prepare_main_csv_file()
	df = pd.read_csv(length_dependent_csv_file)
	
	# Take the two largest lengths only (otherwise many plots get confusing, as all relationships are length dependent)
	Ls     = sorted(set(df['L']))[-2:]
	omegas = sorted(set(df['omega']))
	
	
	for cis_or_trans in axes.keys(): 
		omega = 180.0
		if cis_or_trans == 'cis':
			omega = 0.0
		for xtype in axes[cis_or_trans].keys():
			# We write to this panel:
			ax = axes[cis_or_trans][xtype]
			#
			# resetting the color palette
			palette = itertools.cycle(relationships_palette)
			#
			for L in Ls:
				c = next(palette)
				# 'cdf' stands for 'current dataframe'
				cdf = df[(df['L'] == L) & (df['omega'] == omega)]
				cdf = cdf.sort(xtype)#, ascending=False)
				
				m = 's'
				if omega == 0:
					m = 'o'
				
				X = cdf[xtype]
				if xtype == 'd':
					X = do_d(X)
				if xtype == 'theta':
					X = do_theta(X)
					
				Y = list(cdf[ytype])
				
				calcdf = pd.DataFrame({'X':X,'Y':Y})
				alpha = 0.5
				
				# GETTING ERROR BARS AND INTERVALS
				xsteps = float(np.max(X)-np.min(X))/200
				Xbins  = np.arange(np.min(X)-xsteps/2, np.max(X)+xsteps/2, xsteps)
				xave   = []
				yave   = []
				ystd   = []
				for i in range(len(Xbins)-1):
					# temp df
					ys_in_range = calcdf[(calcdf['X'] >= Xbins[i]) & (calcdf['X'] < Xbins[i+1])]['Y']
					if len(ys_in_range):
						xave.append((Xbins[i]+Xbins[i+1])/2.)
						yave.append(np.mean(ys_in_range))
						ystd.append(np.std(ys_in_range))
				
				xave = np.array(xave)
				yave = np.array(yave)
				ystd = np.array(ystd)
				
				ax.plot(xave, yave, 'k', color=c, label="$L="+str(L)+r"$")#'#CC4F1B')
				ax.fill_between(xave, yave - ystd, yave + ystd, interpolate=True, facecolor=c, alpha=0.3)
			#
		#
	#
	'''    phi psi  R theta d
	      |---|---|---|---|---|
	cis   |0,0|0,1|0,2|0,3|0,4|
	      |---|---|---|---|---|
	trans |1,0|1,1|1,2|0,3|0,4|
	      |---|---|---|---|---|
	'''
	#setting limits based on what is available
	for cis_or_trans in axes.keys(): # either rg or re
		for xtype in axes[cis_or_trans].keys():
			ax = axes[cis_or_trans][xtype]
			
			omega = 180.0
			if cis_or_trans == 'cis':
				omega = 0.0
						
			cdf = df[(df['omega']==omega)]
			xs = cdf[xtype]
			ys = cdf[ytype]
			
			if xtype == 'd':
				xs = do_d(xs)
				ax.set_xticks(range(-6,6,2))
			elif xtype == 'theta':
				xs = do_theta(xs)
				ax.set_xticks(range(0,361,90))
			elif xtype == 'R':
				ax.set_xticks(np.arange(0,1.001,0.5))
			elif xtype == 'phi' or xtype == 'psi':
				ax.set_xticks(range(-180,181,90))
			
			minx,miny,maxx,maxy = np.min(xs),np.min(ys),np.max(xs),np.max(ys)
			xpad_fraction = 0.15            ; ypad_fraction=0.1               ;
			xpad = (maxx-minx)*ypad_fraction; ypad = (maxy-miny)*ypad_fraction; 
			# X limits
			ax.set_xlim(minx-xpad,maxx+xpad)
			# Y limits
			ax.set_ylim(miny-ypad,maxy+ypad)
			
			linestyle = 'dashed'; #'solid','dashed','dashdot','dotted'
			if xtype == 'd':
				for L in Ls:
					allXs = np.abs(df[(df['L'] == L)][xtype])
					if ytype == 'rg':
						x,y = get_rg_estimates(allXs,show=0,L=L)
						ax.plot(x,y, color='k', linewidth=linewidth, linestyle=linestyle) 
					elif ytype == 're':
						#cdf = df[(df['L'] == L)|(df['omega'] == omega)]
						x,y = get_re_estimates(allXs,show=0,L=L)
						ax.plot(x,y, color='k', linewidth=linewidth, linestyle=linestyle)
			
			if cis_or_trans == 'cis': # then we keep all x-tick labels and such
				#ax.yaxis.get_label().set_verticalalignment("baseline")
				plt.sca(ax);plt.xticks(rotation=45)
				ax.set_xlabel(type_to_label[xtype], size=textsize*1.4)
				
			else:
				plt.setp(ax.get_xticklabels(), visible=False)
				
			if xtype == 'phi':  # then we keep all y-tick labels and such
				ax.set_ylabel(cis_or_trans+'\n\n'+type_to_label[ytype], size=textsize*1.2,style='italic')#weight='bold',
				plt.sca(ax);plt.yticks(rotation=45)
			else:
				plt.setp(ax.get_yticklabels(), visible=False)
			
			ax.tick_params(axis='both',direction='inout',length=linewidth*3)
			#ax.tick_params(axis='y',direction='out',length=linewidth*2)
			
			'''
			# If we are dealing with the left-most y-axes
			if xtype == 'phi':
				ax.tick_params(axis='y',direction='out',length=linewidth*2)
				ax.yaxis.set_ticks_position('left')
				ax.spines['left'].set_linewidth(linewidth)
				ax.spines['left'].set_color('k')
				ax.spines['right'].set_visible(False)
			# If we are dealing with the bottom-most x-axes
			if cis_or_trans == 'cis':
				ax.tick_params(axis='x',direction='out',length=linewidth*2)
				ax.xaxis.set_ticks_position('bottom')
				ax.spines['bottom'].set_linewidth(linewidth)
				ax.spines['bottom'].set_color('k')
				ax.spines['top'].set_visible(False)
			'''
			
			# Only show ticks on the left and bottom spines
			# Hide the right and top spines
			#ax.spines['right'].set_visible(False)
			#ax.spines['top'].set_visible(False)
			#ax.spines['bottom'].set_visible(False)
	# Here, we should reset x label positions, if possible (the tick labels seem to unalign each x-axis's labels)
	#for xtype in axes['cis'].keys():
	#	ax = axes[cis_or_trans][xtype]
	# ======================================================================
	# Create custom legend
	with sns.axes_style("white"):
		legend_axis = plt.subplot(gs[0,1:3])
		legend_axis.spines['left'].set_visible( False)
		legend_axis.spines['right'].set_visible( False)
		legend_axis.spines['top'].set_visible(   False)
		legend_axis.spines['bottom'].set_visible(False)
		legend_axis.yaxis.set_ticks_position('none')
		legend_axis.xaxis.set_ticks_position('none')
		plt.setp(legend_axis.get_xticklabels(), visible=False)
		plt.setp(legend_axis.get_yticklabels(), visible=False)
		#legend_axis.xaxis.set_major_locator(mpl.ticker.NullLocator())
		#legend_axis.yaxis.set_major_locator(mpl.ticker.NullLocator())
		h,l = axes['trans']['d'].get_legend_handles_labels() # get labels and handles from ax1
		present_y = 0.0; padding = 0.1;
		for line,label in zip( h , l ):
			c = line.get_color()
			legend_axis.plot([present_y,present_y+0.5],[0.5,0.5],color=c)
			present_y += 0.5 + padding
			legend_axis.text(present_y, 0.5, label, horizontalalignment='left',verticalalignment='center',
				        fontsize=20, size=textsize*1.2, color=c)
			present_y += padding*4
		legend_axis.axis([0-padding*5,present_y+padding*5,0,1])
		legend_axis.set_axis_bgcolor('none')
	#
	# ======================================================================
	# Setting panel titles as panel labels
	panel_letters = itertools.cycle(list(string.ascii_lowercase))
	for xtype in ['phi','R','theta','d']: #axes['trans'].keys():
		ax = axes['trans'][xtype]
		addto_title = r'('+next(panel_letters)+r')'+"\n"
		ax.set_title(addto_title,size=textsize*1.1,loc='left')
	# ======================================================================
	
	plt.tight_layout()#pad=0.0, w_pad=0.0, h_pad=0.0)
	plt.subplots_adjust(hspace=0,wspace=0) 
	plt.savefig(figname, dpi=180, bbox_inches='tight',transparent=True) #facecolor='w', edgecolor='w',
	if show_graphs:
        	os.system(pdf_viewer+" "+figname)
	# =======================================================================


if 0: # FIGURE_START, shouldbeone
	
	np.random.seed(931123231)
	
	figname = "manuscript/automated_figures/fig_rama_intro.pdf"
	sns.reset_orig()
	
	sns.set_style("ticks")
	#sns.set_style("darkgrid")
	set_grays()
	
	textsize = 12.0
	linewidth = 1.0
	alpha     = 0.5
	
	# VERY IMPORTANT FOR GETTING THE SHAPE OF EACH PANEL CORRECT
	plt.figure(figsize=(16,6))
	
	# Plot distributions of secondary structures in various formats
	sstype = 'segno' # check <ss_codes>
	
	#                      nrows  ncols
	#                          |  |
	gs = mpl.gridspec.GridSpec(1, 3, 
	                         height_ratios = [1],
	                         width_ratios  = [1,1,1])
	
	# axes[0] will host the actual secondary structure distribution
	# axes[1] will display a simulated trajectory in Ramachandran space
	axes = [plt.subplot(gs[0,0]),plt.subplot(gs[0,1]),plt.subplot(gs[0,2])]
	
	# -------------------------------------------------
	# LOADING THE SECONDARY STRUCTURE DATABASE
	# First checking to see if it exists
	if not os.path.isfile(master_ss_file):
		 # Create one if it does not exist:
		prepare_master_ss_file()
	# Loading:
	df = pd.read_csv(master_ss_file)
	# Cleaning:
	df = df[(df['d'] != 999)]
	df = df[(df['phi'] != 999)]
	df = df[(df['psi'] != 999)]
	df = df[(df['omega'] != 999)]
	ss_to_study = ['all','G','E','P','H']
	# -------------------------------------------------
	
	normal_levels  = [0.50] # Levels for secondary structure contours ('G','E','P','H')
	levels_for_all = [0.92] # Levels for all amino acids ('all')
	
	
	panel_letters = itertools.cycle(list(string.ascii_lowercase))
	
	kb = 0.001987204118 # Boltzmann constant in kcal/mol
	def treat_angle(a):
		return np.int(round(a,0))
	
	def calc_energy(vals,pfunc,T=298.0):
		total_energies = []
		new_vals = []
		if type(vals) is int or type(vals) is float:
			new_vals = [vals] 
		else:
			new_vals = vals
		for val in new_vals:
			total_energies.append( kb * T * pfunc[val] * np.log( pfunc[val] ) )
			#                     Boltzmann k
		return total_energies
	
	column = 0
	# Creating a palette through which we cycle for various secondary structures	
	palette_colors = sns.hls_palette(len(ss_to_study), l=.4, s=.7)
	palette = itertools.cycle(palette_colors)
	
	
	phis,psis = df['phi'],df['psi']
	
	'''
	# Obtaining probability functions using Kernel Density Estimation (KDE)
	pdfs = {'phi':{},'psi':{}}
	smooth_pdf_phi = scipy.stats.gaussian_kde(list(phis))
	smooth_pdf_psi = scipy.stats.gaussian_kde(list(psis))
	print "populating PDFs"
	import time
	t1 = time.time()
	for a in np.arange(-180,181,1):
		pdfs['phi'][treat_angle(a)] = smooth_pdf_phi(treat_angle(a))[0]
		pdfs['psi'][treat_angle(a)] = smooth_pdf_psi(treat_angle(a))[0]
	t2 = time.time()
	print "\tdone (took %f seconds)" %(t2-t1)
	print "WRITING THE PDFS TO: 'pdfs.py'"
	f = open('local_imports/pdfs.py','w')
	f.write( 'pdfs = %s' %(str(pdfs).replace('{','{\n').replace('}','}\n')))
	f.close()
	'''
	import pdfs
	
	seq = ''
	seq += 'HHHH'
	seq += 'EEEE'
	seq += 'PPPP'
	
	palette_choices = sns.color_palette('deep',len(seq))
	palette_choices.pop(1)
	palette = itertools.cycle(palette_choices)
	# Collecting statistics on the relevant states
	angles_ave = {}
	angles_std = {}
	for ss in set(seq):
		#
		tdf = df[(df['segno']==ss_codes['segno'][ss])]
		#
		phis = tdf['phi']; phi_ave = np.average(phis); phi_std = np.std(phis); 
		psis = tdf['psi']; psi_ave = np.average(psis); psi_std = np.std(psis); 
		#
		angles_ave[ss] = {'phi':{},'psi':{}}
		#
		angles_ave[ss]['phi']['ave'] = phi_ave
		angles_ave[ss]['phi']['std'] = phi_std/2.
		#
		angles_ave[ss]['psi']['ave'] = psi_ave
		angles_ave[ss]['psi']['std'] = psi_std/2.
	
	# Setting up the initial states of each amino acid of the sequence "seq"
	peptide_vals = {'phi':[],'psi':[]}
	for ss in seq:
		phi = np.random.uniform(angles_ave[ss]['phi']['ave']-angles_ave[ss]['phi']['std'], angles_ave[ss]['phi']['ave']+angles_ave[ss]['phi']['std'])
		psi = np.random.uniform(angles_ave[ss]['psi']['ave']-angles_ave[ss]['psi']['std'],angles_ave[ss]['psi']['ave']+angles_ave[ss]['psi']['std'])
		#
		peptide_vals['phi'].append(treat_angle(phi))
		peptide_vals['psi'].append(treat_angle(psi))
	#
	# Calculating initial energies
	peptide_energies = {}
	for p in ['phi','psi']:
		peptide_energies[p] = calc_energy(vals=peptide_vals[p],pfunc=pdfs.pdfs[p])
	#
	
	#axes[1].scatter(peptide_vals['phi'],peptide_vals['psi'],c='b',s=100)
	#axes[1].axis([-180,180,-180,180])
	
	# Making perturbations:
	phi_or_psi       = ['phi','psi']
	residue_position = range(len(seq))
	all_phi_vals = []
	all_psi_vals = []
	steps = 100000
	for trial in range(steps+1):
		p  = np.random.choice(phi_or_psi)
		ri = np.random.choice(residue_position)
		
		old_energy = peptide_energies[p][ri]
		dval       = np.random.choice([-1,1])
		
		old_val = peptide_vals[p][ri]
		
		# Calculating new val
		new_val = old_val+dval
		# Making sure that newval is within bounds
		new_val = (new_val + 180)%360 - 180
		# treating new val
		new_val = treat_angle(new_val)
		
		new_energy = calc_energy(vals=new_val,pfunc=pdfs.pdfs[p])[0]
		
		#print old_energy, new_energy, old_energy < new_energy
		keep = 0
		
		if trial % (steps/10) == 0:
			print trial
			all_phi_vals.append(copy.deepcopy(peptide_vals['phi']))
			all_psi_vals.append(copy.deepcopy(peptide_vals['psi']))
		
		
		Tmc = 200.0
		ediff = new_energy - old_energy # Negative if the new energy is lower than the old energy
		#ediff = pdfs.pdfs[p][old_val] - pdfs.pdfs[p][new_val]
		if ediff < 0:
			keep = 1
		else:
			# Metropolis criterion
			if np.random.random() < np.exp(-ediff/(kb*Tmc)):
				keep = 1
		
		if keep: #pdfs.pdfs[p][new_val] > pdfs.pdfs[p][old_val]:# keep:
			peptide_vals[p][ri]     = new_val
			peptide_energies[p][ri] = new_energy
		
		#print trial #, np.sum(peptide_energies['phi']+peptide_energies['psi'])
		
	all_phi_vals = np.array(all_phi_vals)
	all_psi_vals = np.array(all_psi_vals)
	
	marked_residues = [3,8] # These are the lines to color and thicken
	for resid in range(all_phi_vals.shape[1]):
		current_phis = list(all_phi_vals[:,resid])
		current_psis = list(all_psi_vals[:,resid])
		
		
		boundaries = []
		
		if 1: # Then only draw contiguous regions with one line
			for i in range(1,len(current_phis)):
				if np.abs(current_phis[i-1] - current_phis[i]) > 180: # checking if an abrubt jump happned
					boundaries.append(i)
				if np.abs(current_psis[i-1] - current_psis[i]) > 180: # checking if an abrubt jump happned
					boundaries.append(i)
		
		if not len(boundaries):
			boundaries = [-1]
			
		boundaries = sorted(set(boundaries))
		
		c = 'Gray'; alpha = 0.4; lw = 2;
		#c = next(palette);
		if resid in marked_residues:
			c = next(palette); alpha = 1.0; lw = 3;
		
		# Draw a starting point, if desired
		axes[1].scatter(current_phis[1],current_psis[1],c=c,alpha=alpha,s=100)
		
		# Drawing the 'trajectory'
		previous_boundary = 1 # Starting from the second state of the residue (first state is too 'artificial')
		for boundary in boundaries:
			axes[2].plot(current_phis[previous_boundary:boundary],
			             current_psis[previous_boundary:boundary],c=c,alpha=alpha,lw=lw)
			previous_boundary = boundary
	if 0:
		axes[1].set_aspect(1)
		axes[2].set_aspect(1)	
		axes[1].axis([-180,180,-180,180])
		axes[2].axis([-180,180,-180,180])
		axes[1].set_yticks([-180,0,180])
		axes[2].set_yticks([-180,0,180])
		axes[1].set_xticks([-180,0,180])
		axes[2].set_xticks([-180,0,180])
		plt.show()
		exit()
	
	#
	axis_to_average_x = {}; # axis_to_average_x['theta']=[]; axis_to_average_x['d']=[]; 
	axis_to_average_x['R']=[]
	for ss in ss_to_study:
		levels    = normal_levels
		if ss == 'all':
			levels = levels_for_all
		
		xtype = 'phi'; ytype = 'psi';
		
		cdf = df.copy()
		fill_contour_lines = 1
		draw_contour_lines = 0
		smooth             = 2
		cmap = "#D3D3D3"  # This is the default color for an 'all' map
		c    = cmap # Uses only a color if cmap is a string
		if ss != "all":
			# ----------------
			# NEW COLOR SCHEME:
			c1 = [1,1,1];	c2 = next(palette);
			cdict = {#                c1                   c2
				'red':   ((0.00,  c1[0], c1[0]), (1.0, c2[0], c2[0])), 
				'green': ((0.00,  c1[1], c1[1]), (1.0, c2[1], c2[1])),
				'blue':  ((0.00,  c1[2], c1[2]), (1.0, c2[2], c2[2])) }
			current_cmap = LinearSegmentedColormap('tmpcmap', cdict)
			c = current_cmap(0.99999)
			# ----------------
			draw_contour_lines = 1
			smooth             = 2
			cmap=current_cmap
			# Here, if so inclined, you could put a check to transform the phi psi values
			if ss in ss_codes[sstype]:
				cdf = cdf[(cdf[sstype] == ss_codes[sstype][ss])]
		
		phis = list(cdf['phi']); psis = list(cdf['psi']); 
		
		ss_to_fill_collor =  {'all':'#e5e5e5','H':'#c5c5c5','E':'#c5c5c5','G':'#c5c5c5','P':'#c5c5c5'}
		ss_to_line_collor =  {'all':'#838383','H':'#838383','E':'#838383','G':'#838383','P':'#838383'}
		
		cmap = ss_to_fill_collor[ss]
		c    = ss_to_fill_collor[ss]
		
		all_size     = len(df)
		current_size = len(cdf)
		print ss,'\t',current_size,'\t',100.0*float(current_size)/all_size
		
		# Plotting contours on the ramachandran plot
		# Fill 
		locallib.plot_ramachandran_histogram(phis,psis,levels=levels,cmap=cmap,
			                    fill=1,lines=0, ax=axes[0], 
			                    smooth=smooth, linestyles='dashed',linewidths=1)
		
		if ss != 'all':
			# Drawing Lines
			line_width_modifier = 1.5
			#
			locallib.plot_ramachandran_histogram(phis,psis,levels=[sorted(levels)[-1]],cmap=ss_to_line_collor[ss],
					            fill=0,lines=1, ax=axes[0], 
					            smooth=smooth, linestyles=['solid'], linewidths=linewidth*line_width_modifier)
			'''
			locallib.plot_ramachandran_histogram(phis,psis,levels=[sorted(levels)[-1]],cmap=ss_to_line_collor[ss],
					            fill=0,lines=1, ax=axes[1], 
					            smooth=smooth, linestyles=['dotted'], linewidths=linewidth*line_width_modifier)
			'''
			#
			# Annotate 
			xave = np.median(phis)
			yave = np.median(psis)
			ss_label = str(ss)
			weight = 'normal'
			if ss_label in ss_name_to_label:
				ss_label = ss_name_to_label[ss_label]
			
			axes[0].annotate(ss_label, xy=(xave, yave), xytext=(xave+30, yave-30), weight=weight, 
				color=ss_to_line_collor[ss], xycoords='data', fontsize=textsize*1, 
				arrowprops=dict(arrowstyle="-",color=ss_to_line_collor[ss], lw=linewidth))
		else:
			xave = -90
			yave =  90
			weight = 'normal'
			axes[0].annotate('All', xy=(xave, yave), xytext=(xave+30, yave-30), weight=weight, 
				color=ss_to_line_collor[ss], xycoords='data', fontsize=textsize*1, ha='left',
				arrowprops=dict(arrowstyle="-",color=ss_to_line_collor[ss], lw=linewidth))
	
	
	# (Re)setting labels and titles
	for column in [0,1,2]:
		ax = axes[column]
		
		ax.axis([-180,180,-180,180])
		ax.set_yticks([-180,0,180])
		ax.set_xticks([-180,0,180])
		ax.set_yticks([-180,0,180])
		ax.set_xlabel('$\phi$')
		ax.set_ylabel(' $\psi$')
		plt.sca(ax);plt.xticks(rotation=45)
		plt.sca(ax);plt.yticks(rotation=45)
		ax.set_aspect(1)
		
		for spine in ax.spines.keys():
			ax.spines[spine].set_linewidth(linewidth*1.6)
		
		addto_title = r'('+next(panel_letters)+r')'+"\n"
		t = ax.set_title(addto_title, size=textsize*1.4, loc='left',horizontalalignment='right')
		t.set_y(1.03) # Shifting the title up a little
		
		# axis tick labels
		plt.setp(ax.get_yticklabels(), size=textsize*1)#, rotation="vertical")
		plt.setp(ax.get_xticklabels(), size=textsize*1)
		# axis label (name of the axis)
		ax.set_xlabel(ax.get_xlabel(), fontsize=textsize*2)
		ax.set_ylabel(ax.get_ylabel(), fontsize=textsize*2, rotation="horizontal", ha='center', va='center',labelpad=20)
		#
						
	plt.tight_layout()#pad=0.0, w_pad=0.0, h_pad=)
	#gs.update(wspace=1, hspace=1)
	#plt.subplots_adjust(wspace=.201)
	#plt.subplots_adjust(hspace=0.001)
	plt.savefig(figname, dpi=180, bbox_inches='tight',transparent=True) #facecolor='w', edgecolor='w',
	if show_graphs:
        	os.system(pdf_viewer+" "+figname)

#
#
# 
if 1: # FIGURE_START, shouldbeone
	figname = "manuscript/automated_figures/fig_ss_2d_1d.pdf"
	sns.reset_orig()
	
	#sns.set_context("notebook")#, font_scale=1)
	sns.set_style("ticks")
	set_grays()
	
	textsize = 12.0
	linewidth = 1.0
	alpha     = 0.5
	annotation_color = 'teal'
	# VERY IMPORTANT FOR GETTING THE SHAPE OF EACH PANEL CORRECT
	plt.figure(figsize=(6.5,9))
	
	# Plot distributions of secondary structures in various formats
	sstype = 'segno' # check <ss_codes>
	
	# ----------------
	#                      nrows  ncols
	#                          |  |
	gs = mpl.gridspec.GridSpec(2, 1, 
	                         height_ratios = [5,0.5],#,1,1],
	                         width_ratios  = [1])
	axes = {}
	#    Column
	#    |
	axes[0] = {'rama': plt.subplot(gs[0,0]),
	              'R': plt.subplot(gs[1,0])
	          #'theta': plt.subplot(gs[2,0]),
	          #    'd': plt.subplot(gs[3,0])
	          }
	'''
	axes[1] = {'rama': plt.subplot(gs[0,1]),
	              'R': plt.subplot(gs[1,1])
	          #'theta': plt.subplot(gs[2,1]),
	          #    'd': plt.subplot(gs[3,1])
	          }
	'''
	draw_d_contours = 0
	if draw_d_contours:
		# we will draw a d=0 contour line in each plot that has the name 'rama'
		d0df = pd.read_csv(length_dependent_csv_file)
		d0df = d0df[(d0df['L']==8.0) & (d0df['omega']==180.0)]
	
	mark_important_secondary_structures = [4,8,6,11]		
	
	panel_letters = itertools.cycle(list(string.ascii_lowercase))
	for column in [0]: #,1]:
		
		# -----------------
		# Loading up the SS database
		if column == 0:
			normal_levels = [0.3333,0.6666] # For ramachandran plots, given most ss keys
			levels_for_all = [0.92]         # Same, but in case the 'all' key is used for ss
			
			# Loading:
			if not os.path.isfile(master_ss_file):
				prepare_master_ss_file()
			df = pd.read_csv(master_ss_file)
			# Cleaning:
			df = df[(df['d'] != 999)]
			ss_to_study = ['all','G','E','P','H']
			
		
		if column == 1:
			normal_levels = [0.6666,0.9] # For ramachandran plots, given most ss keys
			if not os.path.isfile(fake_ss_file):
				write_fake_secondary_structure_database(sstype = 'segno',omega=180.0,sample_points=50000, sigma=5)
			df = pd.read_csv(fake_ss_file)
			# Lets study ALL fake secondary structures
			ss_to_study = sorted(set(list(df[sstype])))
			
		'''
		ds_     = np.abs(df[(df[sstype] == ss_codes[sstype]['E'])]['d'])
		thetas_ = df['theta']
		Rs_     = df['R']
		mind    ,maxd     =     ds_.min() ,     ds_.max()
		minR    ,maxR     =     Rs_.min() ,     Rs_.max()
		minTheta,maxTheta = thetas_.min() , thetas_.max()
		yextents_dict = {'d':[1,maxd],'R':[0.3,0.67],'R2':[0.3,0.67],'R2':[0.3,0.67],'theta':[minTheta,maxTheta]}
		'''
		
		# Creating a palette through which we cycle for various secondary structures	
		#palette_colors = sns.dark_palette("red", len(ss_to_study))
		palette_colors = sns.color_palette("cubehelix", len(ss_to_study))
		#palette_colors = sns.color_palette("RdBu", len(ss_to_study))
		#if column == 0:
		#palette_colors = sns.color_palette("colorblind",len(ss_to_study))
		palette_colors = sns.hls_palette(len(ss_to_study), l=.4, s=.7)
		
		palette = itertools.cycle(palette_colors)
		
		if draw_d_contours:
			# draw d=0 contour lines
			X = d0df['phi']
			Y = d0df['psi']
			Z = d0df['d']
			
			X,Y,Z    = locallib.xyz_to_ndarrays(X,Y,Z)
			Z = scipy.ndimage.filters.gaussian_filter(Z, 2)#, order=0)
			CS = axes[column]['rama'].contour(X, Y, Z, [0], colors='gray', linewidths=linewidth, linestyles='dashdot')
			if 0: # Draw contour labels
				# Recast levels
				fmt = {}
				for l in CS.levels:
					s = "%1.3f" %(l)
					if float(l) == float(int(l)):
						s = "d=%d" %(l)
					fmt[l] = s
				#
				axes[column]['rama'].clabel(CS, fontsize=textsize*0.9, inline=1, fmt=fmt) #fmt = contour_label_format)
			#
		#
		axis_to_average_x = {}; axis_to_average_x['theta']=[]; axis_to_average_x['d']=[]; axis_to_average_x['R']=[];
		for ss in ss_to_study:
			levels    = normal_levels
			if ss == 'all':
				levels = levels_for_all
			
			xtype = 'phi'; ytype = 'psi';
			
			cdf = df.copy()
			fill_contour_lines = 1
			draw_contour_lines = 0
			smooth             = 2
			cmap = "#D3D3D3"  # This is the default color for an 'all' map
			c    = cmap # Uses only a color if cmap is a string
			if ss != "all":
				# ----------------
				# NEW COLOR SCHEME:
				c1 = [1,1,1];	c2 = next(palette);
				cdict = {#                c1                   c2
					'red':   ((0.00,  c1[0], c1[0]), (1.0, c2[0], c2[0])), 
					'green': ((0.00,  c1[1], c1[1]), (1.0, c2[1], c2[1])),
					'blue':  ((0.00,  c1[2], c1[2]), (1.0, c2[2], c2[2])) }
				current_cmap = LinearSegmentedColormap('tmpcmap', cdict)
				c = current_cmap(0.99999)
				# ----------------
				draw_contour_lines = 1
				smooth             = 2
				cmap=current_cmap
				# Here, if so inclined, you could put a check to transform the phi psi values
				if ss in ss_codes[sstype]:
					cdf = cdf[(cdf[sstype] == ss_codes[sstype][ss])]
				else:
					# Then this may be a 'fake' secondary structure. No translation needed, just use the ss code.
					cdf = cdf[(cdf[sstype] == ss)]
				
			all_size     = len(df)
			current_size = len(cdf)
			print ss,'\t',current_size,'\t',100.0*float(current_size)/all_size
				
			
			# Plotting contours on the ramachandran plot
			# Fill
			locallib.plot_ramachandran_histogram(cdf['phi'],cdf['psi'],levels=levels,cmap=cmap,
				                    fill=1,lines=0, ax=axes[column]['rama'], 
				                    smooth=smooth, linestyles=['solid', 'dotted'],linewidths=0,
				                    show_diag=1)
			
			if ss != 'all':
				# Lines
				line_width_modifier = 1.3
				if ss in mark_important_secondary_structures:
					line_width_modifier = 2.3
				#
				locallib.plot_ramachandran_histogram(cdf['phi'],cdf['psi'],levels=[sorted(levels)[-1]],cmap='k',#cmap,
						            fill=0,lines=1, ax=axes[column]['rama'], 
						            smooth=smooth, linestyles=['solid'], linewidths=linewidth*line_width_modifier)
				
				#
				# Annotate 
				xave = np.median(list(cdf['phi']))
				yave = np.median(list(cdf['psi']))
				ss_label = str(ss)
				weight = 'normal'
				if ss in mark_important_secondary_structures:
					weight = 'bold'
				if ss_label in ss_name_to_label:
					ss_label = ss_name_to_label[ss_label]
				axes[column]['rama'].annotate(ss_label, xy=(xave, yave), xytext=(xave+30, yave-30), weight=weight, color=c,
				            xycoords='data', arrowprops=dict(arrowstyle="-",color=c, lw=linewidth), fontsize=textsize*0.9)
				     
			def get_ave(vals):
				X,Y = get_his(vals,norm=0); X=list(X); Y=list(Y)
				return X[Y.index(np.max(Y))]
			
			
			# We study the ss distributions as a function of R, theta, d
			if ss != "all":
				bins = 150
				if column == 1:
					bins = 50
				
				# r
				vals = cdf['R']
				X,Y = get_his(vals,bins=bins)
				axes[column]['R'].fill_between(X,Y,0,facecolor=c,alpha=alpha)
				axes[column]['R'].plot(X,Y,c='k',linewidth=linewidth)
				axis_to_average_x['R'].append([ss,get_ave(vals),c])
				

				'''
				# theta
				vals = do_theta(cdf['theta'])
				#vals = abs(cdf['d'])+(cdf['phi']+cdf['psi'] + 360.)/(720.)
				if 0: # rescale values wrt min and max vals
					mind = min(vals)
					#maxd = max(vals)
					#for i in range(len(vals)):
					#	if vals[i] < 0:
					#	vals[i] = vals + mind 
					vals = vals - mind
				
				X, Y = get_his(vals,bins=bins)
				axes[column]['theta'].fill_between(X,Y,0,facecolor=c,alpha=alpha)
				axes[column]['theta'].plot(X,Y,c='k',linewidth=linewidth)
				axis_to_average_x['theta'].append([ss,get_ave(vals),c])
				
				# d
				tdf = cdf.copy()
				vals = do_d(tdf['d'])
				if 0: # rescale values wrt min and max vals
					mind = min(vals)
					maxd = max(vals)
					tdf.ix[tdf.d<=0, 'd'] = maxd+(tdf.ix[tdf.d <= 0, 'd']-mind)
					vals = tdf['d']*-1*np.sign(np.cos(np.radians(tdf['theta'])/2.0))
				
				#vals = tdf['d']*-1.0*np.cos(np.radians(tdf['theta'])/2.0)
				#vals = np.sign(tdf['d'])*np.sin(np.radians(tdf['theta'])) # this is the h number from Mannige, PeerJ, 2017
				
				X,Y = get_his(vals,bins=bins)
				axes[column]['d'].fill_between(X,Y,0,facecolor=c,alpha=alpha) 
				axes[column]['d'].plot(X,Y,c='k',linewidth=linewidth)
				axis_to_average_x['d'].append([ss,get_ave(vals),c])
			
				# Example of how to plot vertically
				#axes[column]['d'].fill_betweenx(X,Y,0,facecolor=c,alpha=alpha) 
				#axes[column]['d'].plot(Y,X,c='k',linewidth=linewidth)
				'''

		# setting the ramachandran plot aspect ratio
		if 1:#for name in ['d','theta','R']: #axes.items():	
			axes[column]['rama'].set_aspect(1)
		
		# Making some tweaks to all but the 
		for name in ['R']: #['d','theta','R']: #axes.items():	
			ax = axes[column][name]
		
			xmin,xmax = ax.get_xlim(); ymin,ymax = ax.get_ylim()  # getting axes min max values
			
			ax.set_yticks([0,1.5])
			
			ax.set_xlabel(type_to_label[name])
		
			# Only show ticks on the left and bottom spines
			ax.xaxis.set_ticks_position('bottom')
			ax.yaxis.set_ticks_position('none')
			# Hide the right and top spines
			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)
			ax.spines['left'].set_visible(False)
			
			xpad =  (xmax-xmin)*0.0; ypad =  (ymax-ymin)*0.0;    # setting padding
			
			ymin = 0.0
			extent = [xmin, xmax+xpad, ymin, ymax+ypad+0.6] # storing (padded) axis min max values
			xscaling = .094                                          # setting how much larger must the x axis be compared to the y axis
			ax.axis(extent)                                       # setting (padded) axis min max values
			#ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/xscaling)
			
			# Dissapearing the axis
			plt.setp(ax.get_yticklabels(), visible=False)
		
		# (Re)setting labels and titles
		
		for name in ['rama','R']:#['rama','R','theta','d']:
			ax = axes[column][name]
			
			for spine in ax.spines.keys():
				ax.spines[spine].set_linewidth(linewidth*1.2)
			
			addto_title = r'('+next(panel_letters)+r')'+"\n"
			t = ax.set_title(addto_title, size=textsize*1.4, loc='left',horizontalalignment='right')
			#t.set_y(1.03)  # Shifting the title up a little
			t.set_x(-0.12)  # Shifting the title to the left
		
			# axis tick labels
			plt.setp(ax.get_yticklabels(), size=textsize*1)#, rotation="vertical")
			plt.setp(ax.get_xticklabels(), size=textsize*1)
			# axis label (name of the axis)
			ax.set_xlabel(ax.get_xlabel(), fontsize=textsize*2)
			ax.set_ylabel(ax.get_ylabel(), fontsize=textsize*2, rotation="horizontal", ha='center', va='center',labelpad=20)
		
		# Annotating distributions
		for name in ['R']:#axis_to_average_x.keys():
			ax = axes[column][name]
			#
			ss_to_avex  = {}
			ss_to_color = {}
			for ss,avex,color in axis_to_average_x[name]:
				
				ss_label = str(ss)
				if ss_label in ss_name_to_label:
					ss_label = ss_name_to_label[ss_label]
				#avex = round(avex,1)
				ss_to_avex[ss]         =       avex
				ss_to_color[ss]        =       color
			#
			ss_names = ss_to_avex.keys()
			ss_avex  = ss_to_avex.values()
			#
			# sort names by avex:
			ss_names_sorted = [x for (y,x) in sorted(zip(ss_avex,ss_names))]
			#
			min_avex = np.min(ss_avex)
			max_avex = np.max(ss_avex)
			
			unique_avexes = sorted(set(ss_avex))
			#
			paddings = itertools.cycle([0.0,0.0]) # If you want to cycle through some of vertical heights of the labels
			for ss in ss_names_sorted:
				# 0 for first sorted ss, 1 for last
				ss_frac = float(ss_names_sorted.index(ss))/float(len(ss_names_sorted)-1)
				#
				text_x = min_avex + ss_frac*(max_avex-min_avex)
				text_y = 1.9+next(paddings)
				#
				point_x = ss_to_avex[ss]
				point_y = 1.05
				#
				total_height = 15
				# length of vertical leg emanating from the label 
				armA = 14
				# length of vertical leg emanating from the pointed region
				label_frac = float(unique_avexes.index(point_x))/(len(unique_avexes)-1)
				armB = total_height*label_frac
				#
				current_arrow_color = ss_to_color[ss]
				current_text_color  = ss_to_color[ss]
				
				angleB = 90			
				if 0: # Also offset the angle 
					arange = 45
					offset_weight = -1*(1-abs(2.0*(label_frac-0.5)))*np.sign(label_frac-0.5)
					'''
					offset_weight  0  .2  .4  .6  .8   1 -.8 -.6 -.4 -.2   0
						       |   |   |   |   |   |   |   |   |   |   |
					label_frac     0  .1  .2  .3  .4  .5  .6  .7  .8  .9   1
					'''
					angleB_offset = offset_weight*arange # when label_frac == 0 or 1 armB_offset is +arange or -arange.
					angleB        = 90 + angleB_offset
				#
				connection_style = "arc,angleA=-90,angleB=%f,armA=%f,armB=%f,rad=0"%(angleB,armA,armB)
				
				
				weight = 'normal'
				if ss in mark_important_secondary_structures:
					weight = 'bold'
				
				ax.annotate(str(ss), xy=(point_x, point_y), xytext=(text_x, text_y),
					    xycoords='data', weight=weight,
					    arrowprops=dict(arrowstyle="-",color=current_arrow_color, 
					                    connectionstyle=connection_style, lw=linewidth),
					    horizontalalignment='center', verticalalignment='bottom',
					    fontsize=textsize*0.9, color=ss_to_color[ss])
		# Add some artificial limits:
		#axes[column]['theta'].set_xticks(range(-360,361,90))
		#axes[column]['theta'].set_xlim(50,275)
		
						
	plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
	#gs.update(wspace=1, hspace=1)
	#plt.subplots_adjust(wspace=.201)
	#plt.subplots_adjust(hspace=0.001)
	plt.savefig(figname, dpi=180, bbox_inches='tight',transparent=True) #facecolor='w', edgecolor='w',
	if show_graphs:
        	os.system(pdf_viewer+" "+figname)

# 
if 0: # FIGURE_START, shouldbeone
	figname = "manuscript/automated_figures/fig_ss_compact_form.pdf"
	sns.reset_orig()
	
	#sns.set_context("notebook")#, font_scale=1)
	sns.set_style("ticks")
	set_grays()
	
	textsize = 12.0
	linewidth = 1.0
	alpha     = 0.5
	annotation_color = 'teal'
	# VERY IMPORTANT FOR GETTING THE SHAPE OF EACH PANEL CORRECT
	plt.figure(figsize=(13,10))
	
	# Plot distributions of secondary structures in various formats
	sstype = 'segno' # check <ss_codes>
	
	# ----------------
	#                      nrows  ncols
	#                          |  |
	gs = mpl.gridspec.GridSpec(2, 2, 
	                         height_ratios = [5,0.5],#,1,1],
	                         width_ratios  = [1,1])
	axes = {}
	#    Column
	#    |
	axes[0] = {'rama': plt.subplot(gs[0,0]),
	              'R': plt.subplot(gs[1,0])
	          #'theta': plt.subplot(gs[2,0]),
	          #    'd': plt.subplot(gs[3,0])
	          }
	axes[1] = {'rama': plt.subplot(gs[0,1]),
	              'R': plt.subplot(gs[1,1])
	          #'theta': plt.subplot(gs[2,1]),
	          #    'd': plt.subplot(gs[3,1])
	          }
	
	draw_d_contours = 1
	if draw_d_contours:
		# we will draw a d=0 contour line in each plot that has the name 'rama'
		d0df = pd.read_csv(length_dependent_csv_file)
		d0df = d0df[(d0df['L']==8.0) & (d0df['omega']==180.0)]
	
	mark_important_secondary_structures = [4,8,6,11]		
	
	panel_letters = itertools.cycle(list(string.ascii_lowercase))
	for column in [0,1]:
		
		# -----------------
		# Loading up the SS database
		if column == 0:
			normal_levels = [0.3333,0.6666] # For ramachandran plots, given most ss keys
			levels_for_all = [0.92]         # Same, but in case the 'all' key is used for ss
			
			# Loading:
			if not os.path.isfile(master_ss_file):
				prepare_master_ss_file()
			df = pd.read_csv(master_ss_file)
			# Cleaning:
			df = df[(df['d'] != 999)]
			ss_to_study = ['all','G','E','P','H']
			
		
		if column == 1:
			normal_levels = [0.6666,0.9] # For ramachandran plots, given most ss keys
			if not os.path.isfile(fake_ss_file):
				write_fake_secondary_structure_database(sstype = 'segno',omega=180.0,sample_points=50000, sigma=5)
			df = pd.read_csv(fake_ss_file)
			# Lets study ALL fake secondary structures
			ss_to_study = sorted(set(list(df[sstype])))
			
		'''
		ds_     = np.abs(df[(df[sstype] == ss_codes[sstype]['E'])]['d'])
		thetas_ = df['theta']
		Rs_     = df['R']
		mind    ,maxd     =     ds_.min() ,     ds_.max()
		minR    ,maxR     =     Rs_.min() ,     Rs_.max()
		minTheta,maxTheta = thetas_.min() , thetas_.max()
		yextents_dict = {'d':[1,maxd],'R':[0.3,0.67],'R2':[0.3,0.67],'R2':[0.3,0.67],'theta':[minTheta,maxTheta]}
		'''
		
		# Creating a palette through which we cycle for various secondary structures	
		#palette_colors = sns.dark_palette("red", len(ss_to_study))
		palette_colors = sns.color_palette("cubehelix", len(ss_to_study))
		#palette_colors = sns.color_palette("RdBu", len(ss_to_study))
		#if column == 0:
		#palette_colors = sns.color_palette("colorblind",len(ss_to_study))
		palette_colors = sns.hls_palette(len(ss_to_study), l=.4, s=.7)
		
		palette = itertools.cycle(palette_colors)
		
		if draw_d_contours:
			# draw d=0 contour lines
			X = d0df['phi']
			Y = d0df['psi']
			Z = d0df['d']
			
			X,Y,Z    = locallib.xyz_to_ndarrays(X,Y,Z)
			Z = scipy.ndimage.filters.gaussian_filter(Z, 2)#, order=0)
			CS = axes[column]['rama'].contour(X, Y, Z, [0], colors='gray', linewidths=linewidth, linestyles='dashdot')
			if 0: # Draw contour labels
				# Recast levels
				fmt = {}
				for l in CS.levels:
					s = "%1.3f" %(l)
					if float(l) == float(int(l)):
						s = "d=%d" %(l)
					fmt[l] = s
				#
				axes[column]['rama'].clabel(CS, fontsize=textsize*0.9, inline=1, fmt=fmt) #fmt = contour_label_format)
			#
		#
		axis_to_average_x = {}; axis_to_average_x['theta']=[]; axis_to_average_x['d']=[]; axis_to_average_x['R']=[];
		for ss in ss_to_study:
			levels    = normal_levels
			if ss == 'all':
				levels = levels_for_all
			
			xtype = 'phi'; ytype = 'psi';
			
			cdf = df.copy()
			fill_contour_lines = 1
			draw_contour_lines = 0
			smooth             = 2
			cmap = "#D3D3D3"  # This is the default color for an 'all' map
			c    = cmap # Uses only a color if cmap is a string
			if ss != "all":
				# ----------------
				# NEW COLOR SCHEME:
				c1 = [1,1,1];	c2 = next(palette);
				cdict = {#                c1                   c2
					'red':   ((0.00,  c1[0], c1[0]), (1.0, c2[0], c2[0])), 
					'green': ((0.00,  c1[1], c1[1]), (1.0, c2[1], c2[1])),
					'blue':  ((0.00,  c1[2], c1[2]), (1.0, c2[2], c2[2])) }
				current_cmap = LinearSegmentedColormap('tmpcmap', cdict)
				c = current_cmap(0.99999)
				# ----------------
				draw_contour_lines = 1
				smooth             = 2
				cmap=current_cmap
				# Here, if so inclined, you could put a check to transform the phi psi values
				if ss in ss_codes[sstype]:
					cdf = cdf[(cdf[sstype] == ss_codes[sstype][ss])]
				else:
					# Then this may be a 'fake' secondary structure. No translation needed, just use the ss code.
					cdf = cdf[(cdf[sstype] == ss)]
				
			all_size     = len(df)
			current_size = len(cdf)
			print ss,'\t',current_size,'\t',100.0*float(current_size)/all_size
				
			
			# Plotting contours on the ramachandran plot
			# Fill
			locallib.plot_ramachandran_histogram(cdf['phi'],cdf['psi'],levels=levels,cmap=cmap,
				                    fill=1,lines=0, ax=axes[column]['rama'], 
				                    smooth=smooth, linestyles=['solid', 'dotted'],linewidths=0)
			if ss != 'all':
				# Lines
				line_width_modifier = 1.3
				if ss in mark_important_secondary_structures:
					line_width_modifier = 2.3
				#
				locallib.plot_ramachandran_histogram(cdf['phi'],cdf['psi'],levels=[sorted(levels)[-1]],cmap='k',#cmap,
						            fill=0,lines=1, ax=axes[column]['rama'], 
						            smooth=smooth, linestyles=['solid'], linewidths=linewidth*line_width_modifier)
				#
				# Annotate 
				xave = np.median(list(cdf['phi']))
				yave = np.median(list(cdf['psi']))
				ss_label = str(ss)
				weight = 'normal'
				if ss in mark_important_secondary_structures:
					weight = 'bold'
				if ss_label in ss_name_to_label:
					ss_label = ss_name_to_label[ss_label]
				
				axes[column]['rama'].annotate(ss_label, xy=(xave, yave), xytext=(xave+30, yave-30), weight=weight, color=c,
				            xycoords='data', arrowprops=dict(arrowstyle="-",color=c, lw=linewidth), fontsize=textsize*0.9)
				        
			def get_ave(vals):
				X,Y = get_his(vals,norm=0); X=list(X); Y=list(Y)
				return X[Y.index(np.max(Y))]
			
			
			# We study the ss distributions as a function of R, theta, d
			if ss != "all":
				bins = 150
				if column == 1:
					bins = 50
				
				# r
				vals = cdf['R']
				X,Y = get_his(vals,bins=bins)
				axes[column]['R'].fill_between(X,Y,0,facecolor=c,alpha=alpha)
				axes[column]['R'].plot(X,Y,c='k',linewidth=linewidth)
				axis_to_average_x['R'].append([ss,get_ave(vals),c])
				

				'''
				# theta
				vals = do_theta(cdf['theta'])
				#vals = abs(cdf['d'])+(cdf['phi']+cdf['psi'] + 360.)/(720.)
				if 0: # rescale values wrt min and max vals
					mind = min(vals)
					#maxd = max(vals)
					#for i in range(len(vals)):
					#	if vals[i] < 0:
					#	vals[i] = vals + mind 
					vals = vals - mind
				
				X, Y = get_his(vals,bins=bins)
				axes[column]['theta'].fill_between(X,Y,0,facecolor=c,alpha=alpha)
				axes[column]['theta'].plot(X,Y,c='k',linewidth=linewidth)
				axis_to_average_x['theta'].append([ss,get_ave(vals),c])
				
				# d
				tdf = cdf.copy()
				vals = do_d(tdf['d'])
				if 0: # rescale values wrt min and max vals
					mind = min(vals)
					maxd = max(vals)
					tdf.ix[tdf.d<=0, 'd'] = maxd+(tdf.ix[tdf.d <= 0, 'd']-mind)
					vals = tdf['d']*-1*np.sign(np.cos(np.radians(tdf['theta'])/2.0))
				
				#vals = tdf['d']*-1.0*np.cos(np.radians(tdf['theta'])/2.0)
				#vals = np.sign(tdf['d'])*np.sin(np.radians(tdf['theta'])) # this is the h number from Mannige, PeerJ, 2017
				
				X,Y = get_his(vals,bins=bins)
				axes[column]['d'].fill_between(X,Y,0,facecolor=c,alpha=alpha) 
				axes[column]['d'].plot(X,Y,c='k',linewidth=linewidth)
				axis_to_average_x['d'].append([ss,get_ave(vals),c])
			
				# Example of how to plot vertically
				#axes[column]['d'].fill_betweenx(X,Y,0,facecolor=c,alpha=alpha) 
				#axes[column]['d'].plot(Y,X,c='k',linewidth=linewidth)
				'''

		# setting the ramachandran plot aspect ratio
		if 1:#for name in ['d','theta','R']: #axes.items():	
			axes[column]['rama'].set_aspect(1)
		
		# Making some tweaks to all but the 
		for name in ['R']: #['d','theta','R']: #axes.items():	
			ax = axes[column][name]
		
			xmin,xmax = ax.get_xlim(); ymin,ymax = ax.get_ylim()  # getting axes min max values
			ax.set_yticks([0,1.5])
			
			ax.set_xlabel(type_to_label[name])
		
			# Only show ticks on the left and bottom spines
			ax.xaxis.set_ticks_position('bottom')
			ax.yaxis.set_ticks_position('none')
			# Hide the right and top spines
			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)
			ax.spines['left'].set_visible(False)
		
			xpad =  (xmax-xmin)*0.0; ypad =  (ymax-ymin)*0.0;    # setting padding
		
			extent = [xmin, xmax+xpad, ymin, ymax+ypad+0.6] # storing (padded) axis min max values
			xscaling = .094                                          # setting how much larger must the x axis be compared to the y axis
			ax.axis(extent)                                       # setting (padded) axis min max values
			#ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/xscaling)
			
			# Dissapearing the axis
			plt.setp(ax.get_yticklabels(), visible=False)
		
		# (Re)setting labels and titles
		
		for name in ['rama','R']:#['rama','R','theta','d']:
			ax = axes[column][name]
			
			for spine in ax.spines.keys():
				ax.spines[spine].set_linewidth(linewidth*1.2)
			
			addto_title = r'('+next(panel_letters)+r')'+"\n"
			t = ax.set_title(addto_title, size=textsize*1.4, loc='left',horizontalalignment='right')
			#t.set_y(1.03) # Shifting the title up a little
		
			# axis tick labels
			plt.setp(ax.get_yticklabels(), size=textsize*1)#, rotation="vertical")
			plt.setp(ax.get_xticklabels(), size=textsize*1)
			# axis label (name of the axis)
			ax.set_xlabel(ax.get_xlabel(), fontsize=textsize*2)
			ax.set_ylabel(ax.get_ylabel(), fontsize=textsize*2, rotation="horizontal", ha='center', va='center',labelpad=20)
		
		# Annotating distributions
		for name in ['R']:#axis_to_average_x.keys():
			ax = axes[column][name]
			#
			ss_to_avex  = {}
			ss_to_color = {}
			for ss,avex,color in axis_to_average_x[name]:
				
				ss_label = str(ss)
				if ss_label in ss_name_to_label:
					ss_label = ss_name_to_label[ss_label]
				#avex = round(avex,1)
				ss_to_avex[ss]         =       avex
				ss_to_color[ss]        =       color
			#
			ss_names = ss_to_avex.keys()
			ss_avex  = ss_to_avex.values()
			#
			# sort names by avex:
			ss_names_sorted = [x for (y,x) in sorted(zip(ss_avex,ss_names))]
			#
			min_avex = np.min(ss_avex)
			max_avex = np.max(ss_avex)
			
			unique_avexes = sorted(set(ss_avex))
			#
			paddings = itertools.cycle([0.0,0.0]) # If you want to cycle through some of vertical heights of the labels
			for ss in ss_names_sorted:
				# 0 for first sorted ss, 1 for last
				ss_frac = float(ss_names_sorted.index(ss))/float(len(ss_names_sorted)-1)
				#
				text_x = min_avex + ss_frac*(max_avex-min_avex)
				text_y = 1.9+next(paddings)
				#
				point_x = ss_to_avex[ss]
				point_y = 1.05
				#
				total_height = 15
				# length of vertical leg emanating from the label 
				armA = 14
				# length of vertical leg emanating from the pointed region
				label_frac = float(unique_avexes.index(point_x))/(len(unique_avexes)-1)
				armB = total_height*label_frac
				#
				current_arrow_color = ss_to_color[ss]
				current_text_color  = ss_to_color[ss]
				
				angleB = 90			
				if 0: # Also offset the angle 
					arange = 45
					offset_weight = -1*(1-abs(2.0*(label_frac-0.5)))*np.sign(label_frac-0.5)
					'''
					offset_weight  0  .2  .4  .6  .8   1 -.8 -.6 -.4 -.2   0
						       |   |   |   |   |   |   |   |   |   |   |
					label_frac     0  .1  .2  .3  .4  .5  .6  .7  .8  .9   1
					'''
					angleB_offset = offset_weight*arange # when label_frac == 0 or 1 armB_offset is +arange or -arange.
					angleB        = 90 + angleB_offset
				#
				connection_style = "arc,angleA=-90,angleB=%f,armA=%f,armB=%f,rad=0"%(angleB,armA,armB)
				
				
				weight = 'normal'
				if ss in mark_important_secondary_structures:
					weight = 'bold'
				
				ax.annotate(str(ss), xy=(point_x, point_y), xytext=(text_x, text_y),
					    xycoords='data', weight=weight,
					    arrowprops=dict(arrowstyle="-",color=current_arrow_color, 
					                    connectionstyle=connection_style, lw=linewidth),
					    horizontalalignment='center', verticalalignment='bottom',
					    fontsize=textsize*0.9, color=ss_to_color[ss])
		# Add some artificial limits:
		#axes[column]['theta'].set_xticks(range(-360,361,90))
		#axes[column]['theta'].set_xlim(50,275)
		
						
	plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
	#gs.update(wspace=1, hspace=1)
	#plt.subplots_adjust(wspace=.201)
	#plt.subplots_adjust(hspace=0.001)
	plt.savefig(figname, dpi=180, bbox_inches='tight',transparent=True) #facecolor='w', edgecolor='w',
	if show_graphs:
        	os.system(pdf_viewer+" "+figname)
	

# 
if 0: # FIGURE_START, shouldbeone
	figname = "manuscript/automated_figures/fig_various_ss_formats.pdf"
	sns.reset_orig()
	
	#sns.set_context("notebook")#, font_scale=1)
	sns.set_style("ticks")
	set_grays()
	
	textsize = 12.0
	linewidth = 1.0
	alpha     = 0.5
	annotation_color = 'teal'
	# VERY IMPORTANT FOR GETTING THE SHAPE OF EACH PANEL CORRECT
	plt.figure(figsize=(13,12))
	
	# Plot distributions of secondary structures in various formats
	sstype = 'segno' # check <ss_codes>
	
	# ----------------
	#                      nrows  ncols
	#                          |  |
	gs = mpl.gridspec.GridSpec(4, 2, 
	                         height_ratios = [5,1,1,1],
	                         width_ratios  = [1,1])
	axes = {}
	#    Column
	#    |
	axes[0] = {'rama': plt.subplot(gs[0,0]),
	              'R': plt.subplot(gs[1,0]),
	          'theta': plt.subplot(gs[2,0]),
	              'd': plt.subplot(gs[3,0])
	          }
	axes[1] = {'rama': plt.subplot(gs[0,1]),
	              'R': plt.subplot(gs[1,1]),
	          'theta': plt.subplot(gs[2,1]),
	              'd': plt.subplot(gs[3,1])
	          }
	
	draw_d_contours = 1
	if draw_d_contours:
		# we will draw a d=0 contour line in each plot that has the name 'rama'
		d0df = pd.read_csv(length_dependent_csv_file)
		d0df = d0df[(d0df['L']==8.0) & (d0df['omega']==180.0)]
	
	mark_important_secondary_structures = [4,8,6,11]		
	
	panel_letters = itertools.cycle(list(string.ascii_lowercase))
	for column in [0,1]:
		
		# -----------------
		# Loading up the SS database
		if column == 0:
			normal_levels = [0.3333,0.6666] # For ramachandran plots, given most ss keys
			levels_for_all = [0.92]         # Same, but in case the 'all' key is used for ss
			
			# Loading:
			if not os.path.isfile(master_ss_file):
				prepare_master_ss_file()
			df = pd.read_csv(master_ss_file)
			# Cleaning:
			df = df[(df['d'] != 999)]
			ss_to_study = ['all','G','E','P','H']
			
		
		if column == 1:
			normal_levels = [0.6666,0.9] # For ramachandran plots, given most ss keys
			if not os.path.isfile(fake_ss_file):
				write_fake_secondary_structure_database(sstype = 'segno',omega=180.0,sample_points=50000, sigma=5)
			df = pd.read_csv(fake_ss_file)
			# Lets study ALL fake secondary structures
			ss_to_study = sorted(set(list(df[sstype])))
			
		'''
		ds_     = np.abs(df[(df[sstype] == ss_codes[sstype]['E'])]['d'])
		thetas_ = df['theta']
		Rs_     = df['R']
		mind    ,maxd     =     ds_.min() ,     ds_.max()
		minR    ,maxR     =     Rs_.min() ,     Rs_.max()
		minTheta,maxTheta = thetas_.min() , thetas_.max()
		yextents_dict = {'d':[1,maxd],'R':[0.3,0.67],'R2':[0.3,0.67],'R2':[0.3,0.67],'theta':[minTheta,maxTheta]}
		'''
		
		# Creating a palette through which we cycle for various secondary structures	
		#palette_colors = sns.dark_palette("red", len(ss_to_study))
		palette_colors = sns.color_palette("cubehelix", len(ss_to_study))
		#palette_colors = sns.color_palette("RdBu", len(ss_to_study))
		#if column == 0:
		#palette_colors = sns.color_palette("colorblind",len(ss_to_study))
		palette_colors = sns.hls_palette(len(ss_to_study), l=.4, s=.7)
		
		palette = itertools.cycle(palette_colors)
		
		if draw_d_contours:
			# draw d=0 contour lines
			X = d0df['phi']
			Y = d0df['psi']
			Z = d0df['d']
			
			X,Y,Z    = locallib.xyz_to_ndarrays(X,Y,Z)
			Z = scipy.ndimage.filters.gaussian_filter(Z, 2)#, order=0)
			CS = axes[column]['rama'].contour(X, Y, Z, [0], colors='gray', linewidths=linewidth, linestyles='dashdot')
			if 0: # Draw contour labels
				# Recast levels
				fmt = {}
				for l in CS.levels:
					s = "%1.3f" %(l)
					if float(l) == float(int(l)):
						s = "d=%d" %(l)
					fmt[l] = s
				#
				axes[column]['rama'].clabel(CS, fontsize=textsize*0.9, inline=1, fmt=fmt) #fmt = contour_label_format)
			#
		#
		axis_to_average_x = {}; axis_to_average_x['theta']=[]; axis_to_average_x['d']=[]; axis_to_average_x['R']=[];
		for ss in ss_to_study:
			levels    = normal_levels
			if ss == 'all':
				levels = levels_for_all
			
			xtype = 'phi'; ytype = 'psi';
			
			cdf = df.copy()
			fill_contour_lines = 1
			draw_contour_lines = 0
			smooth             = 2
			cmap = "#D3D3D3"  # This is the default color for an 'all' map
			c    = cmap # Uses only a color if cmap is a string
			if ss != "all":
				# ----------------
				# NEW COLOR SCHEME:
				c1 = [1,1,1];	c2 = next(palette);
				cdict = {#                c1                   c2
					'red':   ((0.00,  c1[0], c1[0]), (1.0, c2[0], c2[0])), 
					'green': ((0.00,  c1[1], c1[1]), (1.0, c2[1], c2[1])),
					'blue':  ((0.00,  c1[2], c1[2]), (1.0, c2[2], c2[2])) }
				current_cmap = LinearSegmentedColormap('tmpcmap', cdict)
				c = current_cmap(0.99999)
				# ----------------
				draw_contour_lines = 1
				smooth             = 2
				cmap=current_cmap
				# Here, if so inclined, you could put a check to transform the phi psi values
				if ss in ss_codes[sstype]:
					cdf = cdf[(cdf[sstype] == ss_codes[sstype][ss])]
				else:
					# Then this may be a 'fake' secondary structure. No translation needed, just use the ss code.
					cdf = cdf[(cdf[sstype] == ss)]
				
			all_size     = len(df)
			current_size = len(cdf)
			print ss,'\t',current_size,'\t',100.0*float(current_size)/all_size
				
			
			# Plotting contours on the ramachandran plot
			# Fill
			locallib.plot_ramachandran_histogram(cdf['phi'],cdf['psi'],levels=levels,cmap=cmap,
				                    fill=1,lines=0, ax=axes[column]['rama'], 
				                    smooth=smooth, linestyles=['solid', 'dotted'],linewidths=0)
			if ss != 'all':
				# Lines
				line_width_modifier = 1.3
				if ss in mark_important_secondary_structures:
					line_width_modifier = 2.3
				#
				locallib.plot_ramachandran_histogram(cdf['phi'],cdf['psi'],levels=[sorted(levels)[-1]],cmap='k',#cmap,
						            fill=0,lines=1, ax=axes[column]['rama'], 
						            smooth=smooth, linestyles=['solid'], linewidths=linewidth*line_width_modifier)
				#
				# Annotate 
				xave = np.median(list(cdf['phi']))
				yave = np.median(list(cdf['psi']))
				ss_label = str(ss)
				weight = 'normal'
				if ss in mark_important_secondary_structures:
					weight = 'bold'
				if ss_label in ss_name_to_label:
					ss_label = ss_name_to_label[ss_label]
				
				axes[column]['rama'].annotate(ss_label, xy=(xave, yave), xytext=(xave+30, yave-30), weight=weight, color=c,
				            xycoords='data', arrowprops=dict(arrowstyle="-",color=c, lw=linewidth), fontsize=textsize*0.9)
				        
			def get_ave(vals):
				X,Y = get_his(vals,norm=0); X=list(X); Y=list(Y)
				return X[Y.index(np.max(Y))]
			
			
			# We study the ss distributions as a function of R, theta, d
			if ss != "all":
				bins = 150
				if column == 1:
					bins = 50
				
				# r
				vals = cdf['R']
				X,Y = get_his(vals,bins=bins)
				axes[column]['R'].fill_between(X,Y,0,facecolor=c,alpha=alpha)
				axes[column]['R'].plot(X,Y,c='k',linewidth=linewidth)
				axis_to_average_x['R'].append([ss,get_ave(vals),c])
				
				# theta
				vals = do_theta(cdf['theta'])
				#vals = abs(cdf['d'])+(cdf['phi']+cdf['psi'] + 360.)/(720.)
				if 0: # rescale values wrt min and max vals
					mind = min(vals)
					#maxd = max(vals)
					#for i in range(len(vals)):
					#	if vals[i] < 0:
					#	vals[i] = vals + mind 
					vals = vals - mind
				
				X, Y = get_his(vals,bins=bins)
				axes[column]['theta'].fill_between(X,Y,0,facecolor=c,alpha=alpha)
				axes[column]['theta'].plot(X,Y,c='k',linewidth=linewidth)
				axis_to_average_x['theta'].append([ss,get_ave(vals),c])
				
				# d
				tdf = cdf.copy()
				vals = do_d(tdf['d'])
				if 0: # rescale values wrt min and max vals
					mind = min(vals)
					maxd = max(vals)
					tdf.ix[tdf.d<=0, 'd'] = maxd+(tdf.ix[tdf.d <= 0, 'd']-mind)
					vals = tdf['d']*-1*np.sign(np.cos(np.radians(tdf['theta'])/2.0))
				
				#vals = tdf['d']*-1.0*np.cos(np.radians(tdf['theta'])/2.0)
				#vals = np.sign(tdf['d'])*np.sin(np.radians(tdf['theta'])) # this is the h number from Mannige, PeerJ, 2017
				
				X,Y = get_his(vals,bins=bins)
				axes[column]['d'].fill_between(X,Y,0,facecolor=c,alpha=alpha) 
				axes[column]['d'].plot(X,Y,c='k',linewidth=linewidth)
				axis_to_average_x['d'].append([ss,get_ave(vals),c])
			
				# Example of how to plot vertically
				#axes[column]['d'].fill_betweenx(X,Y,0,facecolor=c,alpha=alpha) 
				#axes[column]['d'].plot(Y,X,c='k',linewidth=linewidth)
				
		# setting the ramachandran plot aspect ratio
		if 1:#for name in ['d','theta','R']: #axes.items():	
			axes[column]['rama'].set_aspect(1)
		
		# Making some tweaks to all but the 
		for name in ['d','theta','R']: #axes.items():	
			ax = axes[column][name]
		
			xmin,xmax = ax.get_xlim(); ymin,ymax = ax.get_ylim()  # getting axes min max values
			ax.set_yticks([0,1.5])
			
			ax.set_xlabel(type_to_label[name])
		
			# Only show ticks on the left and bottom spines
			ax.xaxis.set_ticks_position('bottom')
			ax.yaxis.set_ticks_position('none')
			# Hide the right and top spines
			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)
			ax.spines['left'].set_visible(False)
		
			xpad =  (xmax-xmin)*0.0; ypad =  (ymax-ymin)*0.0;    # setting padding
		
			extent = [xmin, xmax+xpad, ymin, ymax+ypad+0.6] # storing (padded) axis min max values
			xscaling = .094                                          # setting how much larger must the x axis be compared to the y axis
			ax.axis(extent)                                       # setting (padded) axis min max values
			#ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/xscaling)
			
			# Dissapearing the axis
			plt.setp(ax.get_yticklabels(), visible=False)
		
		# (Re)setting labels and titles
		for name in ['rama','R','theta','d']:
			ax = axes[column][name]
			
			for spine in ax.spines.keys():
				ax.spines[spine].set_linewidth(linewidth*1.2)
			
			addto_title = r'('+next(panel_letters)+r')'+"\n"
			t = ax.set_title(addto_title, size=textsize*1.4, loc='left',horizontalalignment='right')
			#t.set_y(1.03) # Shifting the title up a little
		
			# axis tick labels
			plt.setp(ax.get_yticklabels(), size=textsize*1)#, rotation="vertical")
			plt.setp(ax.get_xticklabels(), size=textsize*1)
			# axis label (name of the axis)
			ax.set_xlabel(ax.get_xlabel(), fontsize=textsize*2)
			ax.set_ylabel(ax.get_ylabel(), fontsize=textsize*2, rotation="horizontal", ha='center', va='center',labelpad=20)
		
		# Annotating distributions
		for name in axis_to_average_x.keys():
			ax = axes[column][name]
			#
			ss_to_avex  = {}
			ss_to_color = {}
			for ss,avex,color in axis_to_average_x[name]:
				
				ss_label = str(ss)
				if ss_label in ss_name_to_label:
					ss_label = ss_name_to_label[ss_label]
				#avex = round(avex,1)
				ss_to_avex[ss]         =       avex
				ss_to_color[ss]        =       color
			#
			ss_names = ss_to_avex.keys()
			ss_avex  = ss_to_avex.values()
			#
			# sort names by avex:
			ss_names_sorted = [x for (y,x) in sorted(zip(ss_avex,ss_names))]
			#
			min_avex = np.min(ss_avex)
			max_avex = np.max(ss_avex)
			
			unique_avexes = sorted(set(ss_avex))
			#
			paddings = itertools.cycle([0.0,0.0]) # If you want to cycle through some of vertical heights of the labels
			for ss in ss_names_sorted:
				# 0 for first sorted ss, 1 for last
				ss_frac = float(ss_names_sorted.index(ss))/float(len(ss_names_sorted)-1)
				#
				text_x = min_avex + ss_frac*(max_avex-min_avex)
				text_y = 1.9+next(paddings)
				#
				point_x = ss_to_avex[ss]
				point_y = 1.05
				#
				total_height = 15
				# length of vertical leg emanating from the label 
				armA = 14
				# length of vertical leg emanating from the pointed region
				label_frac = float(unique_avexes.index(point_x))/(len(unique_avexes)-1)
				armB = total_height*label_frac
				#
				current_arrow_color = ss_to_color[ss]
				current_text_color  = ss_to_color[ss]
				
				angleB = 90			
				if 0: # Also offset the angle 
					arange = 45
					offset_weight = -1*(1-abs(2.0*(label_frac-0.5)))*np.sign(label_frac-0.5)
					'''
					offset_weight  0  .2  .4  .6  .8   1 -.8 -.6 -.4 -.2   0
						       |   |   |   |   |   |   |   |   |   |   |
					label_frac     0  .1  .2  .3  .4  .5  .6  .7  .8  .9   1
					'''
					angleB_offset = offset_weight*arange # when label_frac == 0 or 1 armB_offset is +arange or -arange.
					angleB        = 90 + angleB_offset
				#
				connection_style = "arc,angleA=-90,angleB=%f,armA=%f,armB=%f,rad=0"%(angleB,armA,armB)
				
				
				weight = 'normal'
				if ss in mark_important_secondary_structures:
					weight = 'bold'
				
				ax.annotate(str(ss), xy=(point_x, point_y), xytext=(text_x, text_y),
					    xycoords='data', weight=weight,
					    arrowprops=dict(arrowstyle="-",color=current_arrow_color, 
					                    connectionstyle=connection_style, lw=linewidth),
					    horizontalalignment='center', verticalalignment='bottom',
					    fontsize=textsize*0.9, color=ss_to_color[ss])
		# Add some artificial limits:
		#axes[column]['theta'].set_xticks(range(-360,361,90))
		#axes[column]['theta'].set_xlim(50,275)
		
						
	plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
	#gs.update(wspace=1, hspace=1)
	#plt.subplots_adjust(wspace=.201)
	#plt.subplots_adjust(hspace=0.001)
	plt.savefig(figname, dpi=180, bbox_inches='tight',transparent=True) #facecolor='w', edgecolor='w',
	if show_graphs:
        	os.system(pdf_viewer+" "+figname)
	
