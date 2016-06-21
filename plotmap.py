helpme = """
============================================
        _       _     __  __    _    ____  
  _ __ | | ___ | |_  |  \/  |  / \  |  _ \ 
 | '_ \| |/ _ \| __| | |\/| | / _ \ | |_) |
 | |_) | | (_) | |_  | |  | |/ ___ \|  __/ 
 | .__/|_|\___/ \__| |_|  |_/_/   \_\_|      v 0.0.0.0.0.0.0...
 |_|                 (Multi-angle Picture)

This tool provides easily readable "pictures" of protein conformations, 
ensembles, and trajectories saved as either a combined protein databank 
(PDB) structure file, or a directory of such files, and produces graphs.
-----
Usage
-----
python plotmap.py -pdb ProteinDatabankStructureFilename.pdb
python plotmap.py -pdb /directory/containing/pdbs/
------
Output (the x-axis always represents the models/structures listed in the PDB)
------
filename.rcode.eps      (y-axis: residue #; color: R number based on "-signed" and <rcode_cmap>)
filename.rcode.his.eps  (y-axis: Ramachandran number (R); color: frequency of R in model)
filename.rcode.rmsf.eps (y-axis: residue #; color: RMSF in R from the previous model)
---------------
Additional tags
---------------
-h       -     Prints this message
-ss      -     Color the ramachandran number codes (R-codes) by 
               secondary structure (default: color by chirality and sign)
-signed  -     Use the signed version of the ramachandran number
-rmsd    -     Also producee "filename.rcode.rmsd.eps"
               (y-axis: residue #; color: RMSD in R from first model)
---------------
Each graph is also accompanied by "_colorbar.eps", which are keys.
---------------
The Ramachandran number concept is discussed in the manuscript:
Mannige, Kundu, Whitelam (2016) "The Ramachandran number: an order parameter for protein geometry" 
Preprint at: http://arxiv.org/abs/1511.03011
============================================
"""

#print helpme
#exit()
import sys,string

signed = 0
rrange = [0,1]
colortype = "Chirality" # can be SecondaryStructure

showeps = 0
dofilter = 0

showrcode = 1
showhis   = 1
showrmsf  = 1
showrmsd  = 0
do_vmd_etc = 1

bins = 100
pdbfn = ""

# python plotmap.py -pdb /home/ranjan/Desktop/old/pairing_functions/for_sharing/structures/nanosheet_birth_U7.pdb
# python plotmap.py -pdb /home/ranjan/Desktop/old/pairing_functions/for_sharing/structures/nanosheet_traj.pdb
# python plotmap.py -pdb /home/ranjan/Desktop/old/pairing_functions/for_sharing/structures/class_a_alpha_1MBA.pdb
# python plotmap.py -pdb /home/ranjan/Desktop/old/pairing_functions/for_sharing/structures/class_c_a_plus_b_2ACY.pdb
import glob
import time
from matplotlib import colors
import copy, re,os
import matplotlib.pyplot as plt 
import matplotlib
import math,numpy
import numpy as np
from scipy import ndimage
#import pylab
import os
import Geometry
import PeptideBuilder
import Bio.PDB
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


forcedmax = False
forcedmin = False

show_graphs = 1
default_fontsize = 22
colorbarXscaling = 0.08
defaultXscaling  = 2.0

SCALE = 10.0 # For the postscript output


# ===================================================================================
# SETTING UP SOME COLORMAPS

COLORSWITCH = 0.45 # THIS IS THE POINT, FOR THE RED/BLUE AND RED/BLUE/YELLOW/BLACK 
                   # COLOR SCHEME, WHERE THE SWITCH IN COLOR HAPPENS (NAIVELY 0.5, 
                   # BUT BETA SHEETS SPILL TO THE "D" PORTION OF THE PLOT, SO IT 
                   # IS 0.45

# First, some definitions:
# DEFINING COLORS BY CHIRALITY:
# c stands for color, bc stands for background color             
             #    when R ranges from   [-1,1]  ("Signed R")     [0,1] (Traditional R)
             #                         ------------             -----------
c1 = [0,0,0] # black                   | \_ c4  / |             |\_    c4 |
c2 = [1,1,0] # yellow                  |   \_ /   |             |  \_     |
c3 = [1,0,0] # red                 psi |c3  /\_c2 |         psi |    \_   |
c4 = [0,0,1] # blue                    |  /    \_ |             |      \_ |
bc = [1,1,1] # white                   |/  c1    \|             |c3      \|
             #                         ------------             -----------
             #                             phi                      phi
# DEFINING POSITIONS AND COLORS BY SECONDARY STRUCTURE:
# POSITIONS
helix_start = 0.31 # the start of the helical region (all assuming R in [0,1])
helix_end   = 0.39 # the end of the helical region
sheet_start = 0.45 # the start of the sheet region
sheet_end   = 0.62 # the end of the sheet region
polyproline_end = 0.66 # the end of the polyprolineII region 
                       # (the start coincides with the sheet region, 
                       # so it just begins after the sheet region ends)
# COLORS
helixR      = (1,0,0)
sheet       = (0,0,1)
polyproline = (0,1,1)

# ----------------
# NEW COLOR SCHEME: color by backbone twist (expected range: R=[0,1])
# ----------------
# This lets you later on get the cmap by name 'TwoColor': cmap = plt.get_cmap('TwoColor')
# POSITION: 0        COLORSWITCH         1
#    COLOR: | white - red | blue - white |
cdict = {
#                         white  white          red    blue          white  white
	'red':   ((0.00,  bc[0], bc[0]), (COLORSWITCH,  c3[0], c4[0]), (1.0, bc[0], bc[0])), 
	'green': ((0.00,  bc[1], bc[1]), (COLORSWITCH,  c3[1], c4[1]), (1.0, bc[1], bc[1])),
	'blue':  ((0.00,  bc[2], bc[2]), (COLORSWITCH,  c3[2], c4[2]), (1.0, bc[2], bc[2])) 
}
cmap = LinearSegmentedColormap('ChiralityTwoColor', cdict)
plt.register_cmap(cmap=cmap)
# ----------------
# NEW COLOR SCHEME: color by backbone twist, variant (expected range: R=[0,1])
# ----------------
cdict = {
#                         white  white          blue   blue
	'red':   ((0.00,  bc[0], bc[0]), (1.0,  c4[0], c4[0])), 
	'green': ((0.00,  bc[1], bc[1]), (1.0,  c4[1], c4[1])),
	'blue':  ((0.00,  bc[2], bc[2]), (1.0,  c4[2], c4[2]))
}
cmap = LinearSegmentedColormap('deleteme', cdict)
plt.register_cmap(cmap=cmap)
cdict = {
#                                       white  white         blue   blue
	'red':   ((0.00,  1, 1), (0.5,  bc[0], bc[0]), (1.0, c4[0], c4[0])), 
	'green': ((0.00,  0, 0	), (0.5,  bc[1], bc[1]), (1.0, c4[1], c4[1])),
	'blue':  ((0.00,  1, 1), (0.5,  bc[2], bc[2]), (1.0, c4[2], c4[2])) 
}
cmap = LinearSegmentedColormap('deletemeSigned', cdict)
plt.register_cmap(cmap=cmap)
# ----------------
# NEW COLOR SCHEME: color by backbone twist, variant (expected range: R=[0,1])
# ----------------
# This lets you later on get the cmap by name 'TwoColorInverted': cmap = plt.get_cmap('TwoColorInverted')
# POSITION: 0              0.25             0.5           0.75            1
#    COLOR: | white - black | yellow - white | white - red | blue - white |
cdict = {
#                         red    red                    white  white         blue   blue
	'red':   ((0.00,  c3[0], c3[0]), (COLORSWITCH,  bc[0], bc[0]), (1.0, c4[0], c4[0])), 
	'green': ((0.00,  c3[1], c3[1]), (COLORSWITCH,  bc[1], bc[1]), (1.0, c4[1], c4[1])),
	'blue':  ((0.00,  c3[2], c3[2]), (COLORSWITCH,  bc[2], bc[2]), (1.0, c4[2], c4[2])) 
}
cmap = LinearSegmentedColormap('ChiralityTwoColor_r', cdict)
plt.register_cmap(cmap=cmap)
# ----------------
# NEW COLOR SCHEME: color by backbone twist (expected range: R=[-1,1])
# ----------------
# This lets you later on get the cmap by name 'FourColor': cmap = plt.get_cmap('FourColor')
# POSITION: 0              0.25             0.5           0.75            1
#    COLOR: | white - black | yellow - white | white - red | blue - white |
cdict = {
#                         white  white           black  yellow         white  white           white  white         blue   blue
	'red':   ((0.00,  bc[0], bc[0]), (0.25,  c1[0], c2[0]), (0.50, bc[0], bc[0]), (0.75,  c3[0], c4[0]), (1.0, bc[0], bc[0])), 
	'green': ((0.00,  bc[1], bc[1]), (0.25,  c1[1], c2[1]), (0.50, bc[1], bc[1]), (0.75,  c3[1], c4[1]), (1.0, bc[1], bc[1])),
	'blue':  ((0.00,  bc[2], bc[2]), (0.25,  c1[2], c2[2]), (0.50, bc[2], bc[2]), (0.75,  c3[2], c4[2]), (1.0, bc[2], bc[2])) 
}
cmap = LinearSegmentedColormap('ChiralityFourColor', cdict)
plt.register_cmap(cmap=cmap)
# ----------------
# NEW COLOR SCHEME: color by backbone twist, variant (expected range: R=[-1,1])
# ----------------
# This lets you later on get the cmap by name 'FourColorInverted': cmap = plt.get_cmap('FourColorInverted')
# POSITION: 0              0.25             0.5           0.75            1
#    COLOR: | black - white | white - yellow | red - white | white - blue |
cdict = {
#                         black  black           white  white         yellow  red             white  white         blue   blue
	'red':   ((0.00,  c1[0], c1[0]), (0.25,  bc[0], bc[0]), (0.50, c2[0], c3[0]), (0.75,  bc[0], bc[0]), (1.0, c4[0], c4[0])), 
	'green': ((0.00,  c1[1], c1[1]), (0.25,  bc[1], bc[1]), (0.50, c2[1], c3[1]), (0.75,  bc[1], bc[1]), (1.0, c4[1], c4[1])),
	'blue':  ((0.00,  c1[2], c1[2]), (0.25,  bc[2], bc[2]), (0.50, c2[2], c3[2]), (0.75,  bc[2], bc[2]), (1.0, c4[2], c4[2])) 
}
cmap = LinearSegmentedColormap('ChiralityFourColor_r', cdict)
plt.register_cmap(cmap=cmap)
# -------------------------
# NEW COLOR SCHEME: secondary structure (expected range: R=[0,1])
# ----------------
# This lets you later on get the cmap by name 'SecondaryStructure': cmap = plt.get_cmap('SecondaryStructure')
# POSITION: 0          helix_start       helix_end       sheet_start     sheet_end             polyproline_end            1
#    COLOR: | white - white | helixR - helixR | white - white | sheet - sheet | polyproline - polyproline | white - white |
cdict = {  
           'red': ((0.00,  bc[0], bc[0]), (helix_start,  bc[0], helixR[0]), (helix_end,  helixR[0], bc[0]), (sheet_start,  bc[0], sheet[0]), (sheet_end,  sheet[0], polyproline[0]), (polyproline_end, polyproline[0], bc[0]), (1, bc[0],bc[0])), 
         'green': ((0.00,  bc[1], bc[1]), (helix_start,  bc[1], helixR[1]), (helix_end,  helixR[1], bc[1]), (sheet_start,  bc[1], sheet[1]), (sheet_end,  sheet[1], polyproline[1]), (polyproline_end, polyproline[1], bc[1]), (1, bc[1],bc[1])),
          'blue': ((0.00,  bc[2], bc[2]), (helix_start,  bc[2], helixR[2]), (helix_end,  helixR[2], bc[2]), (sheet_start,  bc[2], sheet[2]), (sheet_end,  sheet[2], polyproline[2]), (polyproline_end, polyproline[2], bc[2]), (1, bc[2],bc[2]))
        }
cmap = LinearSegmentedColormap('SecondaryStructureTwoColor', cdict)
plt.register_cmap(cmap=cmap)
# ----------------
# NEW COLOR SCHEME: color by secondary structure (expected range: R=[-1,1])
# ----------------
# POSITION (MIRRORRED AROUND 0): 0          helix_start       helix_end       sheet_start     sheet_end             polyproline_end            1
#                         COLOR: | white - white | helixR - helixR | white - white | sheet - sheet | polyproline - polyproline | white - white |
cdict = {  
           'red': [[-1,  bc[0], bc[0]], [polyproline_end*-1, bc[0],polyproline[0]], [sheet_end*-1,  polyproline[0],sheet[0]], [sheet_start*-1,  sheet[0], bc[0]], [helix_end*-1, bc[0],helixR[0]], [helix_start*-1, helixR[0],bc[0]], [helix_start,  bc[0], helixR[0]], [helix_end,  helixR[0], bc[0]], [sheet_start,  bc[0], sheet[0]], [sheet_end,  sheet[0], polyproline[0]], [polyproline_end, polyproline[0], bc[0]], [1, bc[0],bc[0]]],
         'green': [[-1,  bc[1], bc[1]], [polyproline_end*-1, bc[1],polyproline[1]], [sheet_end*-1,  polyproline[1],sheet[1]], [sheet_start*-1,  sheet[1], bc[1]], [helix_end*-1, bc[1],helixR[1]], [helix_start*-1, helixR[1],bc[1]], [helix_start,  bc[1], helixR[1]], [helix_end,  helixR[1], bc[1]], [sheet_start,  bc[1], sheet[1]], [sheet_end,  sheet[1], polyproline[1]], [polyproline_end, polyproline[1], bc[1]], [1, bc[1],bc[1]]], 
          'blue': [[-1,  bc[2], bc[2]], [polyproline_end*-1, bc[2],polyproline[2]], [sheet_end*-1,  polyproline[2],sheet[2]], [sheet_start*-1,  sheet[2], bc[2]], [helix_end*-1, bc[2],helixR[2]], [helix_start*-1, helixR[2],bc[2]], [helix_start,  bc[2], helixR[2]], [helix_end,  helixR[2], bc[2]], [sheet_start,  bc[2], sheet[2]], [sheet_end,  sheet[2], polyproline[2]], [polyproline_end, polyproline[2], bc[2]], [1, bc[2],bc[2]]]  
        } 
# this cdict is not normalized from 0 to 1, which is required for the line following the "for" loop.
minpos = False
maxpos = False
for color in  cdict.keys():
	for i in range(len(cdict[color])):
		if minpos == False:
			minpos = cdict[color][i][0]
		if maxpos == False:
			maxpos = cdict[color][i][0]
		if minpos > cdict[color][i][0]:
			minpos = cdict[color][i][0]
		if maxpos < cdict[color][i][0]:
			maxpos = cdict[color][i][0]
for color in  cdict.keys():
	for i in range(len(cdict[color])):
		cdict[color][i][0] = float(cdict[color][i][0]-minpos)/(maxpos-minpos)
cmap = LinearSegmentedColormap('SecondaryStructureFourColor', cdict)
plt.register_cmap(cmap=cmap)


#rcode_cmap = plt.get_cmap('ChiralityTwoColor')
#rcode_cmap = plt.get_cmap('ChiralityTwoColor_r')
rcode_cmap = plt.get_cmap('ChiralityFourColor')
#rcode_cmap = plt.get_cmap('ChiralityFourColor_r')
#rcode_cmap = plt.get_cmap('SecondaryStructureTwoColor')
#rcode_cmap = plt.get_cmap('SecondaryStructureFourColor')


#secondarystructure_cmap = my_cmap
# ===================================================================================

f = open("templates/postscript.ps","r")
postscript_template = f.read()
f.close()

# ===================================================================================

# Simple smoothing function
def median_filter(vals,nearest_neighbors=1):
	new_vals = []
	len_vals = len(vals)
	for i in range(len_vals):
		val = vals[i]
		if i-nearest_neighbors >= 0 and i+nearest_neighbors < len_vals:
			print vals[i-nearest_neighbors:i+nearest_neighbors+1]
			val = np.median(vals[i-nearest_neighbors:i+nearest_neighbors+1])
		new_vals.append(val)
	return new_vals

# important if you want to take our XYZ data and use it to create images using MATPLOTLIB
def xyz_to_image(X,Y,Z):
	xset = sorted(set(X))
	yset = sorted(set(Y))
	imagex = np.ones((len(yset)+1,len(xset)+1))
	imagey = np.ones((len(yset)+1,len(xset)+1))
	imagez = np.ones((len(yset)  ,len(xset)  ))
	#
	xstep = float(xset[0])
	if len(xset) > 1:
		xstep = float(xset[1])-float(xset[0])
	ystep = float(yset[0])
	if len(yset) > 1:
		ystep = float(yset[1])-float(yset[0])
	#
	for x,y,z in zip(X,Y,Z):
		xi = xset.index(x)
		yi = yset.index(y)
		imagex[yi+1][xi+1] = x - xstep/2
		imagey[yi+1][xi+1] = y
		imagez[yi][xi] = z
	#
	for x in xset:
		imagex[0][xset.index(x)+1] = x #+ xstep
	for y in yset:
		imagey[yset.index(y)+1][0] = y
	#
	#imagex[0][0] = xset[-1] + xstep
	#imagey[0][0] = yset[-1] 
	#
	return imagex,imagey,imagez

# SOMETHING THAT PICKS DECENT VALUES FOR TICKMARKS WITHIN A GRAPH'S AXIS GIVEN THE RANGE IN DATA
# A modified version of the alpha version of the Talbot, Lin, Hanrahan tick mark generator for matplotlib.
# Described in "An Extension of Wilkinson's Algorithm for Positioning Tick Labels on Axes"
# by Justin Talbot, Sharon Lin, and Pat Hanrahan, InfoVis 2010. http://vis.stanford.edu/files/2010-TickLabels-InfoVis.pdf
# Implementation modified by Ranjan Mannige from that provided by Justin Talbot
def Extended(vmin, vmax, density_val = 1, steps = None):
	if steps is None:
		steps = [1, 5, 2, 2.5, 4, 3]
	
	def coverage(dmin, dmax, lmin, lmax):
		range = dmax-dmin
		return 1 - 0.5 * (math.pow(dmax-lmax, 2)+math.pow(dmin-lmin, 2)) / math.pow(0.1 * range, 2)
	
	def coverage_max(dmin, dmax, span):
		range = dmax-dmin
		if span > range:
			half = (span-range)/2.0
			return 1 - math.pow(half, 2) / math.pow(0.1*range, 2)
		else:
			return 1
	
	def density(k, m, dmin, dmax, lmin, lmax):
		r = (k-1.0) / (lmax-lmin)
		rt = (m-1.0) / (max(lmax, dmax) - min(lmin, dmin))
		return 2 - max( r/rt, rt/r )
	
	def density_max(k, m):
		if k >= m:
			return 2 - (k-1.0)/(m-1.0)
		else:
			return 1
	
	def simplicity(q, Q, j, lmin, lmax, lstep):
		eps = 1e-10
		n = len(Q)
		i = Q.index(q)+1
		v = 1 if ((lmin % lstep < eps or (lstep - lmin % lstep) < eps) and lmin <= 0 and lmax >= 0) else 0
		return (n-i)/(n-1.0) + v - j
	
	def simplicity_max(q, Q, j):
		n = len(Q)
		i = Q.index(q)+1
		v = 1
		return (n-i)/(n-1.0) + v - j
	
	def legibility(lmin, lmax, lstep):
		return 1
	
	def legibility_max(lmin, lmax, lstep):
		return 1
	
	def extended(dmin, dmax, m, Q=[1,5,2,2.5,4,3], only_inside=False, w=[0.25,0.2,0.5,0.05]):
		n = len(Q)
		best_score = -2.0

		j = 1.0
		while j < float('infinity'):
			for q in Q:
				sm = simplicity_max(q, Q, j)

				if w[0] * sm + w[1] + w[2] + w[3] < best_score:
					j = float('infinity')
					break

				k = 2.0
				while k < float('infinity'):
					dm = density_max(k, m)

					if w[0] * sm + w[1] + w[2] * dm + w[3] < best_score:
						break

					delta = (dmax-dmin)/(k+1.0)/j/q
					z = math.ceil(math.log(delta, 10))
		
					while z < float('infinity'):
						step = j*q*math.pow(10,z)
						cm = coverage_max(dmin, dmax, step*(k-1.0))

						if w[0] * sm + w[1] * cm + w[2] * dm + w[3] < best_score:
							break

						min_start = math.floor(dmax/step)*j - (k-1.0)*j
						max_start = math.ceil(dmin/step)*j

						if min_start > max_start:
							z = z+1
							break

						for start in range(int(min_start), int(max_start)+1):
							lmin = start * (step/j)
							lmax = lmin + step*(k-1.0)
							lstep = step

							s = simplicity(q, Q, j, lmin, lmax, lstep)
							c = coverage(dmin, dmax, lmin, lmax)
							d = density(k, m, dmin, dmax, lmin, lmax)
							l = legibility(lmin, lmax, lstep)

							score = w[0] * s + w[1] * c + w[2] * d + w[3] * l

							if score > best_score and (not only_inside or (lmin >= dmin and lmax <= dmax)):
								best_score = score
								best = (lmin, lmax, lstep, q, k)
						z = z+1
					k = k+1
			j = j+1
		return best

	#vmin, vmax = axis.get_view_interval()
	size = 5 #_figure.get_size_inches()[_which]
	# density * size gives target number of intervals,
	# density * size + 1 gives target number of tick marks,
	# the density function converts this back to a density in data units (not inches)
	# should probably make this cleaner.
	best = extended(vmin, vmax, density_val * size + 1.0, only_inside=True, w=[0.25, 0.2, 0.5, 0.05])
	locs = np.arange(best[4]) * best[2] + best[0]
	return locs

# MAKES PS CODE THAT DRAWS A POLYGON STARTING FORM (X[0],Y[0]) AND PROCEEDING TO (X[-1],Y[-1]).
def ps_draw_shape(X,Y,linewidth=0,linecolor=(0,0,0),fill=1,fillcolor=(0.6,0.6,0.6)):
	ps_text = ""
	
	# Converting single values (presumably grayscale) to RGB values 
	if not ( isinstance(linecolor, list) or isinstance(linecolor, tuple) ):
		linecolor = (linecolor,linecolor,linecolor)
	if not ( isinstance(fillcolor, list) or isinstance(fillcolor, tuple) ):
		fillcolor = (fillcolor,fillcolor,fillcolor)
	
	# starting to draw placing the cursor at the first point
	ps_text+= "%.10f %.10f moveto\n"%(X[0],Y[0])
	
	for x,y in zip(X[1:],Y[1:]):
		ps_text+= "%.10f %.10f lineto\n"%(x,y)
	ps_text+= "closepath\n"
	#	
	if fill:
		ps_text+= "gsave\n"
		ps_text+= "%.3f %.3f %.3f setrgbcolor\n"%fillcolor
		ps_text+= "fill\n"
		ps_text+= "grestore\n"
	
	#if linewidth:
	ps_text+= str(linewidth)+" setlinewidth\n"
	ps_text+= "%.3f %.3f %.3f setrgbcolor\n"%linecolor
	ps_text+= "stroke"
	
	return ps_text+"\n"

# To see how our PS/EPS system works, set 0 to 1 here
if 0:
	pstext = copy.deepcopy(postscript_template)
	pstext = pstext.replace("SCALE",str(SCALE))
	pstext = pstext.replace("XMAX",str(400*SCALE))
	pstext = pstext.replace("YMAX",str(250*SCALE))
	usermaterial = ""
	# Drawing a red box
	usermaterial+= ps_draw_shape([50,200,200,50],[50,50,200,200],0.4,linecolor = 0.001,fillcolor=(1,0,0))
	# Drawing a blue box
	usermaterial+= ps_draw_shape([50+150,200+150,200+150,50+150],[50,50,200,200],0,linecolor = 0,fillcolor=(0,0,1))
	pstext = pstext.replace("USERMATERIAL",usermaterial)
	f = open("deme.eps","w")
	f.write(pstext)
	f.close()
	os.system("evince deme.eps")
	exit()

# DRAWS POSTSCRIPT/EPS FILES USING INFORMATION IN Xoriginal, Yoriginal, Zoriginal
def make2Dfigure(Xoriginal,Yoriginal,Zoriginal,fn=0,xlim=[],ylim=[],zlim=[],cmap=plt.get_cmap('gray_r'),xscaling=defaultXscaling,xtitle="",ytitle="",xticks=[],xlabels=False,yticks=[],ylabels=False,horizontallines=[],verticallines=[],title="",showeps=0):
	# Drawing the frame
	frame_thickness = 0.06
	
	if fn == 0 or str(fn)[-len(".eps"):] !=".eps":
		print "NO EPS FILENAME GIVEN. EXITING."
		exit()
	
	if len(horizontallines):
		newadditionlines = {}
		if type(horizontallines) is list:
			for i in horizontallines:
				newadditionlines[i] = ""
			horizontallines = copy.deepcopy(newadditionlines)
	
	if len(verticallines):
		newadditionlines = {}
		if type(verticallines) is list:
			for i in verticallines:
				newadditionlines[i] = ""
			verticallines = copy.deepcopy(newadditionlines)
	
	# Making sure that we do not modify the original XYZ
	X = copy.deepcopy(Xoriginal)
	Y = copy.deepcopy(Yoriginal)
	Z = copy.deepcopy(Zoriginal)
	
	if 1: # PREPROCESSING
		# FIRST WE FIND AND THEN DELETE ALL SETS OF VALUES THAT FALL 
		# OUTSIDE OF xlim or ylim
		deleteindices = []
		if len(xlim):
			for i in range(len(X)):
				if X[i] < xlim[0] or X[i] > xlim[1]:
					deleteindices.append(i)
		if len(ylim):
			for i in range(len(Y)):
				if Y[i] < ylim[0] or Y[i] > ylim[1]:
					deleteindices.append(i)
		deleteindices = sorted(list(set(deleteindices)))
		deleteindices.reverse()
		if len(deleteindices):
			for i in deleteindices:
				del X[i]
				del Y[i]
				del Z[i]
		
		# IF Z VALUES ARE LARGER OR SMALLER THAN THE zmin ALLOWS, THEN WE JUST 
		# RESET THOSE VALUES TO THE BORDERS
		if len(zlim):
			zmin = sorted(zlim)[0]
			zmax = sorted(zlim)[-1]
			for i in range(len(Z)):
				if Z[i] < zmin:
					Z[i] = zmin
				elif Z[i] > zmax:
					Z[i] = zmax
				else:
					Z[i] = float(Z[i]-zmin)/float(zmax-zmin)
	
	# Getting the sorted and non-redundant values for X
	Xsorted = sorted(set(X))
	Ysorted = sorted(set(Y))
	
	# Working on setting up the dimensions of the drawing's frame 
	pageMinY = 2.5  # Arbitrary values
	pageMaxY = 14.5 #16.0 # Arbitrary values
	pageMinX = 3.0
	pageMaxX = pageMinX + float(pageMaxY-pageMinY)*xscaling
	
	xmin = Xsorted[0]
	xmax = Xsorted[-1]
	ymin = Ysorted[0]
	ymax = Ysorted[-1]
	
	if len(xlim):
		xmin=xlim[0]
		xmax=xlim[-1]
	if len(ylim):
		ymin=ylim[0]
		ymax=ylim[-1]
	
	title_alignment_type = "cbshow"
	if len(Xsorted) == 1:
		if xscaling == defaultXscaling: # then no new value was provided by the user
			pageMaxX = 6.0
		title_alignment_type = "rbshow"
	
	# Actual values of each position
	pageX = [] # this will be filled later
	pageY = [] # this will be filled later
	pageZ = [] # this will be filled now:
	for z in Z:
		pageZ.append(z)
	
	pageSpacingX = float(pageMaxX - pageMinX)
	if len(Xsorted) > 1:
		# Inferring the expected unit of spacing ...
		
		# Calculate single step distance in x direction
		xzero = Xsorted[0]
		xone = Xsorted[1]
		xone_minus_xzero = xone - xzero
		
		# Calculate total "distance" in x direction
		#xmin = Xsorted[0]
		#xmax = Xsorted[-1]
		xmax_minus_xmin  = xmax - xmin
		
		# Finally, calculating the spacing!
		pageSpacingX = float((xone_minus_xzero)*(pageMaxX - pageMinX))/(xone_minus_xzero+xmax_minus_xmin)
		
		# Converting real world "x" to on page "x"
		for x in X:
			pagex = pageMinX + pageSpacingX*0.5 + float(pageMaxX - pageMinX - pageSpacingX)*(x - xmin)/(xmax_minus_xmin)
			pageX.append(pagex)
	else:
		for x in X:
			pagex = pageMinX + pageSpacingX*0.5
			pageX.append(pagex)
	
	pageSpacingY = float(pageMaxY-pageMinY)
	if len(Ysorted)>1:
		# Inferring the expected unit of spacing ...
		
		# Calculate single step distance in y direction
		yzero = Ysorted[-2]
		yone = Ysorted[-1]
		
		# Calculate total "distance" in y direction
		#ymin = Ysorted[0]
		#ymax = Ysorted[-1]
		
		# Finally, calculating the y spacing!
		pageSpacingY = float((yone - yzero)*(pageMaxY - pageMinY))/(yone-yzero+ymax-ymin)
		
		# Converting real world "y" to on page "y"
		for y in Y:
			pagey = pageMinY + pageSpacingY*0.5 + float(pageMaxY - pageMinY - pageSpacingY)*(y - ymin)/(ymax-ymin)
			pageY.append(pagey)
	else:
		for y in Y:
			pagey = pageMinY + pageSpacingY*0.5
			pageY.append(pagey)	
	
	# Generating labels and ticks
	page_xticks = []
	page_yticks = []
	
	if len(set(X)) > 1:
		if xlabels == False:
			# Then calculate xlabels
			xmin_for_extended = xmin #Xsorted[0]
			xmax_for_extended = xmax #Xsorted[-1] #ax.xaxis.get_data_interval()[-1]
			xrange = xmax_for_extended-xmin_for_extended
			xmin_for_extended, xmax_for_extended = (xmin_for_extended - xrange * 0.05, xmax_for_extended + xrange * 0.05)
			xlabels = Extended(xmin_for_extended,xmax_for_extended,density_val=1.2)
			xticks  = copy.deepcopy(xlabels)			
		else:
			if len(xlabels) == len(xticks):
				# First we check if the user provided xticks with the xlabels
				pass
			elif len(xlabels) == len(Xsorted):
				# If this is the case, then we expect the two to be 
				xticks = copy.deepcopy(Xsorted)
			else:
				xticks = []
				for x in xlabels:
					if isinstance(x,float) or isinstance(x,int):
						xticks.append(x)
					else:
						print "XLABELS ARE NEITHER MATCHED BY 'XTICKS' NOR BY THE NUMBER OF X VALUES. EXITING."
						exit()
		for x in xticks:
			#print [pageMinX, pageSpacingX, pageMaxX, x, xmin, xmax]
			pagex = pageMinX + pageSpacingX*0.5 + float(pageMaxX - pageMinX - pageSpacingX)*(x - xmin)/(xmax-xmin)
			page_xticks.append(pagex)
		
			
	#
	if len(set(Y)) > 1:
		if ylabels == False:
			ymin_for_extended = ymin #Ysorted[0]
			ymax_for_extended = ymax #Ysorted[-1] #ax.yaxis.get_data_interval()[-1]
			yrange = ymax_for_extended-ymin_for_extended
			ymin_for_extended, ymax_for_extended = (ymin_for_extended - yrange * 0.05, ymax_for_extended + yrange * 0.05)		
			ylabels = Extended(ymin_for_extended,ymax_for_extended,density_val=1.2)
			yticks  = copy.deepcopy(ylabels)			
		else:
			if len(ylabels) == len(yticks):
				# First we check if the user provided xticks with the xlabels
				pass
			elif len(ylabels) == len(Ysorted):
				# If this is the case, then we expect the two to be 
				yticks = copy.deepcopy(Ysorted)
			else:
				yticks = []
				for y in ylabels:
					if isinstance(y,float) or isinstance(y,int):
						yticks.append(y)
					else:
						print "YLABELS ARE NEITHER MATCHED BY 'YTICKS' NOR BY THE NUMBER OF Y VALUES. EXITING."
						exit()
		for y in yticks:
			#print [pageMinY, pageSpacingY, pageMaxY, y, ymin, ymax]
			pagey = pageMinY + pageSpacingY*0.5 + float(pageMaxY - pageMinY - pageSpacingY)*(y - ymin)/(ymax-ymin)
			page_yticks.append(pagey)
	
	usermaterial = ""
	
	# DRAWING THE FRAME
	xmin = round(pageMinX - frame_thickness/2.0,3)
	xmax = round(pageMaxX + frame_thickness/2.0,3)
	ymin = round(pageMinY - frame_thickness/2.0,3)
	ymax = round(pageMaxY + frame_thickness/2.0,3)
	Xs = [xmin,xmax,xmax,xmin]
	Ys = [ymin,ymin,ymax,ymax]
	usermaterial += ps_draw_shape(Xs,Ys,linewidth=frame_thickness,linecolor=0,fill=1,fillcolor=1)
	
	
	zmin = round(sorted(pageZ)[0],0)
	zmax = round(sorted(pageZ)[-1],0)
	for x,y,z in zip(pageX,pageY,pageZ):
		
		round_number = 10
		
		xmin = round(x-pageSpacingX/2.0,round_number)
		xmax = round(x+pageSpacingX/2.0,round_number)
		ymin = round(y-pageSpacingY/2.0,round_number)
		ymax = round(y+pageSpacingY/2.0,round_number)
		
		Xs = [xmin,xmax,xmax,xmin]
		Ys = [ymin,ymin,ymax,ymax]
		
		#heat = cmap(z)
		zmax_minus_zmin = float(zmax-zmin)
		heat = cmap(0)
		if zmax_minus_zmin:
			heat = cmap(float(z-zmin)/zmax_minus_zmin)
		usermaterial += ps_draw_shape(Xs,Ys,linewidth=frame_thickness/10000,linecolor=heat[:3],fillcolor=heat[:3])
		#                                             |
		#                                          This extra padding is added 
		#                                          so that the some PDF interpreters
		#                                          dot not interpret the almost-zero 
		#                                          space between squares to be white.	
		
	
	for i in range(len(page_yticks)):
		if len(str(ylabels[i])):
			l = ylabels[i]
			pagey = page_yticks[i]
			y = yticks[i]
			
			if not len(ylim):
				ylim = [Ysorted[0],Ysorted[-1]]
			
			if y <= ylim[1] and y >= ylim[0]:
				Xs = [pageMinX -frame_thickness/2, pageMinX + frame_thickness*2]
				Ys = [pagey                      , pagey                       ]
				
				usermaterial += ps_draw_shape(Xs,Ys,linewidth=frame_thickness*0.8,linecolor=0,fill=0,fillcolor=1)
				
				l = str(l)
				if "." in l:
					l = l.rstrip("0").rstrip(".")
				
				textx = pageMinX - pageMinX*0.05
				texty = pagey
				usermaterial += "%.10f %.10f moveto\n"%(textx,texty)
				usermaterial += "("+str(l)+")  rcshow\n"
	
	if len(xtitle):
		# Adding the x label!
		textx = pageMinX + float(pageMaxX - pageMinX)*0.5
		texty = float(pageMinY)/2.0
		usermaterial += "%.10f %.10f moveto\n"%(textx,texty)
		usermaterial += "("+str(xtitle)+")  ctshow\n"
	#
	if len(ytitle):
		# Adding the y label!
		textx = float(pageMinX)/3.0
		texty = pageMinY + float(pageMaxY-pageMinY)/2.0
		usermaterial += "%.10f %.10f moveto\n"%(textx,texty)
		usermaterial += "("+str(ytitle)+")  ccshowrotate\n"
	#
	if len(title):
		# Adding the title!
		textx = pageMinX + float(pageMaxX-pageMinX)/2.0
		texty = float(pageMaxY)+float(pageMaxY)*0.02
		if title_alignment_type == "rbshow":
			textx = pageMaxX
		usermaterial += "%.10f %.10f moveto\n"%(textx,texty)
		usermaterial += "("+str(title)+") "+title_alignment_type+"\n"
	# 
	for i in range(len(page_xticks)):
		if len(str(xlabels[i])):
			l = xlabels[i]
			pagex = page_xticks[i]
			x = xticks[i]

			if not len(xlim):
				xlim = [Xsorted[0],Xsorted[-1]]
			
			if x <= xlim[1] and x >= xlim[0]:
				Xs = [pagex                         , pagex                       ]
				Ys = [pageMinY - frame_thickness/2.0, pageMinY + frame_thickness*2]
				usermaterial += ps_draw_shape(Xs,Ys,linewidth=frame_thickness*0.8,linecolor=0,fill=0,fillcolor=1)
				
				textx = pagex
				texty = pageMinY - pageMinY*0.1 
				
				l = str(l)
				if "." in l:
					l = l.rstrip("0").rstrip(".")
				
				usermaterial += "%.10f %.10f moveto\n"%(textx,texty)
				usermaterial += "("+str(l)+")  ctshow\n"
	
	if len(verticallines):
		for x,xlabel in verticallines.items():
			x = float(x)
			if x >= xlim[0] or x <= xlim[1]:
				pagex = pageMinX + pageSpacingX*0.5 + float(pageMaxX - pageMinX - pageSpacingX)*(x - Xsorted[0])/(Xsorted[-1]-Xsorted[0])
				
				yend = pageMaxY
				
				Xs = [pagex   , pagex]
				Ys = [pageMinY, yend]
				
				usermaterial += ps_draw_shape(Xs,Ys,linewidth=frame_thickness*0.4,linecolor=0,fill=0,fillcolor=1)
				
				#c.text(xend, pagey, r"\Large "+str(ylabel), [text.halign.boxleft,text.valign.middle])
				
				if len(xlabel):
					textx = pagex
					texty = yend + frame_thickness*6
					usermaterial += "%.10f %.10f moveto\n"%(textx,texty)
					usermaterial += "("+str(xlabel)+")  cbshow\n"

	if len(horizontallines):
		for y,ylabel in horizontallines.items():
			y = float(y)
			if y >= ylim[0] or y <= ylim[1]:
				pagey = pageMinY + pageSpacingY*0.5 + float(pageMaxY - pageMinY - pageSpacingY)*(y - Ysorted[0])/(Ysorted[-1]-Ysorted[0])
				
				xstart = pageMaxX #- frame_thickness/2
				xend   = pageMaxX + frame_thickness*3
				
				Xs = [xstart, xend ]
				Ys = [pagey , pagey]
				
				usermaterial += ps_draw_shape(Xs,Ys,linewidth=frame_thickness*0.4,linecolor=0,fill=0,fillcolor=1)
				
				if len(ylabel):
					textx = xend + frame_thickness*5
					texty = pagey
					usermaterial += "%.10f %.10f moveto\n"%(textx,texty)
					usermaterial += "("+str(ylabel)+")  lcshow\n"
	
	
	# DRAWING THE FRAME
	xmin = round(pageMinX - frame_thickness/2.0,3)
	xmax = round(pageMaxX + frame_thickness/2.0,3)
	ymin = round(pageMinY - frame_thickness/2.0,3)
	ymax = round(pageMaxY + frame_thickness/2.0,3)
	Xs = [xmin,xmax,xmax,xmin]
	Ys = [ymin,ymin,ymax,ymax]
	usermaterial += ps_draw_shape(Xs,Ys,linewidth=frame_thickness,linecolor=0,fill=0,fillcolor=1)
	
	
	
	pstext = copy.deepcopy(postscript_template)
	pstext = pstext.replace("SCALE",str(SCALE))
	pstext = pstext.replace("XMAX" ,str(SCALE*(pageMaxX+1)))
	pstext = pstext.replace("YMAX" ,str(SCALE*(pageMaxY+1)))
	pstext = pstext.replace("USERMATERIAL",usermaterial)
	
	f = open(fn,"w")
	f.write(pstext)
	f.close()
	if showeps:#show_graphs:
		os.system("evince "+fn)
		raw_input()
		
	#print "#WRITTEN TO:",fn
	
	if 0: #Then also do a matplotlib version of the graph (This is in case someone wants to use matplolib in stead)
		plt.clf() 
		imagex,imagey,imagez = xyz_to_image(X,Y,Z)
		xset = sorted(set(X))
		yset = sorted(set(Y))
		extent = [xset[0]-0.5, xset[-1], yset[0], yset[-1]]
		plt.pcolor(imagex,imagey,imagez, cmap=cmap)
		#plt.yticks(np.arange(0,100),range(0,10))
		#plt.xticks(np.arange(0.5,10.5),range(0,10))
		plt.rc("font", size=default_fontsize)
		plt.axis(extent)
		plt.axes().set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/xscaling)
		plt.show()
	
	return 1


"""
KEYS FOR DSSP:
G = 3-turn helix (3_10-helix). Min length 3 residues.
H = 4-turn helix (alpha-helix). Min length 4 residues.
I = 5-turn helix (pi-helix). Min length 5 residues.
E = extended strand in parallel and/or anti-parallel beta-sheet conformation. Min length 2 residues.
----- Less "exciting" motifs
T = hydrogen bonded turn (3, 4 or 5 turn)
B = residue in isolated beta-bridge (single pair beta-sheet hydrogen bond formation)
S = bend (the only non-hydrogen-bond based assignment).
C = coil (residues which are not in any of the above conformations).
"""
def get_resid_to_dssp(fn):
	# Create the dssp report
	dssp_output = fn+".dssp"
	os.system("./dssp-2.0.4-linux-amd64 -i "+str(fn)+" -o "+dssp_output)
	
	f = open(dssp_output)
	dssp_block = f.read()
	f.close()
	
	# populating the DSSP secondary structure records for this domain
	resno_to_sstype = {}
	dssp_block = dssp_block.split("#  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA")[-1].rstrip()
	for line in dssp_block.split("\n"):
		line = line.rstrip()
		#['  304        !              0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0'] 
		if len(line):
			residue_str = line[5:10].lstrip().rstrip()
			if len(residue_str):
				#print [line],dsspfn
				#print line[5:10].lstrip().rstrip()
				residue = int(residue_str)
				ss_type = line[16]
				resno_to_sstype[residue]=ss_type
	
	return resno_to_sstype

# PARAMTERS FOR THE RAMACHANDRAN NUMBER
# The ranges for phi and psi are [-180,180]. 
# Any other value will be garbled (so, remember to 
# convert your angles so that it fits this range.
import math
bound = 360.0 # This does not chang
sigma = 10.0 # This is sigma in the manuscript (Mannige, Kundu, Whitelam, 2016)
multiplier = int(round(sigma*(bound*(2.0**0.5)),0)) # For internal reference
multiplier_by_two = round(sigma*(bound*(2.0**0.5)),0)/2.0 # For internal reference

# To get the raw Ramachandran number from phi and psi:
def raw_ramachandran_number_collapse(phi,psi):
	phi = float(phi)
	psi = float(psi)
	a = round(sigma*(phi-psi + bound)/math.sqrt(2.0),0)
	b = round(sigma*(phi+psi + bound)/math.sqrt(2.0),0)
	return a + b*multiplier
#
# Original phi and psi values are easily obtained from the *raw* Ramachandran number
def raw_ramachandran_number_expand(z):
	z = float(z)
	phi = (math.sqrt(2.0)*np.mod(z,multiplier)/sigma     - bound + math.sqrt(2.0)*np.floor(z / multiplier)/sigma - bound )/2.0
	psi = (math.sqrt(2.0)*np.floor(z / multiplier)/sigma - bound - math.sqrt(2.0)*np.mod(z, multiplier)/sigma + bound)/2.0
	return phi,psi

# First getting the lowest and highest possible unnormalized R numbers
raw_R_min = raw_ramachandran_number_collapse(-180,-180)
raw_R_max = raw_ramachandran_number_collapse(180,180)

# However, we normally does not 
# If signed == 0, then we have the normal R number ranging from 0 to 1.
# If signed != 0, then the R number ranges from -1 to 1. Here, those R
#                 numbers associated with regions to the RIGHT of the 
#                 positive-sloped diagonal are multipled by -1.
def normalized_ramachandran_number(phi,psi, signed=0):
	# was this: 
	#raw_R  = raw_ramachandran_number_collapse(phi,psi)
	#final_r = float(raw_R - raw_R_min)/(raw_R_max-raw_R_min)
	phi = float(phi)
	psi = float(psi)
	a = round(sigma*(phi-psi + bound)/math.sqrt(2.0),0)
	b = round(sigma*(phi+psi + bound)/math.sqrt(2.0),0)
	raw_r = a + b*multiplier
	final_r = float(raw_r - raw_R_min)/float(raw_R_max - raw_R_min)
	if signed:
		if a >= multiplier_by_two:
			final_r = final_r * -1.0
	return final_r

#To see how the Ramachandran number function works, set "if 0:" to "if 1:"
if 0:
	current_phi = -55
	current_psi = -60
	# unsigned R
	R = normalized_ramachandran_number(current_phi,current_psi,signed=0)
	print "       R(phi="+str(current_phi)+",psi="+str(current_psi)+") = "+str(R)
	# signed R
	R = normalized_ramachandran_number(current_phi,current_psi,signed=1)
	print "signed R(phi="+str(current_phi)+",psi="+str(current_psi)+") = "+str(R)
	print
	current_phi = -65
	current_psi = -60
	# unsigned R
	R = normalized_ramachandran_number(current_phi,current_psi,signed=0)
	print "       R(phi="+str(current_phi)+",psi="+str(current_psi)+") = "+str(R)
	# signed R
	R = normalized_ramachandran_number(current_phi,current_psi,signed=1)
	print "signed R(phi="+str(current_phi)+",psi="+str(current_psi)+") = "+str(R)
	exit()

#
def histogram2d(X_vals,Y_vals,cmap=plt.cm.Blues,xyrange=[],title="",fn="", pairtype='rho'):
	print "fn:",fn
	norm = False
	normed = True
	
	if not len(xyrange):
		xyrange=[[min(X_vals),max(X_vals)],[min(Y_vals),max(Y_vals)]]
	# symmetrize
	if 0:
		for i in range(len(X_vals)):
			X_vals.append(Y_vals[i])
			Y_vals.append(X_vals[i])
	
	
	
	if forcedmax is not False:
		xyrange[0][1] = forcedmax
		xyrange[1][1] = forcedmax
	
	
	if forcedmin is not False:
		xyrange[0][0] = forcedmin
		xyrange[1][0] = forcedmin
	
	plt.clf() #show_graphs
	
	print "#fn:",fn
	print "#xyrange:",xyrange
	
	bins = 40.0
	lowest_fraction = 0.04
	
	h, x, y, p = plt.hist2d(X_vals, Y_vals, bins = bins, cmap=cmap,range=xyrange, normed=normed)#norm=norm)
	
	extent = [xyrange[0][0], xyrange[0][1], xyrange[1][0], xyrange[1][1]]
	if 1:
		plt.clf() #show_graphs
		if 0:
			a = plt.imshow(h, interpolation = "gaussian",cmap=cmap)
		else:
			plt.clf()
			step = float(h.max()-h.min())/bins
			
			#contour_values = list(np.arange(h.min(),h.max()+step/2,step)) # for Fig 5
			contour_values = [-0.00001]+list(np.arange(float(h.max())*lowest_fraction,h.max()+step/2,step)) # for Fig 5
			#h.max()/3.0
			#contour_values = np.array([-0.00001,h.max()/3.0,h.max()+0.0001]) # for Fig 5
			
			norm = colors.BoundaryNorm(contour_values, cmap.N)
			a = plt.contourf(h, contour_values, zdir='z', cmap=cmap, linewidth=0.0, norm=norm, extent=extent,antialiased=False)#,offset=h.min)
		
		if norm: #= 0:
			a.set_norm(norm)
		plt.axis(extent)
		
		xyrangestep = float(max(xyrange[0])-min(xyrange[0]))/2
		plt.xticks(np.arange(min(xyrange[0]), max(xyrange[0])+xyrangestep/2, xyrangestep))
		xyrangestep = float(max(xyrange[1])-min(xyrange[1]))/2
		plt.yticks(np.arange(min(xyrange[1]), max(xyrange[1])+xyrangestep/2, xyrangestep),rotation='vertical')
		
		plt.rc("font", size=default_fontsize)
		plt.title(title)
		plt.xlabel('$\\'+pairtype+'_1$', fontsize=35)
		plt.ylabel('$\\'+pairtype+'_2$', fontsize=35)
		cbar = plt.colorbar(ticks=[])
		#cbar.ax.set_xlabel("$\pi$", fontsize=22)
		if len(fn):
			plt.savefig(fn+"_raw.png", format='png',bbox_inches='tight')
			plt.savefig(fn+"_raw.pdf", format='pdf',bbox_inches='tight')
			plt.savefig(fn+"_raw.eps", format='eps',bbox_inches='tight')
		if show_graphs:
			plt.show()
	
	
	if 1:
		plt.clf() #show_graphs
		plt.rc("font", size=default_fontsize)
		plt.title(title)
		#contour_values = np.array([-0.00001,0.0001,0.1,h.max()/3.0,h.max()+0.0001]) # for Fig 5
		#contour_values = np.array([-0.00001,h.max()*lowest_fraction,h.max()+0.0001]) # for Fig 5
		contour_values = np.array([-0.00001,h.max()*0.05,h.max()*0.3,h.max()+0.0001]) # for Fig 5
		norm = colors.BoundaryNorm(contour_values, cmap.N)
		#h = scipy.ndimage.zoom(h, 2)
		
		cextent = [xyrange[0][0], xyrange[0][1], xyrange[1][0], xyrange[1][1]]
		
		cset = plt.contourf(h, contour_values, zdir='z', cmap=cmap, linewidth=0.0, norm=norm, extent=cextent)#,offset=h.min)
		cset.set_norm(norm)
		plt.contour(h, contour_values, zdir='z', colors="k", linewidths=1.2, extent=cextent)#,offset=h.min) #offset=Zmin,  #, norm=norm)
		CB = plt.colorbar(cset,spacing='proportional',ticks=[])#, shrink=0.8, extend='both')
		plt.xlabel('$\\'+pairtype+'_1$', fontsize=35)
		plt.ylabel('$\\'+pairtype+'_2$', fontsize=35)
		
		xyrangestep = float(max(xyrange[0])-min(xyrange[0]))/2
		plt.xticks(np.arange(min(xyrange[0]), max(xyrange[0])+xyrangestep/2, xyrangestep))
		xyrangestep = float(max(xyrange[1])-min(xyrange[1]))/2
		plt.yticks(np.arange(min(xyrange[1]), max(xyrange[1])+xyrangestep/2, xyrangestep),rotation='vertical')
		
		#plt.axis([xyrange[0][0], xyrange[0][1], xyrange[1][0], xyrange[1][1]])
		#cbar.ax.set_xlabel("$\pi$", fontsize=22)
		if len(fn):
			plt.savefig(fn+"_contour.png", format='png',bbox_inches='tight')
			plt.savefig(fn+"_contour.pdf", format='pdf',bbox_inches='tight')
			plt.savefig(fn+"_contour.eps", format='eps',bbox_inches='tight')
		if show_graphs:
			plt.show()
	#exit()
	#plt.savefig(output_fn_base+".eps")
	#plt.savefig(output_fn_base+".png",dpi=200)
	#plt.savefig(output_fn_base+".eps",)

# radius of gyration
def calculate_rg(ref_structure,atoms_to_be_used=[]):
	ref_model = ref_structure[0]
	#
	ref_atoms = []
	# Iterate of all chains in the model in order to find all residues
	for ref_chain in ref_model:
		# Iterate of all residues in each model in order to find proper atoms
		for ref_res in ref_chain:
			# Check if residue number ( .get_id() ) is in the list
			use_residue = 0
			if not len(atoms_to_be_used):
				use_residue = 1
			else:
				if ref_res.get_id()[1] in atoms_to_be_used:
					use_residue = 1
			if use_residue:
				ref_atoms.append(ref_res['CA'])
				#ref_atoms.append(ref_res['C'])
				#ref_atoms.append(ref_res['N'])
	
	atoms = []#numpy.zeros((len(ref_atoms),3))
	for i in ref_atoms:
		atoms.append(i.get_coord())
	atoms = numpy.array(atoms)
	CoM = sum(atoms)/len(atoms)
	q = 0
	atoms = atoms - CoM 
	return numpy.sqrt((atoms ** 2).sum() / len(atoms))
	#return math.sqrt(q / len(atoms)) #(2*(len(atoms)+1)**2)
#

# end-to-end distance
def calculate_end_to_end_distance(ref_structure,atoms_to_be_used=[]):
	ref_model = ref_structure[0]
	#
	distances = []
	ref_atoms = []
	# Iterate of all chains in the model in order to find all residues
	for ref_chain in ref_model:
		# Iterate of all residues in each model in order to find proper atoms
		for ref_res in ref_chain:
			# Check if residue number ( .get_id() ) is in the list
			use_residue = 0
			if not len(atoms_to_be_used):
				use_residue = 1
			else:
				if ref_res.get_id()[1] in atoms_to_be_used:
					use_residue = 1
			if use_residue:
				ref_atoms.append(ref_res['CA'])
				distances.append(ref_res['CA'].get_coord())
				#ref_atoms.append(ref_res['C'])
				#ref_atoms.append(ref_res['N'])
	
	distance = numpy.linalg.norm(numpy.array(distances[-1])-numpy.array(distances[0]))
	return distance
	"""
	atoms = []#numpy.zeros((len(ref_atoms),3))
	for i in ref_atoms:
		atoms.append(i.get_coord())
	atoms = numpy.array(atoms)
	CoM = sum(atoms)/len(atoms)
	q = 0
	atoms = atoms - CoM 
	return numpy.sqrt((atoms ** 2).sum() / len(atoms))
	#return math.sqrt(q / len(atoms)) #(2*(len(atoms)+1)**2)
	"""
#
def calculate_rmsd(ref_structure,sample_structure,atoms_to_be_aligned1=[],atoms_to_be_aligned2=[]):
	ref_model = ref_structure[0]
	sample_model = sample_structure[0]
	#ref_structure = pdb_parser.get_structure("reference", "1D3Z.pdb")
	#sample_structure = pdb_parser.get_structure("samle", "1UBQ.pdb")
	# Make a list of the atoms (in the structures) you wish to align.
	# In this case we use CA atoms whose index is in the specified range
	ref_atoms = []
	sample_atoms = []
	
	#if not len(atoms_to_be_aligned1):
	#	atoms_to_be_aligned1 = range(1, len(ref_model)+1)
	#print atoms_to_be_aligned1
	#exit()
	
	if not len(atoms_to_be_aligned2):
		atoms_to_be_aligned2 = atoms_to_be_aligned1
	
	# Iterate of all chains in the model in order to find all residues
	for ref_chain in ref_model:
		# Iterate of all residues in each model in order to find proper atoms
		for ref_res in ref_chain:
			# Check if residue number ( .get_id() ) is in the list
			if ref_res.get_id()[1] in atoms_to_be_aligned1:
				# Append CA atom to list
				ref_atoms.append(ref_res['CA'])
				ref_atoms.append(ref_res['N'])
				ref_atoms.append(ref_res['C'])
	
	# Do the same for the sample structure
	for sample_chain in sample_model:
		for sample_res in sample_chain:
			if sample_res.get_id()[1] in atoms_to_be_aligned2:
				sample_atoms.append(sample_res['CA'])
				sample_atoms.append(sample_res['N'])
				sample_atoms.append(sample_res['C'])	
	# Now we initiate the superimposer:
	super_imposer = Bio.PDB.Superimposer()
	super_imposer.set_atoms(ref_atoms, sample_atoms)
	super_imposer.apply(sample_model.get_atoms())
	
	return super_imposer.rms,sample_structure


#
def calculate_dihedral_angle(p):
	b = p[:-1] - p[1:]
	b[0] *= -1
	v = numpy.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
	# Normalize vectors
	v /= numpy.sqrt(numpy.einsum('...i,...i', v, v)).reshape(-1,1)
	b1 = b[1] / numpy.linalg.norm(b[1])
	x = numpy.dot(v[0], v[1])
	m = numpy.cross(v[0], b1)
	y = numpy.dot(m, v[1])
	d = numpy.degrees(numpy.arctan2( y, x ))
	return d

aa_three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

from Bio import PDB
def read_pdb(fn,signed=0):
	p=PDB.PDBParser() #(PERMISSIVE=1)
	structure=p.get_structure(fn[:-len(".pdb")], fn)
	#for model in structure:
	#	print [model.id]
	model_to_chain_to_resno_atom_to_vals = {}
	# structure (models) -> model -> chain -> residue -> atom
	for model in structure:
		model_number = model.id
		#
		if not model_number in model_to_chain_to_resno_atom_to_vals:
			model_to_chain_to_resno_atom_to_vals[model_number] = {}
		#
		for chain in model:
			segname = chain.id
			if not segname in model_to_chain_to_resno_atom_to_vals[model_number]:
				model_to_chain_to_resno_atom_to_vals[model_number][segname] = {}
			
			for residue in chain:
				resname = residue.resname
				resno   = residue.id[1]
				
				#
				i = resno
				im = i-1
				ip = i+1
				
				neighbors_found = 1
				try:
					a = structure[model_number][segname][im]["C"].coord
					b = structure[model_number][segname][i]["N"].coord
					c = structure[model_number][segname][i]["CA"].coord
					d = structure[model_number][segname][i]["C"].coord
					e = structure[model_number][segname][ip]["N"].coord
					
					if not resno in model_to_chain_to_resno_atom_to_vals[model_number][segname]:
						model_to_chain_to_resno_atom_to_vals[model_number][segname][resno] = {}
					
					model_to_chain_to_resno_atom_to_vals[model_number][segname][resno]["resname"] = resname
					singleaa = resname
					if resname in aa_three_to_one:
						singleaa = aa_three_to_one[resname]
					model_to_chain_to_resno_atom_to_vals[model_number][segname][resno]["aa"] = singleaa
					model_to_chain_to_resno_atom_to_vals[model_number][segname][i]["n"]  = b
					model_to_chain_to_resno_atom_to_vals[model_number][segname][i]["ca"] = c
					model_to_chain_to_resno_atom_to_vals[model_number][segname][i]["c"]  = d
				
				except:
					neighbors_found = 0
				
				if neighbors_found: #im in resids and ip in resids:
					#phi = calculate_dihedral_angle(numpy.array([a,b,c,d]))
					#psi = calculate_dihedral_angle(numpy.array([b,c,d,e]))
					#rho = normalized_ramachandran_number(phi,psi,signed)
					
					phi_atoms = numpy.array([a,b,c,d])
					psi_atoms = numpy.array([b,c,d,e])
					
					# CALCULATING PHI
					# identical to: phi = calculate_dihedral_angle(numpy.array([a,b,c,d]))
					p = phi_atoms
					b = p[:-1] - p[1:]
					b[0] *= -1
					v = numpy.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
					# Normalize vectors
					v /= numpy.sqrt(numpy.einsum('...i,...i', v, v)).reshape(-1,1)
					b1 = b[1] / numpy.linalg.norm(b[1])
					x = numpy.dot(v[0], v[1])
					m = numpy.cross(v[0], b1)
					y = numpy.dot(m, v[1])
					phi = numpy.degrees(numpy.arctan2( y, x ))
					
					# CALCULATING PSI
					# identical to: psi = calculate_dihedral_angle(numpy.array([b,c,d,e]))
					p = psi_atoms
					b = p[:-1] - p[1:]
					b[0] *= -1
					v = numpy.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
					# Normalize vectors
					v /= numpy.sqrt(numpy.einsum('...i,...i', v, v)).reshape(-1,1)
					b1 = b[1] / numpy.linalg.norm(b[1])
					x = numpy.dot(v[0], v[1])
					m = numpy.cross(v[0], b1)
					y = numpy.dot(m, v[1])
					psi = numpy.degrees(numpy.arctan2( y, x ))
					
					# CALCULATING RHO
					# identical (in outcome) to: rho = normalized_ramachandran_number(phi,psi,signed)
					sqrttwo        = math.sqrt(2.0)
					psi_plus_bound = psi + bound
					a = round(sigma*(phi-psi_plus_bound)/sqrttwo,0)
					b = round(sigma*(phi+psi_plus_bound)/sqrttwo,0)
					raw_r   = a + b*multiplier
					rho = float(raw_r - raw_R_min)/(raw_R_max - raw_R_min)
					if signed:
						if a >= multiplier_by_two:
							rho = rho * -1.0
					#
					model_to_chain_to_resno_atom_to_vals[model_number][segname][i]["phi"] = phi
					model_to_chain_to_resno_atom_to_vals[model_number][segname][i]["psi"] = psi
					model_to_chain_to_resno_atom_to_vals[model_number][segname][i]["r"] = rho
					#
				#
			#
		#
		if not len(model_to_chain_to_resno_atom_to_vals[model_number]):
			del model_to_chain_to_resno_atom_to_vals[model_number]
	return model_to_chain_to_resno_atom_to_vals


# OLD VERSION (IN HOUSE). IT IS FASTER THAN THE CURRENT "read_pdb", WHICH IS BIOPDB RUN, BUT IT IS NOT 
# AS WELL TESTED.
def read_pdb_old(fn,signed=0):
	"""
	ATOM     10 1H   LYS A   1       0.763   3.548  -0.564
	ATOM     11 2H   LYS A   1       1.654   2.664   0.488
	ATOM    482  N   PRO A  61      27.194  -5.761  14.684  1.00  9.09           N  
	ATOM      2  CA  BLYSX   1     -77.937 -26.325   6.934  1.00  0.00      U1    
	ATOM      3  CB  BLYSX   1     -79.612 -24.499   7.194  1.00  0.00      U1    
	ATOM      4  CE  BLYSX   1     -80.894 -24.467   8.039  1.00  0.00      U1    
	ATOM      5  NZ  BLYSX   1     -80.687 -24.160   9.434  1.00  0.00      U1    
	ATOM      2  HT1 MET U   1       0.208   0.762 -12.141  0.00  0.00      UBIQ  
	ATOM      3  HT2 MET U   1      -1.052  -0.551 -12.281  0.00  0.00      UBIQ  
	          |   |   |  |   |        |       |       |                     |
	     atomno   |   |  |   |        x       y       z                 segname
	       atom type  |  |   |                                          (CHAIN)
	            restype  |   3resno
	                 chainID
	"""
	
	f = open(fn,"r")
	pdbblock = f.read()
	f.close()
	
	getlines = re.compile("ATOM\s+(?P<atomno>\d+)\s+(?P<atomtype>\S+)\s+.(?P<resname>...)..\s+(?P<resno>\d+)\s+(?P<x>\-*\d+\.*\d*)\s+(?P<y>\-*\d+\.*\d*)\s+(?P<z>\-*\d+\.*\d*).{17}(?P<segname>.{5})",re.M)
	getlines_short = re.compile("ATOM\s+(?P<atomno>\d+)\s+(?P<atomtype>\S+)\s+(?P<resname>...).(?P<segname>.)\s+(?P<resno>\d+)\s+(?P<x>\-*\d+\.*\d*)\s+(?P<y>\-*\d+\.*\d*)\s+(?P<z>\-*\d+\.*\d*)",re.M)
	
	resnos = []
	#models = pdbblock.split("\nEND\n")
	models = re.split("\nEND|\nMODEL|\nTER",pdbblock)
	
	model_number = 0
	model_to_chain_to_resno_atom_to_vals = {}
	# structure (models) -> model -> chain -> residue -> atom
	
	#t0 = time.time()
	#print "#\treading...",
	for model_index in range(len(models)):
		model = models[model_index]
		if len(model.rstrip()) > 1:
			model_number+=1
			if not model_number in model_to_chain_to_resno_atom_to_vals:
				model_to_chain_to_resno_atom_to_vals[model_number] = {}
			
			segname_exists = 1
			currentlines = getlines.finditer(model)
			if not getlines.search(model):
				currentlines = getlines_short.finditer(model)
				segname_exists = 0
			
			for i in currentlines:
				vals = i.groupdict()
				atomtype = vals["atomtype"] #line[11:17].lstrip().rstrip()
				
				if atomtype=="CA" or atomtype =="N" or atomtype =="C":
					resno = int(vals["resno"]) #int(resno) #int(line[22:26].lstrip().rstrip())
					xyz = numpy.array([float(vals["x"]),float(vals["y"]),float(vals["z"])])
					
					segname = "A"
					if segname_exists:
						segname = vals["segname"].lstrip().rstrip()
					
					if not segname in model_to_chain_to_resno_atom_to_vals[model_number]:
						model_to_chain_to_resno_atom_to_vals[model_number][segname] = {}
					
					if not resno in model_to_chain_to_resno_atom_to_vals[model_number][segname]:
						model_to_chain_to_resno_atom_to_vals[model_number][segname][resno] = {}
					
					model_to_chain_to_resno_atom_to_vals[model_number][segname][resno][atomtype.lower()] = xyz
					model_to_chain_to_resno_atom_to_vals[model_number][segname][resno]["resname"] = vals["resname"]
			
			if not len(model_to_chain_to_resno_atom_to_vals[model_number]):
				del model_to_chain_to_resno_atom_to_vals[model_number]
				model_number-=1
	#
	for model in sorted(model_to_chain_to_resno_atom_to_vals.keys()):
		for chain in sorted(model_to_chain_to_resno_atom_to_vals[model].keys()):
			for resno in sorted(model_to_chain_to_resno_atom_to_vals[model][chain].keys()):
				triplet_found = 0
				if "ca" in model_to_chain_to_resno_atom_to_vals[model][chain][resno]:
					triplet_found+=1
				if "n" in model_to_chain_to_resno_atom_to_vals[model][chain][resno]:
					triplet_found+=1
				if "c" in model_to_chain_to_resno_atom_to_vals[model][chain][resno]:
					triplet_found+=1
				if triplet_found == 3:
					i = resno
					im = i-1
					ip = i+1
					
					neighbors_found = 0
					if im in model_to_chain_to_resno_atom_to_vals[model][chain]:
						if "c" in model_to_chain_to_resno_atom_to_vals[model][chain][im]:
							neighbors_found += 1
					if ip in model_to_chain_to_resno_atom_to_vals[model][chain]:
						if "n" in model_to_chain_to_resno_atom_to_vals[model][chain][ip]:
							neighbors_found += 1
					
					if neighbors_found == 2: #im in resids and ip in resids:
						a = model_to_chain_to_resno_atom_to_vals[model][chain][im]["c"] # resno_to_coordC[before]
						b = model_to_chain_to_resno_atom_to_vals[model][chain][i]["n"] # resno_to_coordN[current]
						c = model_to_chain_to_resno_atom_to_vals[model][chain][i]["ca"] #resno_to_coordCA[current]
						d = model_to_chain_to_resno_atom_to_vals[model][chain][i]["c"] # resno_to_coordC[current]
						e = model_to_chain_to_resno_atom_to_vals[model][chain][ip]["n"]  # resno_to_coordN[after]
						
						phi = calculate_dihedral_angle(numpy.array([a,b,c,d]))
						psi = calculate_dihedral_angle(numpy.array([b,c,d,e]))
						rho = normalized_ramachandran_number(phi,psi,signed)
						#print rho
						#if rho < 0.5:
						#	print (phi,psi)
						
						model_to_chain_to_resno_atom_to_vals[model][chain][i]["phi"] = phi
						model_to_chain_to_resno_atom_to_vals[model][chain][i]["psi"] = psi
						model_to_chain_to_resno_atom_to_vals[model][chain][i]["r"] = rho
	
	return model_to_chain_to_resno_atom_to_vals

# OLD VERSION (IN HOUSE). IT IS FASTER THAN THE CURRENT "read_pdb", WHICH IS BIOPDB RUN, BUT IT IS NOT 
# AS WELL TESTED.
def check_pdb(fn):
	"""
	ATOM     10 1H   LYS A   1       0.763   3.548  -0.564
	ATOM     11 2H   LYS A   1       1.654   2.664   0.488
	ATOM    482  N   PRO A  61      27.194  -5.761  14.684  1.00  9.09           N  
	ATOM      2  CA  BLYSX   1     -77.937 -26.325   6.934  1.00  0.00      U1    
	ATOM      3  CB  BLYSX   1     -79.612 -24.499   7.194  1.00  0.00      U1    
	ATOM      4  CE  BLYSX   1     -80.894 -24.467   8.039  1.00  0.00      U1    
	ATOM      5  NZ  BLYSX   1     -80.687 -24.160   9.434  1.00  0.00      U1    
	ATOM      2  HT1 MET U   1       0.208   0.762 -12.141  0.00  0.00      UBIQ  
	ATOM      3  HT2 MET U   1      -1.052  -0.551 -12.281  0.00  0.00      UBIQ  
	          |   |   |  |   |        |       |       |                     |
	     atomno   |   |  |   |        x       y       z                 segname
	       atom type  |  |   |                                          (CHAIN)
	            restype  |   resno
	                 chainID
	"""
	
	chainIDindex = 21
	chainIDindexMinusOne = chainIDindex-1
	lenATOM = len("ATOM ")
	
	chainIDpossibilities = ""
	chainIDpossibilities+=string.uppercase # 'A' through 'Z'.
	for i in range(10):
		chainIDpossibilities+=str(i)
	chainIDpossibilities+=string.lowercase # 'a' through 'z'.
	lenchainIDpossibilities = len(chainIDpossibilities)
	largestchainIDindex = 0
	
	made_changes = 0
	f = open(fn,"r")
	lines = f.readlines()
	f.close()
	
	segname_to_chainID = {}
	for i in range(len(lines)):
		if len(lines[i]) > 67:
			if lines[i][:lenATOM] == "ATOM ":
				chainID = lines[i][chainIDindex].rstrip()
				chainIDspacebefore = lines[i][chainIDindexMinusOne].rstrip()
				if len(chainIDspacebefore): # This is because some CHARMM sidechains have four letters, and that trips biopython
					pdb_is_possibly_problematic = 1
				
				if len(chainID)==0 or chainID=="X": # CHARMM SOMETIMES SAVES THE CHAINID AS 'X' IRRESPECTIVE OF SEGNAME
					pdb_is_possibly_problematic = 1
				#
			#
		#
	if pdb_is_possibly_problematic:
		return 0
	else:
		return 1

# Draw the mapping between phi,psi and Ramachandran number
def show_ramachandran_mapping(cmap=plt.get_cmap("Blues"),fn = "deme", stepsize=5,signed=0):
	PHI = []
	PSI = []
	RHO = []
	for phi in range(-180,181,stepsize):
		for psi in range(-180,181,stepsize):
			PHI.append(phi)
			PSI.append(psi)
			normalizedr = normalized_ramachandran_number(phi,psi,signed=signed)
			RHO.append(normalizedr)
	X=np.array(PHI)
	Y=np.array(PSI)
	Z=np.array(RHO)
	
	#color_bar_range = np.arange(0,1.01,0.01)
	color_bar_range = np.arange(rrange[0],rrange[-1]+0.005,0.01)
	
	if fn[-len(".eps"):] == ".eps":
		fn = fn[:-len(".eps")]
	
	cbfn = fn+"_colorbar.eps"
	fn = fn+".eps"
	
	print "#WRITING TO:",cbfn
	make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=cmap, fn=cbfn, xscaling=colorbarXscaling,xtitle="Key",showeps=0)
	#
	print "#WRITING TO:",fn
	make2Dfigure(X,Y,Z,fn,cmap=cmap,xscaling=1.0,xtitle="phi",ytitle="psi",xlabels=range(-180,181,90),ylabels=range(-180,181,90),showeps=0)
	#
	return 1

if __name__ == "__main__":
	if not "-pdb" in sys.argv:
		if "-h" in sys.argv or "-help" in sys.argv or "--help" in sys.argv:
			pass
		else:
			print "Must provide '-pdb' parameter. Exiting."
			exit(0)
	for i in range(len(sys.argv)):
		if sys.argv[i] == "-rmsd":
			showrmsd = 1
		if sys.argv[i] == "-signed":
			print "Using the R number with range [-1,1]"
			signed = 1
			rrange = [-1,1]
		if sys.argv[i] == "-ss":
			colortype = "SecondaryStructure" # default: chirality
		if sys.argv[i] == "-h" or sys.argv[i] == "-help" or sys.argv[i] == "--help":
			print helpme
			exit(1)
		if sys.argv[i] == "-pdb":
			if len(sys.argv) <= i+1:
				print helpme
				print "MUST PROVIDE PDB NAME."
				exit(0)
			else:
				pdbfn = str(sys.argv[i+1])
				print "# pdbfn set to:",pdbfn
		elif sys.argv[i] == "-bins":
			if len(sys.argv) <= i+1:
				helpme
				print "When using '-bins', you must provide bin number. Exiting."
				exit(0)
			else:
				if not sys.argv[i+1].isdigit():
					print helpme
					print "The -bin parameter must be a positive integer (provided: "+str(sys.argv[i+1])+") Exiting."
					exit(0)
				else:
					bins = int(sys.argv[i+1])
					print "# bins set to:",bins
					if bins == 0:
						print helpme
						print "Must have greater than 0 bins. Exiting."
						exit(0)
	
	colormap_name = colortype+'TwoColor'
	if signed:
		colormap_name = colortype+'FourColor'
	print "Using color map name:",colormap_name
	rcode_cmap = plt.get_cmap(colormap_name)
	#rcode_cmap = plt.get_cmap("deletemeSigned")
	
	pdbfn = os.path.abspath(pdbfn)
	pdbdir = os.path.dirname(pdbfn)
	pdbfilenames = []
	
	if os.path.isfile(pdbfn):
		# then this pathname leads to a FILE
		# ... so keep as is
		pdbfilenames = [pdbfn]
	elif os.path.isdir(pdbfn):
		pdbdir = pdbfn
		pdbfilenames = sorted(glob.glob(pdbdir+"/*.pdb"))
	else:
		print helpme
		exit("Either filename or directory expected. Exiting.")
	
	target_dir = pdbdir+"/reports/"
	#if len(pdbfilenames)>1:
	#	target_dir = pdbdir+"/report/"
	if not os.path.isdir(target_dir):
		os.makedirs(target_dir)
	NAME = os.path.basename(pdbfilenames[0])[:-len(".pdb")]
	target_base = target_dir.rstrip("/")+"/"+NAME
	
	#
	# To see how the (phi,psi) to Ramachandran number mapping occurs, set "if 0:" to "if 1:" below.
	if 0:
		phi_psi_to_ramachandran_number_filename = "phi_psi_to_r"
		if signed:
			phi_psi_to_ramachandran_number_filename+="_signed"
		if colortype == "SecondaryStructure":
			phi_psi_to_ramachandran_number_filename+="_ss"
		show_ramachandran_mapping(signed=signed, cmap = rcode_cmap, fn=phi_psi_to_ramachandran_number_filename,stepsize=5)
		exit()
	
	# JUST "CLEVERLY" ARRANGING THE FILENAMES, IF WE HAVE A SET OF FILENAMES RATHER THAN ONE
	# (e.g., pdbfilenames = [something2part1,something1part2,something1part1,something10part1]
	# pdbfilenames.sort() this list to: [something1part1,something1part2,something2part1,something10part1]
	REXP = re.compile( r'\d+' )
	def key_function( s ): return map(int, re.findall(REXP, s ))
	pdbfilenames.sort( key=key_function)

	print "# Parsing the PDB (structure) data"
	# structure -> model -> chain -> residue -> atom -> 'x','y','z','phi','psi','r'
	structure = {}
	for pdbfn in pdbfilenames:
		# Check if the PDB has no subunit IDs, and then check if segnames exist (right most column)
		# and renaming the subunit IDs alphabetically and then numerically
		
		# READ PDB
		if check_pdb(pdbfn):
			latest_structure = read_pdb(pdbfn,signed)
		else:
			latest_structure = read_pdb_old(pdbfn,signed)
		
		# AGAIN: structure -> model -> chain -> residue -> atom -> 'x','y','z','phi','psi','r'
		if len(structure.keys()):
			largest_model_number = max(structure.keys())
			# Then we already opened a previous file 
			for model_number in latest_structure.keys():
				# added the "1" since model numbers normally start from 0 (at least for BIOPDB):
				new_model_number = largest_model_number+model_number+1 
				# Finally, we rename the model_number key
				latest_structure[new_model_number] = latest_structure.pop(model_number)
		
		# Adding the current models (renamed if they are the same as a previous files') 
		# to the main structure dictionary
		structure.update(latest_structure)
	print "\t...done"
	
	# FIRST, GETTING THE DIFFERENT CHAIN IDs
	chains = []
	for model in sorted(structure.keys()):
		for chain in sorted(structure[model].keys()):
			chains.append(chain)
	chains = list(sorted(set(chains)))
	
	# FINALLY, WE WILL GO THROUGH EACH CHAIN AND PRODUCE GRAPHS
	batchedfilenames = {}
	for chain in chains:
		filenames = []	
		print "CHAIN:",chain
		X=[]
		Y=[]
		Z=[]
		for model in sorted(structure.keys()):
			if chain in structure[model]:
				for resno in structure[model][chain].keys():
					if "r" in structure[model][chain][resno]:
						X.append(model)
						Y.append(resno)
						Z.append(structure[model][chain][resno]["r"])
		
		for i in range(dofilter):
			Z = median_filter(Z)
		
		
		if showrcode:
			sortedZ = sorted(Z)
			
			fn = target_base+"_chain_"+chain+".rcode.eps"
			cbfn = fn+".colorbar.eps"
			filenames.append(fn); filenames.append(cbfn)
			color_bar_range = np.arange(rrange[0],rrange[-1]+0.005,0.01)
			cmap = rcode_cmap
			# colorbar
			make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=cmap, fn=cbfn, title="R",showeps=0,xscaling=colorbarXscaling)
			# data
			make2Dfigure(X,Y,Z,fn, cmap=cmap, xtitle="Model #",ytitle="Residue #",title="R (ch:`%s')"%chain,showeps=showeps)
		
		# ----------------------------------
		xy_to_z = {}
		for x,y,z in zip(X,Y,Z):
			xy_to_z[(x,y)] = z
		
		# ----------------------------------
		# RMSD from reference time
		if showrmsd:
			reference_index = 0
			reference_time = sorted(set(X))[0]
			xy_to_z_rmsd = {}
			for (x,y) in xy_to_z.keys():
				xy_to_z_rmsd[(x,y)] = abs(xy_to_z[(x,y)]-xy_to_z[(reference_time,y)])
			
			X = []
			Y = []
			Z = []
			for (x,y) in xy_to_z_rmsd.keys():
				z = xy_to_z_rmsd[(x,y)]
				X.append(x)
				Y.append(y)
				Z.append(z)
			
			if 0:#normalize
				if max(Z):
					for i in range(len(Z)):
						Z[i] = Z[i]/max(Z)
			
			fn = target_base+"_chain_"+chain+".rcode.rmsd"+str(reference_index)+".eps"
			cbfn = fn+".colorbar.eps"
			filenames.append(fn); filenames.append(cbfn)
			color_bar_range = np.arange(0,1.01,0.01)
			make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=plt.get_cmap("Reds"), ytitle="", fn=cbfn, showeps=0,xscaling=colorbarXscaling, title="RMSD")
			make2Dfigure(X,Y,Z,fn, cmap=plt.get_cmap("Reds"), xtitle="Model #",ytitle="Residue #",showeps=showeps,title="RMSD (ch:`%s')"%chain)
		
		# ----------------------------------
		# Fluctuations from previous time
		if 0:#showrmsf:
			x_sorted = sorted(set(X))
			xy_to_z_fluc = {}
			for (x,y) in xy_to_z.keys():
				current_x_index = x_sorted.index(x)
				previous_x_index = current_x_index-1
				if previous_x_index < 0:
					previous_x_index = 0
				previous_x = x_sorted[previous_x_index]
				xy_to_z_fluc[(x,y)] = abs(xy_to_z[(x,y)]-xy_to_z[(previous_x,y)])
			
			X = []
			Y = []
			Z = []
			for (x,y) in xy_to_z_fluc.keys():
				z = xy_to_z_fluc[(x,y)]
				X.append(x)
				Y.append(y)
				Z.append(z)
			if 1:#normalize
				if max(Z):
					for i in range(len(Z)):
						Z[i] = Z[i]/max(Z)
			
			fn = target_base+"_chain_"+chain+".rcode.rmsf.eps"
			cbfn = fn+".colorbar.eps"
			filenames.append(fn); filenames.append(cbfn)
			color_bar_range = np.arange(0,1.01,0.01)
			make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=plt.get_cmap("Blues"), title="RMSF", fn=cbfn, showeps=0,xscaling=colorbarXscaling)
			make2Dfigure(X,Y,Z,fn, cmap=plt.get_cmap("Blues"), xtitle="Model #",ytitle="Residue #",showeps=showeps, title="RMSF (ch:`%s')"%chain)
		# ----------------------------------
		if showhis:
			X=[]
			Y=[]
			Z=[]
			for model in sorted(structure.keys()):
				Rs = []
				if chain in structure[model]:
					for resno in structure[model][chain].keys():
						if "r" in structure[model][chain][resno]:
							Rs.append(structure[model][chain][resno]["r"])
				if len(Rs):
					step = 0.005
					bins = np.arange(rrange[0]-step/2,rrange[-1]+step,step)
					if rrange[0] < 0:
						step = 0.005*2
						bins = np.arange(rrange[0]+step/2,rrange[-1]+step,step)
					a,b = np.histogram(Rs,bins=bins)
					for i in range(len(a)):
						X.append(model)
						Y.append(float(b[i]+b[i+1])/2.0)
						Z.append(float(a[i])/np.max(a))
			#
			fn = target_base+"_chain_"+chain+".rcode.his.eps"
			cbfn = fn+".colorbar.eps"
			filenames.append(fn); filenames.append(cbfn)
			color_bar_range = np.arange(0,1.01,0.01)
			make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=plt.get_cmap("gray_r"), title="P(R)", fn=cbfn, showeps=0,xscaling=colorbarXscaling)
			
			for i in range(dofilter):
				Z = median_filter(Z)
			
			ylim = [0.29,0.71]
			if signed:
				ylim = [-0.71,0.71]
			
			make2Dfigure(X,Y,Z,fn, ylim=ylim,cmap=plt.get_cmap("gray_r"), xtitle="Model #",ytitle="R",showeps=showeps, title="P(R) (ch:`%s')"%chain)
		batchedfilenames[chain] = copy.deepcopy(filenames)
	
	if do_vmd_etc:
		
		f = open("templates/template.tex","r")
		latexblock = f.read()
		f.close()
		
		USERMATERIAL  = r"\newcommand{\w}{0.15}"+"\n"
		USERMATERIAL += r"\newcommand{\h}{0.4}"+"\n"
		USERMATERIAL += r"\newcommand{\hc}{0.4}"+"\n"
		USERMATERIAL += r"\newcommand{\centerit}[1]{\noindent\textcolor{white}{.}\hfil#1\hfil\newline}"+"\n"
		
		
		for chain in sorted(batchedfilenames.keys()):
			#USERMATERIAL += r"\centerit{{\huge The ${\mathcal{R}}$eport for \texttt{\url{%s}.pdb}}}" % NAME +"\n"
			USERMATERIAL += r"{\large The ${\mathcal{R}}$eport for \texttt{\url{%s}.pdb, chain (ch):`%s'}" %(NAME,chain) +"\n"
			USERMATERIAL += r"\vfill"+"\n"
			filenames = batchedfilenames[chain]
			for i in range(len(filenames)/2):
				fn1 = filenames[i*2]
				fn2 = filenames[i*2+1]
				USERMATERIAL += r"\noindent"
				USERMATERIAL += r"\includegraphics[height=\h\textwidth]{%s}" % fn1[:-len('.eps')]+"\n"
				USERMATERIAL += r"\includegraphics[height=\hc\textwidth]{%s}" % fn2[:-len('.eps')]+"\n"
				USERMATERIAL += r"\newline"+"\n"
			USERMATERIAL += r"\vfill"+"\n\n~\n"
			if len(batchedfilenames.keys()) > 1:
				USERMATERIAL += r"\clearpage"+"\n"
		
		outputfilename = target_base+".tex"
		print "Writing to:",outputfilename
		latexblock = latexblock.replace("USERMATERIAL",USERMATERIAL)
		f = open(outputfilename,"w")
		f.write(latexblock)
		f.close()
		
		os.system("pdflatex -shell-escape -output-directory="+os.path.dirname(outputfilename)+" "+outputfilename)
		print "OPENING:",outputfilename[:-len(".tex")]+".pdf"
		os.system("evince "+outputfilename[:-len(".tex")]+".pdf")
		os.system("rm "+outputfilename[:-len(".tex")]+".aux")
		os.system("rm "+outputfilename[:-len(".tex")]+".log")
		os.system("rm "+outputfilename[:-len(".tex")]+".out")
		os.system("rm "+outputfilename[:-len(".tex")]+"*converted-to.pdf")
	
