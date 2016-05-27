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

Usage:
python plotmap.py -pdb ProteinDatabankStructureFilename.pdb
python plotmap.py -pdb /directory/containing/pdbs/

Output (the x-axis always represents the models/structures listed in the PDB):
filename.rcode.his.eps  (y-axis: Ramachandran number (R); color: frequency of R in model)
filename.rcode.ss.eps   (y-axis: residue #; color: by secondary structure HELIX: red, SHEET: blue, PPII: cyan)
filename.rcode.raw.eps  (y-axis: residue #; color: by chirality L: Blue, D: Red: Extended: White)
filename.rcode.rmsd.eps (y-axis: residue #; color: RMSD in R from first model)
filename.rcode.rmsf.eps (y-axis: residue #; color: RMSF in R from the previous model)

Additionally, each graph is accompanied by "_colorbar.eps", which are keys.

The Ramachandran number concept is based on the manuscript:
Mannige, Kundu, Whitelam (2016) "The Ramachandran number: an order parameter for protein geometry" 
Preprint at: http://arxiv.org/abs/1511.03011

============================================
"""

#print helpme
#exit()
import sys

# Defining colors
#SCHEME 1
helixR      = (1,0,0) 
sheet       = (0,0,1)
polyproline = (0,1,1)

#helixR      = (160.0/256.0,32.0/256.0,240.0/256.0) # purple
#sheet       = (1,1,0)
#polyproline = (0,1,0)

showeps = 1
dofilter = 0
do_vmd_etc = 0

showrmsd  = 0
showrmsf  = 1
showrcode = 1
showhis   = 1

bins = 100
pdbfn = ""


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


if 0:
	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib.colors import LinearSegmentedColormap
	"""
	
	Example: suppose you want red to increase from 0 to 1 over the bottom
	half, green to do the same over the middle half, and blue over the top
	half.  Then you would use:
	
	cdict = {'red':   ((0.0,  0.0, 0.0),
			(0.5,  1.0, 1.0),
			(1.0,  1.0, 1.0)),
	
		'green': ((0.0,  0.0, 0.0),
			(0.25, 0.0, 0.0),
			(0.75, 1.0, 1.0),
			(1.0,  1.0, 1.0)),

		'blue':  ((0.0,  0.0, 0.0),
			(0.5,  0.0, 0.0),
			(1.0,  1.0, 1.0))}
	
	If, as in this example, there are no discontinuities in the r, g, and b
	components, then it is quite simple: the second and third element of
	each tuple, above, is the same--call it "y".  The first element ("x")
	defines interpolation intervals over the full range of 0 to 1, and it
	must span that whole range.  In other words, the values of x divide the
	0-to-1 range into a set of segments, and y gives the end-point color
	values for each segment.
	
	Now consider the green. cdict['green'] is saying that for
	0 <= x <= 0.25, y is zero; no green.
	0.25 < x <= 0.75, y varies linearly from 0 to 1.
	x > 0.75, y remains at 1, full green.

	If there are discontinuities, then it is a little more complicated.
	Label the 3 elements in each row in the cdict entry for a given color as
	(x, y0, y1).  Then for values of x between x[i] and x[i+1] the color
	value is interpolated between y1[i] and y0[i+1].

	Going back to the cookbook example, look at cdict['red']; because y0 !=
	y1, it is saying that for x from 0 to 0.5, red increases from 0 to 1,
	but then it jumps down, so that for x from 0.5 to 1, red increases from
	0.7 to 1.  Green ramps from 0 to 1 as x goes from 0 to 0.5, then jumps
	back to 0, and ramps back to 1 as x goes from 0.5 to 1.

	row i:   x  y0  y1
			/
		/
	row i+1: x  y0  y1

	Above is an attempt to show that for x in the range x[i] to x[i+1], the
	interpolation is between y1[i] and y0[i+1].  So, y0[0] and y1[-1] are
	never used.

	"""


	cdict1 = {'red':   ((0.0, 0.0, 0.0),
			(0.5, 0.0, 0.1),
			(1.0, 1.0, 1.0)),

		'green': ((0.0, 0.0, 0.0),
			(1.0, 0.0, 0.0)),

		'blue':  ((0.0, 0.0, 1.0),
			(0.5, 0.1, 0.0),
			(1.0, 0.0, 0.0))
		}

	cdict2 = {'red':   ((0.0, 0.0, 0.0),
			(0.5, 0.0, 1.0),
			(1.0, 0.1, 1.0)),

		'green': ((0.0, 0.0, 0.0),
			(1.0, 0.0, 0.0)),

		'blue':  ((0.0, 0.0, 0.1),
			(0.5, 1.0, 0.0),
			(1.0, 0.0, 0.0))
		}

	cdict3 = {'red': ((0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 0.8, 1.0), (0.75, 1.0, 1.0), (1.0, 0.4, 1.0)),
		'green': ((0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 0.9, 0.9), (0.75, 0.0, 0.0), (1.0, 0.0, 0.0)),
		'blue':  ((0.0, 0.0, 0.4), (0.25, 1.0, 1.0), (0.5, 1.0, 0.8), (0.75, 0.0, 0.0), (1.0, 0.0, 0.0))
		}
	helix_start = 0.31
	helix_end   = 0.39
	sheet_start = 0.45
	sheet_end   = 0.62
	polyproline_end = 0.66
	
	# c stands for color, bc stands for background color
	c1 = [0,0,0] # black
	c2 = [1,1,0] # yellow 
	c3 = [1,0,0] # red
	c4 = [0,0,1] # blue
	bc = [1,1,1] # white
	
	cdict3 = {
	#                         black  black           white  white         yellow  red             white  white         blue   blue
		'red':   ((0.00,  c1[0], c1[0]), (0.25,  bc[0], bc[0]), (0.50, c2[0], c3[0]), (0.75,  bc[0], bc[0]), (1.0, c4[0], c4[0])), 
		'green': ((0.00,  c1[1], c1[1]), (0.25,  bc[1], bc[1]), (0.50, c2[1], c3[1]), (0.75,  bc[1], bc[1]), (1.0, c4[1], c4[1])),
		'blue':  ((0.00,  c1[2], c1[2]), (0.25,  bc[2], bc[2]), (0.50, c2[2], c3[2]), (0.75,  bc[2], bc[2]), (1.0, c4[2], c4[2])) 
	}
	cmap = LinearSegmentedColormap('BlueRed2', cdict2)
	plt.register_cmap(cmap=blue_red2)
	cmap = plt.get_cmap('BlueRed2')
	
	# Make a modified version of cdict3 with some transparency
	# in the middle of the range.
	cdict4 = cdict3.copy()
	cdict4['alpha'] = ((0.0, 1.0, 1.0),
			#   (0.25,1.0, 1.0),
			(0.5, 0.3, 0.3),
			#   (0.75,1.0, 1.0),
			(1.0, 1.0, 1.0))


	# Now we will use this example to illustrate 3 ways of
	# handling custom colormaps.
	# First, the most direct and explicit:

	blue_red1 = LinearSegmentedColormap('BlueRed1', cdict1)

	# Second, create the map explicitly and register it.
	# Like the first method, this method works with any kind
	# of Colormap, not just
	# a LinearSegmentedColormap:

	blue_red2 = LinearSegmentedColormap('BlueRed2', cdict2)
	plt.register_cmap(cmap=blue_red2)

	# Third, for LinearSegmentedColormap only,
	# leave everything to register_cmap:

	plt.register_cmap(name='BlueRed3', data=cdict3)  # optional lut kwarg
	plt.register_cmap(name='BlueRedAlpha', data=cdict4)

	# Make some illustrative fake data:

	x = np.arange(0, np.pi, 0.1)
	y = np.arange(0, 2*np.pi, 0.1)
	X, Y = np.meshgrid(x, y)
	Z = np.cos(X) * np.sin(Y) * 10

	# Make the figure:

	plt.figure(figsize=(6, 9))
	plt.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)

	# Make 4 subplots:

	plt.subplot(2, 2, 1)
	plt.imshow(Z, interpolation='nearest', cmap=blue_red1)
	plt.colorbar()

	plt.subplot(2, 2, 2)
	cmap = plt.get_cmap('BlueRed2')
	plt.imshow(Z, interpolation='nearest', cmap=cmap)
	plt.colorbar()

	# Now we will set the third cmap as the default.  One would
	# not normally do this in the middle of a script like this;
	# it is done here just to illustrate the method.

	plt.rcParams['image.cmap'] = 'BlueRed3'

	plt.subplot(2, 2, 3)
	plt.imshow(Z, interpolation='nearest')
	plt.colorbar()
	plt.title("Alpha = 1")

	# Or as yet another variation, we can replace the rcParams
	# specification *before* the imshow with the following *after*
	# imshow.
	# This sets the new default *and* sets the colormap of the last
	# image-like item plotted via pyplot, if any.
	#

	plt.subplot(2, 2, 4)
	# Draw a line with low zorder so it will be behind the image.
	plt.plot([0, 10*np.pi], [0, 20*np.pi], color='c', lw=20, zorder=-1)

	plt.imshow(Z, interpolation='nearest')
	plt.colorbar()

	# Here it is: changing the colormap for the current image and its
	# colorbar after they have been plotted.
	plt.set_cmap('BlueRedAlpha')
	plt.title("Varying alpha")
	#

	plt.suptitle('Custom Blue-Red colormaps', fontsize=16)

	plt.show()
	exit()
	# BAAAAAAAAAAAAAAAAAAAAAAAAA



show_graphs = 1

pairtype = "rho"    # Cantor pairing function
forcedmax = False
forcedmin = False

default_fontsize = 22

SCALE = 10.0

import pylab as m

cdict = {
'red'  :  ((0.00, 1.00, 1.00), (0.50, 0.25, 0.25), (1.00, 1.00, 1.00)),
'green':  ((0.00, 1.00, 1.00), (0.70, 0.00, 0.50), (1.00, 1.00, 1.00)),
'blue' :  ((0.00, 1.00, 1.00), (0.50, 0.00, 0.00), (1.00, 1.00, 1.00))
}
cdict = {'red':   ((0.0,  0.0, 0.0),
                   (0.5,  1.0, 1.0),
                   (1.0,  1.0, 1.0)),

         'green': ((0.0,  0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.75, 1.0, 1.0),
                   (1.0,  1.0, 1.0)),

         'blue':  ((0.0,  0.0, 0.0),
                   (0.5,  0.0, 0.0),
                   (1.0,  1.0, 1.0))}

cdict = {'red':  ((0.0, 0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.5, 0.8, 1.0),
                   (0.75, 1.0, 1.0),
                   (1.0, 0.4, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.5, 0.9, 0.9),
                   (0.75, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.4),
                   (0.25, 1.0, 1.0),
                   (0.5, 1.0, 0.8),
                   (0.75, 0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }


helix_start = 0.31
helix_end   = 0.39
sheet_start = 0.45
sheet_end   = 0.62
polyproline_end = 0.66
cdict = {  
#                          white white                 white                                     white                 white                                                                                      white->
           'red': ((0.00,  1.00, 1.00), (helix_start,  1.00, helixR[0]), (helix_end,  helixR[0], 1.00), (sheet_start,  1.00, sheet[0]), (sheet_end,  sheet[0], polyproline[0]), (polyproline_end, polyproline[0], 1), (1, 1,1)), 
         'green': ((0.00,  1.00, 1.00), (helix_start,  1.00, helixR[1]), (helix_end,  helixR[1], 1.00), (sheet_start,  1.00, sheet[1]), (sheet_end,  sheet[1], polyproline[1]), (polyproline_end, polyproline[1], 1), (1, 1,1)),
          'blue': ((0.00,  1.00, 1.00), (helix_start,  1.00, helixR[2]), (helix_end,  helixR[2], 1.00), (sheet_start,  1.00, sheet[2]), (sheet_end,  sheet[2], polyproline[2]), (polyproline_end, polyproline[2], 1), (1, 1,1))  
        }


colors    = [(1,1,1), (1,1,1), helixR, (1,1,1), sheet, sheet, polyproline, (1,1,1), (1,1,1)]
positions = [ 0.00  , 0.30   ,  0.37 ,  0.44  ,  0.54,  0.57,      0.62  ,  0.70  ,  1.00  ]

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# c stands for color, bc stands for background color
c1 = [0,0,0] # black
c2 = [1,1,0] # yellow 
c3 = [1,0,0] # red
c4 = [0,0,1] # blue
bc = [1,1,1] # white

cdict3 = {
#                         black  black           white  white         yellow  red             white  white         blue   blue
	'red':   ((0.00,  c1[0], c1[0]), (0.25,  bc[0], bc[0]), (0.50, c2[0], c3[0]), (0.75,  bc[0], bc[0]), (1.0, c4[0], c4[0])), 
	'green': ((0.00,  c1[1], c1[1]), (0.25,  bc[1], bc[1]), (0.50, c2[1], c3[1]), (0.75,  bc[1], bc[1]), (1.0, c4[1], c4[1])),
	'blue':  ((0.00,  c1[2], c1[2]), (0.25,  bc[2], bc[2]), (0.50, c2[2], c3[2]), (0.75,  bc[2], bc[2]), (1.0, c4[2], c4[2])) 
}

# More dramatic version
cdict3 = {
#                         white  white           black  yellow         white  white           white  white         blue   blue
	'red':   ((0.00,  bc[0], bc[0]), (0.25,  c1[0], c2[0]), (0.50, bc[0], bc[0]), (0.75,  c3[0], c4[0]), (1.0, bc[0], bc[0])), 
	'green': ((0.00,  bc[1], bc[1]), (0.25,  c1[1], c2[1]), (0.50, bc[1], bc[1]), (0.75,  c3[1], c4[1]), (1.0, bc[1], bc[1])),
	'blue':  ((0.00,  bc[2], bc[2]), (0.25,  c1[2], c2[2]), (0.50, bc[2], bc[2]), (0.75,  c3[2], c4[2]), (1.0, bc[2], bc[2])) 
}
cmapTEST = LinearSegmentedColormap('BlackYellowRedBlue', cdict3)
#plt.register_cmap(cmap=cmapTEST)
#cmap = plt.get_cmap('BlackYellowRedBlue')


#generate the colormap with 1024 interpolated values
#my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('HelixSheet', cdict)

postscript_template = """%!PS-Adobe-3.0 EPSF-3.0
%%Document-Fonts: Times-Roman
%%Pages: 1
%%BoundingBox: 0 0 XMAX YMAX
%%LanguageLevel: 1
%%EndComments
%%BeginProlog
%%EndProlog
SCALE dup scale
% -------------
% Define text
/Times-Roman findfont
0.6 scalefont
setfont
% ----------------------------------------------------------------
% STRING STUFF
% -------------
% Find the heighth of the character string 
/FindHeight { 
gsave 
   (0) true charpath 
   flattenpath pathbbox 
   exch pop exch sub 
   /height exch def pop 
   grestore }
def 
% -------------
% Center the string vertically and horizontally 
/ccshow { 
   FindHeight 
   height 2 div neg 
   exch dup stringwidth pop 
   2 div neg 
   3 -1 roll rmoveto 
   show } 
def 
% -------------
% Center the string vertically and horizontally 
/ccshowrotate { 
   90 rotate
   FindHeight
   height 2 div neg 
   exch dup stringwidth pop 
   2 div neg 
   3 -1 roll rmoveto
   show 
   90 neg rotate
   } 
def 
% -------------
% Left justify and center vertically 
/lcshow { 
   FindHeight 
   0 
   height 2 div neg 
   rmoveto 
   show } 
def 
% -------------
% Right justify and center vertically 
/rcshow { 
   FindHeight 
   height 2 div neg 
   exch dup stringwidth pop neg 
   3 -1 roll rmoveto 
   show } 
def 
% -------------
% Center horizontally and set to top vertically 
/ctshow { 
   FindHeight 
   height neg 
   exch dup stringwidth pop 
   2 div neg 
   3 -1 roll rmoveto 
   show } 
def 
% -------------
% Left justify and set to top vertically 
/ltshow { 
   FindHeight 
   0 
   height neg 
   rmoveto 
   show } 
def 
% -------------
% Right justify and set to top vertically 
/rtshow { 
   FindHeight 
   height neg 
   exch dup stringwidth pop neg 
   3 -1 roll rmoveto 
   show } 
def 
% -------------
% Center horizontally and set to bottom vertically 
/cbshow { 
   0 
   exch dup stringwidth pop 
   2 div neg 
   3 -1 roll rmoveto 
   show } 
def 
% -------------
% Left justify and set to bottom vertically 
/lbshow { show } def 
% Right Justify and set to bottom vertically 
/rbshow { 
   0 
   exch dup stringwidth pop neg 
   3 -1 roll rmoveto 
   show } 
def 
% -------------
USERMATERIAL
% -------------
showpage
%%Trailer
"""

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

def ps_draw_shape(X,Y,linewidth=1,linecolor=(0,0,0),fill=1,fillcolor=(0.6,0.6,0.6)):
	""" 
	% This will create a red lined and blue filled triangle
	0 0 1 setrgbcolor
	100 100 moveto
	200 100 lineto
	150 150 lineto
	closepath
	gsave fill grestore
	1 0 0 setrgbcolor
	2 setlinewidth
	stroke
	"""
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
	
	if fill:
	#	ps_text+= "gsave\n"
		ps_text+= "%.3f %.3f %.3f setrgbcolor\n"%fillcolor
		ps_text+= "fill\n"
	#	ps_text+= "grestore\n"
	
	if linewidth:
		ps_text+= "%.3f %.3f %.3f setrgbcolor\n"%linecolor
		ps_text+= str(linewidth)+" setlinewidth\n"
	ps_text+= "stroke"
	
	return ps_text+"\n"


# To see how our EPS system works, set 0 to 1 here
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

'''
NAME
	Custom Colormaps for Matplotlib
PURPOSE
	This program shows how to implement make_cmap which is a function that
	generates a colorbar.  If you want to look at different color schemes,
	check out https://kuler.adobe.com/create.
PROGRAMMER(S)
	Chris Slocum
REVISION HISTORY
	20130411 -- Initial version created
	20140313 -- Small changes made and code posted online
	20140320 -- Added the ability to set the position of each color
'''
def make_cmap(colors, position=None, bit=False):
	'''
	make_cmap takes a list of tuples which contain RGB values. The RGB
	values may either be in 8-bit [0 to 255] (in which bit must be set to
	True when called) or arithmetic [0 to 1] (default). make_cmap returns
	a cmap with equally spaced colors.
	Arrange your tuples so that the first color is the lowest value for the
	colorbar and the last is the highest.
	position contains values from 0 to 1 to dictate the location of each color.
	'''
	import matplotlib as mpl
	bit_rgb = np.linspace(0,1,256)
	if position == None:
		position = np.linspace(0,1,len(colors))
	else:
		if len(position) != len(colors):
			sys.exit("position length must be the same as colors")
		elif position[0] != 0 or position[-1] != 1:
			sys.exit("position must start with 0 and end with 1")
	if bit:
		for i in range(len(colors)):
			colors[i] = (bit_rgb[colors[i][0]],
						 bit_rgb[colors[i][1]],
						 bit_rgb[colors[i][2]])
	cdict = {'red':[], 'green':[], 'blue':[]}
	for pos, color in zip(position, colors):
		cdict['red'].append((pos, color[0], color[0]))
		cdict['green'].append((pos, color[1], color[1]))
		cdict['blue'].append((pos, color[2], color[2]))

	cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
	return cmap

### An example of how to use make_cmap
fig = plt.figure()
ax = fig.add_subplot(311)
### Create a list of RGB tuples
colors    = [(1,1,1), (1,1,1), helixR, (1,1,1), sheet, sheet, polyproline, (1,1,1), (1,1,1)]
positions = [ 0.00  , 0.30   ,  0.37 ,  0.44  ,  0.54,  0.57,      0.62  ,  0.70  ,  1.00  ]
#                               helix             sheet    loop
### Call the function make_cmap which returns your colormap
#secondarystructure_cmap = make_cmap(colors, position=positions)
secondarystructure_cmap = my_cmap
#secondarystructure_cmap = plt.get_cmap("seismic_r")
#secondarystructure_cmap = plt.get_cmap("Accent")
### Use your colormap

figure_base="./"+str(pairtype)#+"figs"
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
	"-":"Undesignated"
}

#cmap = 
#cmap = plt.get_cmap("YlOrRd")
#cmap = plt.cm.Blues
#code_to_cmap = plt.cm.Reds
#code_to_cmap = plt.cm.Oranges
#code_to_cmap = plt.cm.Greens
#code_to_cmap = plt.cm.Purples
#code_to_cmap = plt.cm.seismic
#plt.cm.bone_r
#plt.cm.gray_r
#cm = plt.cm.gray_r
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
	"-":cm    #"Undesignated"
}


def make2Dfigure(Xoriginal,Yoriginal,Zoriginal,fn=0,xlim=[],ylim=[],zlim=[],cmap=plt.get_cmap('gray_r'),yscaling=0.5,xtitle="",ytitle="",xticks=[],xlabels=[],yticks=[],ylabels=[],horizontallines=[],verticallines=[],title="",showeps=0):
	if fn == 0 or str(fn)[-len(".eps"):] !=".eps":
		print "No EPS filename given. Exiting."
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
	
	# Working on setting up the dimensions of the grawing's frame 
	# X
	xmin = Xsorted[0]
	xmax = Xsorted[-1]
	pageMinX = 3.0  # Arbitrary values
	pageMaxX = 16.0 # Arbitrary values
	
	# Y
	ymin = Ysorted[0]
	ymax = Ysorted[-1]
	pageMinY = 2.5
	pageMaxY = pageMinY + (pageMaxX-pageMinX)*yscaling
	
	if len(Xsorted) == 1:
		pageMaxX = 5.5
	
	
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
		xmin = Xsorted[0]
		xmax = Xsorted[-1]
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
		ymin = Ysorted[0]
		ymax = Ysorted[-1]
		
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
		if len(xlabels) == 0:
			xmin_for_extended = Xsorted[0]
			xmax_for_extended = Xsorted[-1] #ax.xaxis.get_data_interval()[-1]
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
						print "Xlabels are neither matched by 'xticks' nor by the number of x values. Exiting."
						exit()
		for x in xticks:
			#print [pageMinX, pageSpacingX, pageMaxX, x, xmin, xmax]
			pagex = pageMinX + pageSpacingX*0.5 + float(pageMaxX - pageMinX - pageSpacingX)*(x - xmin)/(xmax-xmin)
			page_xticks.append(pagex)
		
			
	#
	if len(set(Y)) > 1:
		if len(ylabels) == 0:
			ymin_for_extended = Ysorted[0]
			ymax_for_extended = Ysorted[-1] #ax.yaxis.get_data_interval()[-1]
			yrange = ymax_for_extended-ymin_for_extended
			ymin_for_extended, ymax_for_extended = (ymin_for_extended - yrange * 0.05, ymax_for_extended + yrange * 0.05)		
			ylabels = Extended(ymin_for_extended,ymax_for_extended,density_val=1.2)
			yticks  = copy.deepcopy(ylabels)			
		else:
			if len(ylabels) == len(yticks):
				# First we check if the user provided yticks with the ylabels
				pass
			elif len(ylabels) == len(Ysorted):
				# If this is the case, then we eypect the two to be 
				yticks = copy.deepcopy(Ysorted)
			else:
				yticks = []
				for y in ylabels:
					if isinstance(y,float) or isinstance(y,int):
						yticks.append(y)
					else:
						print "Ylabels are neither matched by 'yticks' nor by the number of y values. Eyiting."
						eyit()
		for y in yticks:
			#print [pageMinY, pageSpacingY, pageMaxY, y, ymin, ymax]
			pagey = pageMinY + pageSpacingY*0.5 + float(pageMaxY - pageMinY - pageSpacingY)*(y - ymin)/(ymax-ymin)
			page_yticks.append(pagey)
	
	usermaterial = ""
	
	
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
		heat = cmap(float(z-zmin)/float(zmax-zmin))
		
		usermaterial += ps_draw_shape(Xs,Ys,linewidth=0,linecolor=0,fillcolor=heat[:3])
	
	# Drawing the frame
	frame_thickness = 0.05
	xmin = round(pageMinX - frame_thickness/2.0,3)
	xmax = round(pageMaxX + frame_thickness/2.0,3)
	ymin = round(pageMinY - frame_thickness/2.0,3)
	ymax = round(pageMaxY + frame_thickness/2.0,3)
	Xs = [xmin,xmax,xmax,xmin]
	Ys = [ymin,ymin,ymax,ymax]
	usermaterial += ps_draw_shape(Xs,Ys,linewidth=frame_thickness,linecolor=0,fill=0,fillcolor=1)
	
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
		# Adding the x label!
		textx = float(pageMinX)/3.0
		texty = pageMinY + float(pageMaxY-pageMinY)/2.0
		usermaterial += "%.10f %.10f moveto\n"%(textx,texty)
		#usermaterial += "90 rotate\n"
		usermaterial += "("+str(ytitle)+")  ccshowrotate\n"
		#usermaterial += "90 neg rotate\n"
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
		
	print "WRITTEN TO:",fn
	
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
		plt.axes().set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))*yscaling)
		plt.show()
	
	return 1

# KEYS FOR DSSP:
# G --> 3-turn helix
# H --> 4-turn helix (alpha helix)
# I --> 5-turn helix
# E --> strand
"""
G = 3-turn helix (3_10-helix). Min length 3 residues.
H = 4-turn helix (alpha-helix). Min length 4 residues.
I = 5-turn helix (pi-helix). Min length 5 residues.
T = hydrogen bonded turn (3, 4 or 5 turn)
E = extended strand in parallel and/or anti-parallel beta-sheet conformation. Min length 2 residues.
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
rho_scaling = 10.0 # This is sigma in the manuscript (Mannige, Kundu, Whitelam, 2016)
multiplier = int(round(rho_scaling*(bound*(2.0**0.5)),0)) # For internal reference
multiplier_by_two = round(rho_scaling*(bound*(2.0**0.5)),0)/2.0 # For internal reference

# To get the raw Ramachandran number from phi and psi:
def raw_ramachandran_number_collapse(phi,psi):
	phi = float(phi)
	psi = float(psi)
	a = round(rho_scaling*(phi-psi + bound)/math.sqrt(2.0),0)
	b = round(rho_scaling*(phi+psi + bound)/math.sqrt(2.0),0)
	return a + b*multiplier
#
# Original phi and psi values are easily obtained from the *raw* Ramachandran number
def raw_ramachandran_number_expand(z):
	z = float(z)
	phi = (math.sqrt(2.0)*np.mod(z,multiplier)/rho_scaling     - bound + math.sqrt(2.0)*np.floor(z / multiplier)/rho_scaling - bound )/2.0
	psi = (math.sqrt(2.0)*np.floor(z / multiplier)/rho_scaling - bound - math.sqrt(2.0)*np.mod(z, multiplier)/rho_scaling + bound)/2.0
	return phi,psi

# First getting the lowest and highest possible unnormalized R numbers
raw_R_min = raw_ramachandran_number_collapse(-180,-180)
raw_R_max = raw_ramachandran_number_collapse(180,180)
# To get the normalized Ramachandran number from the raw Ramachandran number 
def normalized_ramachandran_number(phi,psi):
	raw_R = raw_ramachandran_number_collapse(phi,psi)
	R = float(raw_R - raw_R_min)/(raw_R_max-raw_R_min)
	return R

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
def normalized_ramachandran_number(phi,psi):
	# with a twist!
	phi = float(phi)
	psi = float(psi)
	a = round(rho_scaling*(phi-psi + bound)/math.sqrt(2.0),0)
	b = round(rho_scaling*(phi+psi + bound)/math.sqrt(2.0),0)
	
	raw_r = a + b*multiplier
	final_r = float(raw_r - raw_R_min)/float(raw_R_max - raw_R_min)
	if a >= multiplier_by_two:
		final_r = final_r * -1.0
	return final_r

#To the Ramachandran number function, set "if 0:" to "if 1:"
if 0:
	current_phi = -60
	current_psi = -60
	# To get the normalized Ramachandran number:
	R = normalized_ramachandran_number(current_phi,current_psi)
	print "R(phi="+str(current_phi)+",psi="+str(current_psi)+") = "+str(R)

#
def histogram2d(X_vals,Y_vals,cmap=plt.cm.Blues,xyrange=[],title="",fn=""):
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

def read_pdb(pdbblock):
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
	
	getlines = re.compile("ATOM\s+(?P<atomno>\d+)\s+(?P<atomtype>\S+)\s+.(?P<resname>...)..\s+(?P<resno>\d+)\s+(?P<x>\-*\d+\.*\d*)\s+(?P<y>\-*\d+\.*\d*)\s+(?P<z>\-*\d+\.*\d*).{17}(?P<segname>.{5})",re.M)
	getlines_short = re.compile("ATOM\s+(?P<atomno>\d+)\s+(?P<atomtype>\S+)\s+(?P<resname>...).(?P<segname>.)\s+(?P<resno>\d+)\s+(?P<x>\-*\d+\.*\d*)\s+(?P<y>\-*\d+\.*\d*)\s+(?P<z>\-*\d+\.*\d*)",re.M)
	
	resnos = []
	resno_to_coordN  = {}
	resno_to_coordC  = {}
	resno_to_coordCA = {}
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
						rho = normalized_ramachandran_number(phi,psi)
						#print rho
						#if rho < 0.5:
						#	print (phi,psi)
						
						model_to_chain_to_resno_atom_to_vals[model][chain][i]["phi"] = phi
						model_to_chain_to_resno_atom_to_vals[model][chain][i]["psi"] = psi
						model_to_chain_to_resno_atom_to_vals[model][chain][i]["r"] = rho
	
	return model_to_chain_to_resno_atom_to_vals

# Draw the mapping between phi,psi and Ramachandran number
def show_ramachandran_mapping(cmap=plt.get_cmap("Paired"),stepsize=10):	
	PHI = []
	PSI = []
	RHO = []
	rho_max = raw_ramachandran_number_collapse(180.0,180.0)
	rho_min = raw_ramachandran_number_collapse(-180.0,-180.0)
	for phi in range(0,361,stepsize):
		for psi in range(0,361,stepsize):
			PHI.append(phi)
			PSI.append(psi)
			unnormalized_rho = raw_ramachandran_number_collapse(phi,psi)
			RHO.append(float(unnormalized_rho - rho_min)/(rho_max - rho_min))
	x=np.array(PHI)
	y=np.array(PSI)
	z=np.array(RHO)
	
	color_bar_range = np.arange(0,1.01,0.01)
	
	fn = "rho"
	fn_bar = fn+"_colorbar.eps"
	print "#WRITING TO:",fn_bar
	make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=cmap, fn=fn_bar, yscaling=0.5,xtitle="Key",ytitle="R",showeps=showeps)
	#
	print "#WRITING TO:",fn+".eps"
	make2Dfigure(x,y,z,fn+".eps",cmap=cmap,yscaling=1,xtitle="phi",ytitle="psi",xlabels=range(-360,361,90),ylabels=range(-360,361,90),showeps=showeps)
	#
	return 1

# To see how the (phi,psi) to Ramachandran number mapping occurs, uncomment this:
#show_ramachandran_mapping()

if __name__ == "__main__":
	if not "-pdb" in sys.argv:
		print "Must provide '-pdb' parameter. Exiting."
		exit(0)
	for i in range(len(sys.argv)):
		if sys.argv[i] == "-rmsd":
			showrmsd = 1
		if sys.argv[i] == "-pdb":
			if len(sys.argv) <= i+1:
				print "MUST PROVIDE PDB NAME."
				exit(0)
			else:
				pdbfn = str(sys.argv[i+1])
				print "# pdbfn set to:",pdbfn
		elif sys.argv[i] == "-bins":
			if len(sys.argv) <= i+1:
				print "When using '-bins', you must provide bin number. Exiting."
				exit(0)
			else:
				if not sys.argv[i+1].isdigit():
					print "The -bin parameter must be a positive integer (provided: "+str(sys.argv[i+1])+") Exiting."
					exit(0)
				else:
					bins = int(sys.argv[i+1])
					print "# bins set to:",bins
					if bins == 0:
						print "Must have greater than 0 bins. Exiting."
						exit(0)
	
	#
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
		exit("Either filename or directory expected. Exiting.")
	
	target_dir = pdbdir+"/"+os.path.basename(pdbfilenames[0])[:-len(".pdb")]+"/"
	if len(pdbfilenames)>1:
		target_dir = pdbdir+"/report/"
	if not os.path.isdir(target_dir):
		os.makedirs(target_dir)

	target_base = target_dir.rstrip("/")

	# JUST "CLEVERLY" ARRANGING THE FILENAMES 
	# (e.g., pdbfilenames = [something2part1,something1part2,something1part1,something10part1]
	# pdbfilenames.sort() this list to: [something1part1,something1part2,something2part1,something10part1]
	REXP = re.compile( r'\d+' )
	def key_function( s ): return map(int, re.findall(REXP, s ))
	pdbfilenames.sort( key=key_function)

	pdbblock = ""
	for pdbfn in pdbfilenames:
		print "# Reading '"+pdbfn+"'"
		f = open(pdbfn,"r")
		pdbblock +="\n\nEND\n\n"+f.read()
		f.close()

	offset = 0
	# structure -> model -> chain -> residue -> atom -> 'x','y','z','phi','psi','r'
	print "# Parsing the PDB (structure) data"
	structure = read_pdb(pdbblock)
	print "\t...done"
	
	# FIRST, GETTING THE DIFFERENT CHAIN IDs
	chains = []
	for model in sorted(structure.keys()):
		for chain in sorted(structure[model].keys()):
			chains.append(chain)
	chains = list(sorted(set(chains)))

	# FINALLY, WE WILL GO THROUGH EACH CHAIN AND PRODUCE GRAPHS
	for chain in chains:
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
			
			#fn = target_base+".rcode.ss.eps"
			#cbfn = fn+".colorbar.eps"
			#color_bar_range = np.arange(round(sortedZ[0],0),round(sortedZ[-1],0)+0.005,0.01)
			#make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=secondarystructure_cmap, fn=cbfn, ytitle="R",showeps=0)#,ylim = [0.3,0.7])#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
			#make2Dfigure(X,Y,Z,fn, cmap=secondarystructure_cmap, xtitle="Model #",ytitle="Residue #",showeps=showeps)#,ylim = [0.3,0.7])#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
			
			fn = target_base+".rcode.eps"
			cbfn = fn+".colorbar.eps"
			color_bar_range = np.arange(round(sortedZ[0],0),round(sortedZ[-1],0)+0.005,0.01)
			cmap = cmapTEST #plt.get_cmap("seismic_r")
			make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=cmap, fn=cbfn, ytitle="R",showeps=0)#,ylim = [0.3,0.7])#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
			make2Dfigure(X,Y,Z,fn, cmap=cmap, xtitle="Model #",ytitle="Residue #",showeps=showeps)#,ylim = [0.3,0.7])#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
		
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
			
			fn = target_base+".rcode.rmsd"+str(reference_index)+".eps"
			cbfn = fn+".colorbar.eps"
			color_bar_range = np.arange(0,1.01,0.01)
			make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=plt.get_cmap("Reds"), ytitle="RMSD (R)", fn=cbfn, yscaling=0.5,showeps=0)#,ylim = [0.3,0.7])#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
			make2Dfigure(X,Y,Z,fn, cmap=plt.get_cmap("Reds"), xtitle="Model #",ytitle="RMSD (R)",showeps=showeps)#, cmap=plt.get_cmap("gray"))#,ylim = [0.3,0.7])#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
		
		# ----------------------------------
		# Fluctuations from previous time
		if showrmsf:
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
			
			fn = target_base+".rcode.rmsf.eps"
			cbfn = fn+".colorbar.eps"
			color_bar_range = np.arange(0,1.01,0.01)
			
			make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=plt.get_cmap("Blues"), ytitle="RMSF (R)", fn=cbfn, yscaling=0.5,showeps=0)#,ylim = [0.3,0.7])#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
			make2Dfigure(X,Y,Z,fn, cmap=plt.get_cmap("Blues"), xtitle="Model #",ytitle="RMSF (R)",showeps=showeps)#, cmap=plt.get_cmap("gray"))#,ylim = [0.3,0.7])#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
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
					a,b = np.histogram(Rs,bins=np.arange(-1.01,1.03,0.02))
					for i in range(len(a)):
						X.append(model)
						Y.append(float(b[i]+b[i+1])/2.0)
						Z.append(float(a[i])/np.max(a))
			#
			fn = target_base+".rcode.his.eps"
			cbfn = fn+".colorbar.eps"
			
			color_bar_range = np.arange(0,1.01,0.01)
			make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=plt.get_cmap("gray_r"), ytitle="P(R)/max(P(R))", fn=cbfn, yscaling=0.5,showeps=0)#,ylim = [0.3,0.7])#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
			
			for i in range(dofilter):
				Z = median_filter(Z)
			
			demeY = []
			for x,y in zip(X,Y):
				if x == X[0]:
					demeY.append(y)
			print sorted(demeY)
			make2Dfigure(X,Y,Z,fn, cmap=plt.get_cmap("gray_r"), xtitle="Model #",ytitle="R",showeps=showeps)#,ylim = [0.3,0.7])#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
	
	
	if do_vmd_etc:
		# writing first and last 
		ms = pdbblock.split("\nEND\n")#|TER)\n",pdbblock,re.MULTILINE)
		pdbmodels = []
		for pdbmodel in ms:
			if len(pdbmodel.rstrip()) > 1:
				pdbmodels.append(pdbmodel)

		lines_par_preamble = """with line
		line on
		line loctype world
		line g0
		line XSTART, 0.3, XSTOP, 0.3
		line linewidth 4.0
		line linestyle 3
		line color 10
		line arrow 0
		line arrow type 0
		line arrow length 1.000000
		line arrow layout 1.000000, 1.000000
		line def
		with line
		line on
		line loctype world
		line g0
		line XSTART, 0.4, XSTOP, 0.4
		line linewidth 4.0
		line linestyle 3
		line color 10
		line arrow 0
		line arrow type 0
		line arrow length 1.000000
		line arrow layout 1.000000, 1.000000
		line def
		with line
		line on
		line loctype world
		line g0
		line XSTART, 0.47, XSTOP, 0.47
		line linewidth 4.0
		line linestyle 3
		line color 5
		line arrow 0
		line arrow type 0
		line arrow length 1.000000
		line arrow layout 1.000000, 1.000000
		line def
		with line
		line on
		line loctype world
		line g0
		line XSTART, 0.58, XSTOP, 0.58
		line linewidth 4.0
		line linestyle 3
		line color 5
		line arrow 0
		line arrow type 0
		line arrow length 1.000000
		line arrow layout 1.000000, 1.000000
		line def
		line def
		with line
		line on
		line loctype world
		line g0
		line XSTART, 0.2, XSTOP, 0.2
		line linewidth 2
		line linestyle 1
		line color 1
		line arrow 0
		line arrow type 0
		line arrow length 1.000000
		line arrow layout 1.000000, 1.000000
		line def
		"""

		sheet_par_template = """with line
		line on
		line loctype world
		line g0
		line XSTART, YSTART, XSTOP, YSTOP
		line linewidth 20.0
		line linestyle 1
		line color 5
		line arrow 2
		line arrow type 0
		line arrow length 1.000000
		line arrow layout 1.000000, 1.000000
		line def
		"""

		turn_par_template = """with line
		line on
		line loctype world
		line g0
		line XSTART, YSTART, XSTOP, YSTOP
		line linewidth 7.0
		line linestyle 1
		line color 1
		line arrow 0
		line arrow type 0
		line arrow length 1.000000
		line arrow layout 1.000000, 1.000000
		line def
		"""

		helix_par_template = """with line
		line on
		line loctype world
		line g0
		line XSTART, YSTART, XSTOP, YSTOP
		line linewidth 20.0
		line linestyle 1
		line color 10
		line arrow 0
		line arrow type 0
		line arrow length 1.000000
		line arrow layout 1.000000, 1.000000
		line def
		"""

		secondary_structure_code_to_parameter_template = {
			# helices
			"G":helix_par_template,
			"H":helix_par_template,
			"I":helix_par_template,
			# sheet
			"T":turn_par_template,
			# sheet
			"E":sheet_par_template,
			"B":sheet_par_template
		}
		
		for model_key in [sorted(structure.keys())[0],sorted(structure.keys())[-1]]:
			start_end = "first"
			if model_key != 1:
				start_end = "last"
			
			chain_to_resno_to_r = {}
			chain_to_resno_to_phipsi = {}
			for chain in structure[model_key].keys():
				chain_to_resno_to_r[chain] = {}
				chain_to_resno_to_phipsi[chain] = {}
				for resno in structure[model_key][chain].keys():
					if "r" in structure[model_key][chain][resno]:
						chain_to_resno_to_r[chain][resno] = structure[model_key][chain][resno]["r"]
						chain_to_resno_to_phipsi[chain][resno] = (structure[model_key][chain][resno]["phi"],structure[model_key][chain][resno]["psi"])
			
			# -------------------->
			pdbmodel = pdbmodels[model_key-1]
			getlines       = re.compile("ATOM\s+(?P<atomno>\d+)\s+(?P<atomtype>\S+)\s+.{5}\s+(?P<resno>\d+)\s+(?P<x>\-*\d+\.*\d*)\s+(?P<y>\-*\d+\.*\d*)\s+(?P<z>\-*\d+\.*\d*).{17}(?P<segname>.{5})",re.M)
			getlines_short = re.compile("ATOM\s+(?P<atomno>\d+)\s+(?P<atomtype>\S+)\s+.{5}\s+(?P<resno>\d+)\s+(?P<x>\-*\d+\.*\d*)\s+(?P<y>\-*\d+\.*\d*)\s+(?P<z>\-*\d+\.*\d*)",re.M)
			if not getlines.search(pdbmodel):
				# then we have a PDB that does not have a segname, which DSSP chokes on
				new_pdb_model = ""
				for l in pdbmodel.split("\n"):
					if getlines_short.match(l):
						new_pdb_model += l+"                 A    \n"
					else:
						new_pdb_model += l+"\n"
				pdbmodel = new_pdb_model
			
			tmp_pdb_fn = target_base+"_pdb_"+str(start_end)+".pdb"
			
			print "#WRITING TO:",tmp_pdb_fn
			
			f = open(tmp_pdb_fn,"w")
			f.write(pdbmodel)
			f.close()
			
			p = tmp_pdb_fn+".dat.png"
			if not os.path.isfile(p): # then open structure file using vmd and draw snapshot
				os.system("vmd -e templates/template.vmd -args "+tmp_pdb_fn)
				# The output PNG from the previous command will be
				os.system("convert "+p+" -trim "+p)
				"tachyon -res 1338 1668 -aasamples 12 class_h_1BKV_pdb_first.pdb.dat -format PNG -o class_h_1BKV_pdb_first.pdb.dat.png"
				
			# G --> 3-turn helix
			# H --> 4-turn helix (alpha helix)
			# I --> 5-turn helix
			# E --> strand
			resno_to_dssp = get_resid_to_dssp(tmp_pdb_fn)
			# <--------------------
			
			
			# WRITING SECONDARY STRUCTURES
			secondary_y_position = 0.2
			# H (alpha-helix)
			
			parfile = lines_par_preamble.replace("XSTART",str(np.min(resno_to_dssp.keys()))).replace("XSTOP",str(np.max(resno_to_dssp.keys())))
			for motifs_to_find in [["E","B"],["G","H","I"],["T"]]:
				resnos = []
				for resno in resno_to_dssp.keys():
					if resno_to_dssp[resno] in motifs_to_find:
						resnos.append(resno)
				
				if len(resnos):
					previous_stretch = resnos[0]
					startstop = []
					for i in range(len(resnos)):
						if i==0:
							startstop.append(resnos[i])
						elif resnos[i]-resnos[i-1] != 1:
							startstop+=[resnos[i-1],resnos[i]]
							if i == len(resnos)-1:
								startstop.append(resnos[i])
						elif i == len(resnos)-1:
							startstop.append(resnos[i])
					startstop_pairs = []
					for i in range(len(startstop)/2):
						startstop_pairs.append((startstop[i*2],startstop[i*2+1]))
					
					for xstart,xstop in startstop_pairs:
						current_par = secondary_structure_code_to_parameter_template[resno_to_dssp[xstart]]
						current_par = current_par.replace("XSTART",str(xstart))
						current_par = current_par.replace("YSTART",str(secondary_y_position))
						current_par = current_par.replace("XSTOP",str(xstop))
						current_par = current_par.replace("YSTOP",str(secondary_y_position))
						parfile += current_par
			
			current_par_fn = target_base+"_resno_vs_r_model_"+str(start_end)+".par"
			f = open(current_par_fn,"w")
			f.write(parfile)
			f.close()
			
			fn = target_base+"_resno_vs_r_model_"+str(start_end)+".dat"
			print "#Writing:",fn
			f = open(fn,"w")
			for chain in sorted(chain_to_resno_to_r.keys()):
				f.write("#CHAIN"+chain+"\n")
				for resno in chain_to_resno_to_r[chain].keys():
					f.write(str(resno)+"\t"+str(chain_to_resno_to_r[chain][resno])+"\n")
				f.write("\n\n")
			f.close()
			agrfile = fn[:-len(".dat")]+".agr"
			psfile  = fn[:-len(".dat")]+".ps"
			epsfile  = fn[:-len(".dat")]+".eps"
			
			#print sorted(resno_to_dssp.keys()),np.min(resno_to_dssp.keys())
			os.system("xmgrace -par templates/residue_vs_r.par "+fn+" -pexec 'title "+'""'+"' -par "+current_par_fn+"  -world "+str(np.min(resno_to_dssp.keys()))+" 0 "+str(np.max(resno_to_dssp.keys()))+" 1  -hardcopy -saveall "+agrfile+" -printfile "+psfile)
			#os.system("xmgrace -world ["+str(np.min(resno_to_dssp.keys()))+" 0 "+str(np.max(resno_to_dssp.keys()))+" 1]   "+agrfile+" -hardcopy -saveall "+agrfile+" -printfile "+psfile)
			
			#os.system("xmgrace "+fn+" -par "+current_par_fn+" -hardcopy -saveall "+agrfile+" -printfile "+psfile)
			os.system("ps2eps -f --rotate=+ "+psfile)
			#os.system("xmgrace "+agrfile)
			#raw_input()
			
			#RAWR = []
			RAWPHI = []
			RAWPSI = []
			
			fn = target_base+"_ramachandran_"+str(start_end)+".dat"
			print "#Writing:",fn
			f = open(fn,"w")
			for chain in sorted(chain_to_resno_to_phipsi.keys()):
				f.write("#CHAIN"+chain+"\n")
				for resno in chain_to_resno_to_phipsi[chain].keys():
					f.write(str(chain_to_resno_to_phipsi[chain][resno][0])+"\t"+str(chain_to_resno_to_phipsi[chain][resno][1])+"\n")
					RAWPHI.append(chain_to_resno_to_phipsi[chain][resno][0])
					RAWPSI.append(chain_to_resno_to_phipsi[chain][resno][1])
				f.write("\n\n")
			f.close()
			agrfile = fn[:-len(".dat")]+".agr"
			psfile  = fn[:-len(".dat")]+".ps"
			epsfile  = fn[:-len(".dat")]+".eps"
			os.system("xmgrace "+fn+" -pexec 'title "+'""'+"' -par templates/ramachandran.par -hardcopy -saveall "+agrfile+" -printfile "+psfile)
			#os.system("xmgrace "+fn+" -par "+current_par_fn+" -hardcopy -saveall "+agrfile+" -printfile "+psfile)
			os.system("ps2eps -f --rotate=+ "+psfile)
			#os.system("xmgrace "+agrfile)
			
			#
			if 1: # Draw individual ramachandran plots
				for chain in structure[model_key].keys():
					chain_to_resno_to_r[chain] = {}
					chain_to_resno_to_phipsi[chain] = {}
					for resno in structure[model_key][chain].keys():
						if "r" in structure[model_key][chain][resno]:
							chain_to_resno_to_r[chain][resno] = structure[model_key][chain][resno]["r"]
							chain_to_resno_to_phipsi[chain][resno] = (structure[model_key][chain][resno]["phi"],structure[model_key][chain][resno]["psi"])
				
				
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
				
				fn = fn[:-len(".dat")]+"_histogram.eps"
				make2Dfigure(PHI,PSI,WEIGHT,fn=fn,yscaling=1,xtitle="phi",ytitle="psi",xlabels=range(-180,181,90),ylabels=range(-180,181,90),showeps=showeps)
				#raw_input()
				
				X=[]
				Y=[]
				Z=[]
				for chain in structure[model_key].keys():
					chain_to_resno_to_r[chain] = {}
					chain_to_resno_to_phipsi[chain] = {}
					for resno in structure[model_key][chain].keys():
						if "r" in structure[model_key][chain][resno]:
							X.append(1)
							Y.append(resno)
							Z.append(structure[model_key][chain][resno]["r"])
				
				color_bar_range = np.arange(0,1.01,0.01)
				fn = fn[:-len(".dat")]+"_r1_colorbar.eps"
				print "#WRITING TO:",fn
				make2Dfigure(np.ones(len(color_bar_range)),color_bar_range,color_bar_range,cmap=secondarystructure_cmap, fn=fn, ytitle="R",showeps=0)#,ylim = [0.3,0.7])#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
				fn = fn[:-len(".dat")]+"_r1.eps"
				print "#WRITING TO:",fn
				
				#make2Dfigure(X,Y,Z,fn, cmap=secondarystructure_cmap, xtitle="Model #",ytitle="Residue #",showeps=1)#,ylim = [0.3,0.7])#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
				for i in range(dofilter):
					Z = median_filter(Z)
				make2Dfigure(X,Y,Z,fn, cmap=secondarystructure_cmap, xtitle="Model #",ytitle="Residue #",showeps=showeps)#,ylim = [0.3,0.7])#,zlim=[0.0,0.3])# ylim = [0.3,0.7])
				
			#raw_input()
			
		exit()
		f = open("templates/latex_template.tex","r")
		latexblock = f.read()
		f.close()

		basename = os.path.basename(pdbfn)[:-len(".pdb")]
		texoutput = os.path.dirname(pdbfn)+"/"+basename+".tex"
		pdfoutput = os.path.dirname(pdbfn)+"/"+basename+".pdf"

		f = open("tmp.tex","w")
		f.write(latexblock.replace("FILEBASE",basename))
		f.close()

		os.system("pdflatex tmp.tex")
		os.system("cp tmp.tex "+texoutput)
		os.system("cp tmp.pdf "+pdfoutput)

