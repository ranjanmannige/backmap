"""A helper file that describes colormaps that are needed to create graphs

These colormaps are used either by either by :py:func:`backmap.backmap.draw_figures` 
or the stand along app run by :py:mod:`backmap.cli`
"""

# matplotlib imports
import matplotlib.pyplot as plt
# A function to create custom color maps for matplotb plots
from matplotlib.colors import LinearSegmentedColormap
# For displaying the cmaps using display_cmaps()
import matplotlib.cm as cm
import matplotlib.colors as colors

import warnings
#warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning)
#warnings.filterwarnings("ignore", category=DeprecationWarning)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SETTING UP SOME COLORMAPS

COLORSWITCH = 0.5  # THIS IS THE POINT, FOR THE RED/BLUE AND RED/BLUE/YELLOW/BLACK 
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
helixR      = (1.,0.,0.)
sheet       = (0.,0.,1.)
polyproline = (0.,1.,1.)


cmaps = {}
"""Mapping of custom colormap names to LinearSegmentedColormap instances. 
Reference the :ref:`Publications <publications>` section of the manual to 
see how the colors are used to understand/represent the properties a protein 
backbone (often, over an ensemble or time). Also, to visualize these cmaps,
refer to the :py:func:`backmap.local_colormaps.display_cmaps` function.

Keys:
	:Chirality: A red-blue coloring method, where red is predominantly 
				right handed and blue is predominantly left handed. 
				Both colors fade to white as they leave the diagonal.
	:Chirality_r: Similar to ``Chirality``, but with the color-to-white 
				gradient reversed (i.e., the diagonal is white, while 
				the bottom left and top right portions are red and blue 
				respectively).
	:ChiralityFourColor: Similar to the properties above, but with two 
				more colors that are needed to represent the negative values.
	:ChiralityFourColor_r: Same as the other above and, but with inverted 
				color-to-white gradients.
	:[SS]Hard: [SS=SecondaryStructure] The helix, sheet, and ppII helices shown as 
				red, blue and cyan. Boundaries are hard (i.e., without gradients)
	:SecondaryStructure: Same as above, but with color-to-white gradients (so, a backbone with a darker 
						 red color is closer to the average helice's position than a lighter red) 
	:[SS]FourColor: [SS=SecondaryStructure] Same as ``SecondaryStructureHard``, but accomodating negative R numbers.
"""

cmap_ranges = {}
"""Listing of custom colormap names to their expected ranges. Not using these ranges will result in inaccurate coloring."""

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
#cmap = LinearSegmentedColormap('Chirality', cdict)
#register_cmap(cmap=cmap)
cmaps['Chirality'] = LinearSegmentedColormap('Chirality', cdict)
cmap_ranges['Chirality'] = [0,1]
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
#cmap = LinearSegmentedColormap('Chirality_r', cdict)
#register_cmap(cmap=cmap)
cmaps['Chirality_r'] = LinearSegmentedColormap('Chirality_r', cdict)
cmap_ranges['Chirality_r'] = [0,1]
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
#cmap = LinearSegmentedColormap('ChiralityFourColor', cdict)
#register_cmap(cmap=cmap)
cmaps['ChiralityFourColor'] = LinearSegmentedColormap('ChiralityFourColor', cdict)
cmap_ranges['ChiralityFourColor'] = [-1,1]
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
#cmap = LinearSegmentedColormap('Chirality_rFourColor', cdict)
#register_cmap(cmap=cmap)
cmaps['Chirality_rFourColor'] = LinearSegmentedColormap('Chirality_rFourColor', cdict)
cmap_ranges['Chirality_rFourColor'] = [-1,1]
# -------------------------
# NEW COLOR SCHEME: secondary structure (expected range: R=[0,1])
# ----------------
# This lets you later on get the cmap by name 'SecondaryStructure': cmap = plt.get_cmap('SecondaryStructure')
#
#                         white  white                  white        red                      red  white                  white      blue                    blue            cyan                               cyan              white white
cdict = {  #                  |      |                      |          |                        |      |                      |         |                       |               |                                  |                  |     |
           'red': ((0.00,  bc[0], bc[0]), (helix_start,  bc[0], helixR[0]), (helix_end,  helixR[0], bc[0]), (sheet_start,  bc[0], sheet[0]), (sheet_end,  sheet[0], polyproline[0]), (polyproline_end, polyproline[0], bc[0]), (1, bc[0],bc[0])), 
         'green': ((0.00,  bc[1], bc[1]), (helix_start,  bc[1], helixR[1]), (helix_end,  helixR[1], bc[1]), (sheet_start,  bc[1], sheet[1]), (sheet_end,  sheet[1], polyproline[1]), (polyproline_end, polyproline[1], bc[1]), (1, bc[1],bc[1])),
          'blue': ((0.00,  bc[2], bc[2]), (helix_start,  bc[2], helixR[2]), (helix_end,  helixR[2], bc[2]), (sheet_start,  bc[2], sheet[2]), (sheet_end,  sheet[2], polyproline[2]), (polyproline_end, polyproline[2], bc[2]), (1, bc[2],bc[2]))
        }
#cmap = LinearSegmentedColormap('SecondaryStructureHard', cdict)
#register_cmap(cmap=cmap)
cmaps['SecondaryStructureHard'] = LinearSegmentedColormap('SecondaryStructureHard', cdict)
cmap_ranges['SecondaryStructureHard'] = [0,1]
# -------------------------
# NEW COLOR SCHEME: secondary structure (expected range: R=[0,1])
def border_mod(v):
	# Old min/max
	#  0               1
	#  |   v           |
	# to:
	# New min/max
	#             0.9  1
	#              |v  |
	old_min = 0.0; old_max=1.0
	new_min = 0.9; new_max=1.0
	return new_min + (new_max-new_min)*(v-old_min)/(old_max-old_min)
#
#                         white  white                  white                   red (ish)                                   red        red                             red(ish)  white                  white                 blue(ish)                                 blue      blue                               blue(ish)        cyan                                      cyan(ish)  white       white white
cdict = {  #                  |      |                      |                     |                                           |          |                                   |       |                      |                    |                                         |         |                                  |                |                                             |       |           |     |                  
           'red': ((0.00,  bc[0], bc[0]), (helix_start,  bc[0], border_mod(helixR[0])), ((helix_start+helix_end)/2.0,  helixR[0], helixR[0]), (helix_end,  border_mod(helixR[0]), bc[0]), (sheet_start,  bc[0], border_mod(sheet[0])), ((sheet_start+sheet_end)/2.0,  sheet[0], sheet[0]), (sheet_end,  border_mod(sheet[0]), polyproline[0]), (polyproline_end, border_mod(polyproline[0]), bc[0]), (1, bc[0],bc[0])),
         'green': ((0.00,  bc[1], bc[1]), (helix_start,  bc[1], border_mod(helixR[1])), ((helix_start+helix_end)/2.0,  helixR[1], helixR[1]), (helix_end,  border_mod(helixR[1]), bc[1]), (sheet_start,  bc[1], border_mod(sheet[1])), ((sheet_start+sheet_end)/2.0,  sheet[1], sheet[1]), (sheet_end,  border_mod(sheet[1]), polyproline[1]), (polyproline_end, border_mod(polyproline[1]), bc[1]), (1, bc[1],bc[1])), 
          'blue': ((0.00,  bc[2], bc[2]), (helix_start,  bc[2], border_mod(helixR[2])), ((helix_start+helix_end)/2.0,  helixR[2], helixR[2]), (helix_end,  border_mod(helixR[2]), bc[2]), (sheet_start,  bc[2], border_mod(sheet[2])), ((sheet_start+sheet_end)/2.0,  sheet[2], sheet[2]), (sheet_end,  border_mod(sheet[2]), polyproline[2]), (polyproline_end, border_mod(polyproline[2]), bc[2]), (1, bc[2],bc[2]))
        }
#cmap = LinearSegmentedColormap('SecondaryStructure', cdict)
#register_cmap(cmap=cmap)
cmaps['SecondaryStructure'] = LinearSegmentedColormap('SecondaryStructure', cdict)
cmap_ranges['SecondaryStructureHard'] = [0,1]
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
for color in  list(cdict.keys()):
	for i in range(len(cdict[color])):
		if minpos == False:
			minpos = cdict[color][i][0]
		if maxpos == False:
			maxpos = cdict[color][i][0]
		if minpos > cdict[color][i][0]:
			minpos = cdict[color][i][0]
		if maxpos < cdict[color][i][0]:
			maxpos = cdict[color][i][0]
for color in  list(cdict.keys()):
	for i in range(len(cdict[color])):
		cdict[color][i][0] = float(cdict[color][i][0]-minpos)/(maxpos-minpos)
#cmap = LinearSegmentedColormap('SecondaryStructureFourColor', cdict)
#register_cmap(cmap=cmap)
cmaps['SecondaryStructureFourColor'] = LinearSegmentedColormap('SecondaryStructureFourColor', cdict)
cmap_ranges['SecondaryStructureFourColor'] = [-1,1]



def display_cmaps():
    """Display all custom colormaps defined in ``cmaps`` as horizontal colorbars.

    The colormaps are shown in insertion order with their names as titles. The figure
    and axes are returned so callers can further customize or save the visualization.

    Returns
    -------
    tuple: Matplotlib figure and axes from ``plt.subplots`` used to render the colorbars.
    """
    cmapkeys = list(cmaps.keys())
    fig, axs = plt.subplots(len(cmapkeys),figsize=(6,len(cmapkeys)*1))

    fig.suptitle(r'Custom CMAPs available in $\it{backmap.local\_colormaps.cmaps}$')
    for keyindex in range(len(cmapkeys)):
        cmap_key  = cmapkeys[keyindex]
        cmap_axes = axs[keyindex]
        cmap = cmaps[cmap_key]
        #
        cmap_axes.set_title(cmap_key)

        # 2. Create a Normalize object (optional, but good for controlling the color range)
        # This maps data values (e.g., 0 to 1) to the 0-1 range of the colormap
        norm = colors.Normalize(vmin=-1, vmax=1)

        # 3. Create a ScalarMappable object
        # This object connects the colormap and the normalization to a "mappable" entity
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([]) # An empty array is sufficient as we are not mapping actual data

        # Draw the colorbar
        cbar = plt.colorbar(sm, cax=cmap_axes, orientation="horizontal")

        # Optional: Set a title for the colorbar
        #cbar.set_label(None)
    plt.tight_layout()
    plt.show()
    return fig, axs


