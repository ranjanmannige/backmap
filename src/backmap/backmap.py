"""Main module. Example usage: ``backmap.R()``, which is equal to ``backmap.backmap.R()``"""
import argparse
import sys
import logging
import pandas as pd
from typing import Union
from . import utils
from . import local_colormaps
import matplotlib.patches as patches # for drawing markers

import warnings
#warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning)
# Using the following warning because plt.gca() give out a DeprecationWarning when pytest runs in Python 3.9
warnings.filterwarnings("ignore", category=DeprecationWarning) 


# Standard imports
import sys,re,os,math

# Commonly available imports
import copy,glob
import numpy as np



# matplotlib imports
import matplotlib.pyplot as plt




_logger = logging.getLogger(__name__)

helpme = r"""
============================================
  ____             _    __  __          _____  
 |  _ \           | |  |  \/  |   /\   |  __ \ 
 | |_) | __ _  ___| | _| \  / |  /  \  | |__) |
 |  _ < / _` |/ __| |/ / |\/| | / /\ \ |  ___/ 
 | |_) | (_| | (__|   <| |  | |/ ____ \| |     
 |____/ \__,_|\___|_|\_\_|  |_/_/    \_\_|     
                       (Multi-angle Picture)                                             

This tool provides easily readable "pictures" of protein conformations, 
ensembles, and trajectories saved as either a combined protein databank 
(PDB) structure file, or a directory of such files, and produces graphs.
-----
Usage
-----
The Ramachandran number concept is discussed in the following manuscripts (this tool is discussed in the first reference):
1. Mannige (2017) "A simpler Ramachandran number can simplify the life of a protein simulator" PeerJ 6:e5745. Full Text: https://doi.org/10.7717/peerj.5745
2. Mannige, Kundu, Whitelam (2016) "The Ramachandran number: an order parameter for protein geometry" PLoS ONE 11(8):e0160023. Full Text: https://doi.org/10.1371/journal.pone.0160023
============================================
"""


#sys.path.insert(0, "./local_imports/") # for the local imports
#import Geometry, PeptideBuilder, locallib


# Default values
signed = 0
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

# python plotmap.py -pdb tests/pdbs/nanosheet_birth_U7.pdb
# python plotmap.py -pdb tests/pdbs/1mba.pdb
# python plotmap.py -pdb tests/pdbs/1xqq.pdb
# python plotmap.py -pdb tests/pdbs/2fft.pdb


forcedmax = False
forcedmin = False

show_graphs = 1
default_fontsize = 22
colorbarXscaling = 0.08
defaultXscaling  = 2.0

SCALE = 10.0 # For the postscript output

#Options: 'Chirality','Chirality_r','ChiralityFourColor','ChiralityFourColor_r'
#          'SecondaryStructure','SecondaryStructureFourColor'
rcode_cmap = local_colormaps.cmaps['ChiralityFourColor']


def normalized_ramachandran_number(phi:float,psi:float,signed=False):
    """Return the normalized Ramachandran number for given backbone angles (``phi`` and ``psi``).

    Args:
        phi (float): Backbone phi angle in degrees of range [-180,180].
        psi (float): Backbone psi angle in degrees of range [-180,180].
        signed (bool): If ``True``, return a signed value where ``psi < phi`` yields a negative number.

    Returns:
        float: Normalized Ramachandran number randing ``[0, 1]``, or ``[-1, 1]`` when ``signed`` is ``True``.
    """
    # Calculating the normalized R number
    r = (phi+psi+360)/720.0
    #
    # Signing it if requested
    if signed:
        if psi < phi:
            r = r * -1.0
        #
    #
    return r
#

def ramachandran_number(phi,psi,signed=False):
    """Different abbreviated ways to request the ramachandran number
    Refer to :py:func:`backmap.backmap.normalized_ramachandran_number` for usage details.
    """
    return normalized_ramachandran_number(phi,psi,signed)
#

def r(phi,psi,signed=False):
    '''
    Different abbreviated ways to request the ramachandran number
    Refer to :py:func:`backmap.backmap.normalized_ramachandran_number` for usage details.
    '''
    return normalized_ramachandran_number(phi,psi,signed)
#

def R(phi,psi,signed=False):
    '''
    Different abbreviated ways to request the ramachandran number
    Refer to :py:func:`backmap.backmap.normalized_ramachandran_number` for usage details.
    '''
    return normalized_ramachandran_number(phi,psi,signed)
#

def write_image(fn_base, figure_object=None, write=True, show=True):
    """Save a figure to disk and optionally display it.

    Args:
        fn_base (str): Base filename (without extension) for the output files.
        figure_object (matplotlib.figure.Figure | None): Figure to save; defaults to the current figure when ``None``.
        write (bool): If ``True``, write ``.eps`` and ``.png`` files to disk.
        show (bool): If ``True``, display the figure after saving.

    Returns:
        None
    """
    if figure_object is None:
        figure_object = plt.gcf()
    if write:
        print("\tSaving to:"+fn_base+'.[eps|png]')
        figure_object.savefig(fn_base+'.eps',dpi=200,bbox_inches='tight')
        figure_object.savefig(fn_base+'.png',dpi=200,bbox_inches='tight')
    if show: 
        plt.show()
    return figure_object
#

def set_nan_color(cmap,color='green'):
    cmap.set_bad(color, 1.)
    return cmap

def draw_xyz(X,Y,Z, ylim=False, cmap='Greys', missing_color="green", xlabel=False,ylabel=False,zlabel=False,title=False,vmin=None,vmax=None):
    """Render a heatmap-style plot for spatially indexed values.

    Args:
        X (Sequence[float]): X-coordinates corresponding to each value in ``Z``.
        Y (Sequence[float]): Y-coordinates corresponding to each value in ``Z``.
        Z (Sequence[float]): Values to plot at each (X, Y) location.
        ylim (tuple[float, float] | bool, optional): Y-axis limits; falsy to leave unchanged.
        cmap (str or matplotlib.colors.Colormap, optional): Colormap name or object to use.
        missing_color (str): A color that will be used when a residue is missing (default: 'green')
        xlabel (str or bool, optional): X-axis label; falsy to skip.
        ylabel (str or bool, optional): Y-axis label; falsy to skip.
        zlabel (str or bool, optional): Colorbar label; falsy to skip.
        title (str or bool, optional): Plot title; falsy to skip.
        vmin (float, optional): Lower bound for colormap normalization.
        vmax (float, optional): Upper bound for colormap normalization.

    Returns:
        tuple[bool, matplotlib.axes.Axes]: ``(True, ax)`` where ``ax`` is the axes containing the plot.
    """
    #
    if cmap in local_colormaps.cmaps:
        cmap = local_colormaps.cmaps[cmap]
    else:
        cmap = plt.get_cmap(cmap)
    #
    cmap = set_nan_color(cmap,color=missing_color)
    
    aspect = 2.
    if len(set(X)) == 1:
        # Some structures have only one model, which is too thin, 
        # so we add another row with everything else being identical 
        # except that all new xs are skewed by 1
        X = list(X) + list(np.array(X)+1)
        Y = list(Y) + list(Y)
        Z = list(Z) + list(Z)
        # Also, the aspect ratio needs to be reset, since we are dealing with only one column
        aspect = .2
    #
    # Getting unique values for X
    setX = list(sorted(set(X)))
    # Getting unique values for Y 
    setY = list(sorted(set(Y)))
    
    # Code that offsets the values
    if 1: 
        # We want whole numbers to be situated at the middle of each column, not at the beginning and end
        # X
        # Getting the grid step size
        xsteps = []
        for i in range(1,len(setX)):
            xsteps.append(setX[i]-setX[i-1])
        xstep = np.median(xsteps)
        # Making the offset
        X = np.array(X)-xstep
        # Y
        # Getting the grid step size
        ysteps = []
        for i in range(1,len(setY)):
            ysteps.append(setY[i]-setY[i-1])
        ystep = np.median(ysteps)
        # Making the offset
        Y = np.array(Y)-ystep
        
        # Resetting the sorted unique values
        # Getting unique values for X
        setX = sorted(set(X))
        # Getting unique values for Y 
        setY = sorted(set(Y))
            
    
    #
    # A dictionary containing X values and their indices once ordered 
    X_to_ix = dict([[setX[ix],ix] for ix in range(len(setX))])
    # Creating a new array of indices instead of values
    Xix     =      [   X_to_ix[v] for  v in              X  ]
    
    # A dictionary containing X values and their indices once ordered 
    Y_to_ix = dict([[setY[ix],ix] for ix in range(len(setY))])
    # Creating a new array of indices instead of values
    Yix     =      [   Y_to_ix[v] for  v in              Y  ]
    
    # Creating an empty array with the right dimensions
    z_array = np.zeros((len(setY),len(setX))) * np.nan
    # Setting values of Z based on their position in the matrix
    z_array[Yix, Xix] = Z

    # Initiating the figure
    ax = plt.gca()
    #fig, ax = plt.subplots()

    # Drawing the main part of the figure
    im = plt.imshow(z_array,origin='lower',cmap=cmap,vmin=vmin,vmax=vmax,interpolation='nearest', extent=[min(X),max(X),min(Y),max(Y)])
    
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    cb = plt.colorbar(im,fraction=0.023, pad=0.04)
    [i.set_linewidth(1.5) for i in ax.spines.values()]
    
    # Setting labels
    if xlabel: plt.xlabel(xlabel, fontsize=15);
    if ylabel: plt.ylabel(ylabel, fontsize=15);
    if zlabel: cb.ax.set_title(zlabel,  rotation=0,fontsize=15);
    
    # Setting title
    if title: plt.title(title,fontsize=16)
    
    # limiting y if specified
    if ylim: plt.ylim(ylim);
    
    # Setting the aspect ratio
    utils.forceAspect(aspect,ax=ax)
    
    # Neatening things out
    plt.tight_layout()
    #
    # To see this graph and quit, uncomment the following:
    #plt.show(); exit();
    return True, ax
#

def obtain_ramachandran_dataframe(coordinates_df, signed=False):
    """Aggregate residue-level dihedral angles from per-atom backbone coordinates.

    Args:
        coordinates_df (pd.DataFrame): Backbone atom coordinates with at least
            ``model``, ``chain``, ``resid``, ``atom``, ``X``, ``Y``, and ``Z`` columns.
        signed (bool): Whether to compute the signed Ramachandran number ``R``.

    Returns: 
        pd.DataFrame: Per-residue records containing model, chain, residue
        identity, phi, psi, and Ramachandran number values. Missing neighbors
        leave phi/psi/R as ``False``.
    """
    keys_to_keep = ['model','chain','resname','resid']
    coordinates_df['unique_model_and_chain'] = coordinates_df['model'].astype(str)+'_'+coordinates_df['chain']
    final_records = []
    for model_and_chain in coordinates_df['unique_model_and_chain'].unique():
        cdf = coordinates_df[coordinates_df['unique_model_and_chain']==model_and_chain]
        
        # converting to dict of resid to dict
        resid_to_records = {}
        for ix, row in cdf.iterrows():
            resid = row['resid']
            if not resid in resid_to_records:
                resid_to_records[resid] = []
            #
            resid_to_records[resid].append(dict(row))

        for resid in sorted(resid_to_records.keys()):
            i = resid
            im = i-1
            ip = i+1
            
            neighbors_found = 0
            
            # 
            # Internal variables: a------b------c------d------e
            # Atom names          C-     N      CA     C      N+
            # Dihedral notions:   |<-------phi-------->|
            #                            |<------psi--------->|
            #
            b = utils.get_coord_by_key_value_match(key='atom',value='N',list_of_dict=resid_to_records[i])
            c = utils.get_coord_by_key_value_match(key='atom',value='CA',list_of_dict=resid_to_records[i])
            d = utils.get_coord_by_key_value_match(key='atom',value='C',list_of_dict=resid_to_records[i])
            
            phi,psi,rho = np.nan, np.nan, np.nan
            #print(len(b)+len(c)+len(d))
            if len(b)+len(c)+len(d) == 9:
                if im in resid_to_records:
                    a = utils.get_coord_by_key_value_match(key='atom',value='C',list_of_dict=resid_to_records[im])
                    if len(a) == 3:
                        phi = utils.calculate_dihedral_angle(np.array([a,b,c,d]))
                if ip in resid_to_records:
                    e = utils.get_coord_by_key_value_match(key='atom',value='N',list_of_dict=resid_to_records[ip])
                    if len(e) == 3:
                        psi = utils.calculate_dihedral_angle(np.array([b,c,d,e]))

                # a = model_to_chain_to_resid_atom_to_vals[model][chain][im]["c"] # resid_to_coordC[before]
                # e = model_to_chain_to_resid_atom_to_vals[model][chain][ip]["n"]  # resid_to_coordN[after]
                
                if not ( np.isnan(phi) and np.isnan(psi)):
                    rho= normalized_ramachandran_number(phi,psi,signed)
            
            ignore_cols = {'atom','X','Y','Z','unique_model_and_chain'}
            if len(resid_to_records[resid]):
                current_row = {key:val for key,val in resid_to_records[resid][0].items()
                               if key not in ignore_cols}
                current_row.update({"phi":phi, "psi": psi, "R":rho} )
            final_records.append(current_row)
    
    
    '''

    Inputs:
        pdbfn:str                    Filepath that points to either a PDB *OR* a directory that contains .pdb files
        signed:bool=False            Switch to using the signed Ramachandran number (default=False)
    Outputs:
        structure_df      A pandas dataframe describing each backbone atom, containing the following keys:
                          'model': PDB chain MODEL number ()
                          'chain': PDB CHAIN
                          'resname': three letter residue name
                          'phi': the backbone dihedral angle phi 
                          'psi': the backbone dihedral angle psi
                          'R': the Ramachandran number
    '''
    

    return pd.DataFrame(final_records)

def read_pdb_from_pdbfn_or_filehandle(pdbfn_or_filehandle, structure_df=pd.DataFrame()):
    """Parse PDB content, compute Ramachandran metrics, and merge with an existing dataframe.

    Args:
        pdbfn_or_filehandle: Path to a PDB file or an open file-like object accepted by ``utils.read_pdb``.
        structure_df (pd.DataFrame): Existing dataframe of processed structures to append to; empty by default.

    Returns:
        pd.DataFrame: Combined residue-level dataframe with continuous model numbering and Ramachandran values.
    """
    
    #
    # READ PDB in the form of a matrix with columns ['model','chain','resid','R']
    raw_pdb_data = utils.read_pdb(filename_or_filehandle=pdbfn_or_filehandle)
    raw_pdb_data = utils.bytecheck(raw_pdb_data)
    #
    latest_structure_df = obtain_ramachandran_dataframe(raw_pdb_data,signed)
    #latest_structure_df, latest_structure = calculate_R_from_raw_pdb_data(raw_pdb_data, signed=signed)
    #
    if len(latest_structure_df)==0:
        print("WARNING: PDB FILE '{}' NOT FOUND")
        return pd.DataFrame()
    #
    # Sorting by model number
    latest_structure_df = latest_structure_df.sort_values(by=['chain','model','resid'],ascending=True)
    #
    # Assigning continuous model numbers 
    old_model_numbers = list(sorted(set(latest_structure_df['model'])))
    new_model_numbers = range(1, len(old_model_numbers)+1)
    model_number_conversion = dict(zip(old_model_numbers,new_model_numbers))
    latest_structure_df['model'] = [model_number_conversion[mn] for mn in latest_structure_df['model']]
    #
    # If this is not the first structure file that is being processed...
    if len(structure_df) > 0:
        structure_df_max_model = max(set(latest_structure_df['model']))
        latest_structure_df['model'] = latest_structure_df['model']+structure_df_max_model
    #
    structure_df = pd.concat([structure_df,latest_structure_df.copy()])
    return structure_df

def process_PDB(pdbfn:str, signed:bool=False):
    """
    Load PDB file(s) and compute residue-level Ramachandran metrics.

    Args:
        pdbfn (str): Path to a PDB file or to a directory containing ``.pdb[.gz|zip|tar]`` files.
        signed (bool): If ``True``, compute signed Ramachandran numbers instead of the default unsigned values.

    Returns:
        pd.DataFrame: Per-residue records with model, chain, residue identity, phi,
        psi, and Ramachandran number columns. Returns an empty dataframe when no
        input PDBs are found.
    """

    # Since the user can either provide one PDB file *or* one 
    # PDB file-containing directory, we need to resolve these PDBs
    list_pdbfilenames,pdbdir = utils.get_pdb_filenames(pdbfn)
    #
    
    # for list_pdb_filenames:
    #     utils.get_filename_and_filehandle(
    if len(list_pdbfilenames) == 0:
        return pd.DataFrame()

    #NAME = os.path.basename(list_pdbfilenames[0])[:-len(".pdb")]
    #target_base = target_dir.rstrip("/")
    #
    structure = np.array([]) # depricated
    structure_df = pd.DataFrame()
    for _pdbfn in list_pdbfilenames:#[:10]:
        # Each PDB filename itself can be an archive of PDB files, 
        #
        # Opening the file, assuming that it could be a list of files (can be, if 
        # the file ends with .tgz .tar.gz or .zip)
        filehandles = utils.return_filehandle_from_gz_zip_or_normal_file(_pdbfn)
        #
        for f in filehandles:
            # Check if the PDB has no subunit IDs, and then check if segnames exist (right most column)
            # and renaming the subunit IDs alphabetically and then numerically
            #
            final_pdbfn = f.name
            final_pdbfh = f
            #
            structure_df = read_pdb_from_pdbfn_or_filehandle(final_pdbfn, structure_df=structure_df)
            
        #
    structure_df.attrs['pdbfn'] = pdbfn
    structure_df.attrs['pdb_name'] = os.path.split(pdbfn)[-1]
    structure_df.attrs['signed'] = signed
    #
    # models, residues, Rs = list(structure_df['model']),list(structure_df['resid']),list(structure_df['R'])
    # print(f'{len(set(models))=} {models=},')
    # print(f'{len(set(residues))=} {residues=}')
    # print(f'{len(set(Rs))=} {Rs=}')
    return structure_df


def fill_in_missing_resids(structure_df, fill_with=False):
    """Insert placeholder rows for missing residue IDs per chain and model.

    Args:
        structure_df (pd.DataFrame): Backbone geometry table containing
            ``chain``, ``model``, and ``resid`` columns.
        fill_with (Any): Value used to populate missing residue rows; defaults
            to ``False``.

    Returns:
        pd.DataFrame: Reindexed dataframe spanning continuous residue ranges for
        each chain and model, filled with ``fill_with`` where data was absent.
    """

    final_structure_df = pd.DataFrame()
    for chain in sorted(structure_df['chain'].unique()):
        _cdf = structure_df[structure_df['chain']==chain]
        # 
        all_residue_ids = set(_cdf['resid'])
        #
        resid_range     = list(range(min(all_residue_ids),max(all_residue_ids)+1))
        for model in sorted(_cdf['model'].unique()):
            _model_df = _cdf[_cdf['model']==model]
            _model_df = _model_df.set_index('resid')
            _model_df = _model_df.reindex(resid_range, fill_value=fill_with)
            _model_df['model'] = model 
            _model_df['chain'] = chain
            _model_df = _model_df.reset_index()
            final_structure_df = pd.concat([final_structure_df,_model_df])
    return final_structure_df

def mark_figure(df,ax=None,mark_column='mark'):
    """Add right-side markers to rows flagged in ``mark_column`` on the current plot.

    Args:
        df (pd.DataFrame): Table containing a ``resid`` column and a boolean
            column used to flag residues for annotation.
        ax (matplotlib.axes.Axes, optional): Axes to annotate; defaults to the
            current matplotlib axes.
        mark_column (str): Column name containing boolean markers; defaults to
            ``'mark'``.

    Returns:
        matplotlib.axes.Axes: Axes with rectangle markers added for flagged rows.
    """
    
    if ax is None:
        ax = plt.gca()
    if mark_column in df.columns:
        for ix,row in df[df[mark_column]==True].iterrows():
            resid = row['resid']
            xmin,xmax = plt.xlim()
            width = (xmax-xmin)*0.1
            rect = patches.Rectangle((xmax, resid), width, 1, #transform=ax.transAxes,
                                linewidth=0, edgecolor='k', facecolor='k',clip_on=False)
            # Add the patch to the Axes
            ax.add_patch(rect)
    return ax
        


def draw_figures(structure_df, output_dir='', write=True, show=True):
    """Generate per-chain Ramachandran visualizations and optionally save/show them.

    Args:
        structure_df (pd.DataFrame): Per-residue Ramachandran data containing at
            least ``model``, ``chain``, ``resid``, and ``R`` columns.
        output_dir (str, optional): Directory to write figures; defaults to
            ``<pdbdir>/reports`` when empty and ``write`` is ``True``.
        write (bool): If ``True``, save ``.eps`` and ``.png`` files for each plot.
        show (bool): If ``True``, display figures after creation.

    Returns:
        tuple[bool, dict]: ``(True, figures)`` where ``figures`` maps plot labels
        to matplotlib figure objects for the value maps, histograms, RMSF, and
        RMSD heatmaps computed per chain.
    """
    # All figure elements (of plt.gcf() type) are returned with the figures dict 
    figures = {}
    #
    structure_df = fill_in_missing_resids(structure_df, fill_with=np.nan)
    #
    unique_chains = list(sorted(set(structure_df['chain'])))
    #
    pdbfn = structure_df.attrs['pdbfn']
    pdbdir = os.path.dirname(pdbfn)
    if write:
        if len(output_dir) == 0:
            output_dir = pdbdir+"/reports/"

        # Making the report directory, if not created already
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
    #
    name = os.path.split(pdbfn)[-1][:-len('.pdb')]
    print(" ---- \t---------")
    #
    vmin           =  0 
    vmax           =  1
    ss_cmap        = 'SecondaryStructure'
    chirality_cmap = 'Chirality'
    if signed:
        vmin           = -1 
        vmax           =  1
        ss_cmap        = ss_cmap        + 'FourColor'
        chirality_cmap = chirality_cmap + 'FourColor'
    #
    structure_df = fill_in_missing_resids(structure_df, fill_with=np.nan)
    #
    print(" 1  \tRamachandran number (PDB: %s)"%(pdbfn))
    different_plots = []
    # setting the name of the colormap
    for cmap in ['Greys',ss_cmap,chirality_cmap]: #, 'Chirality_r', 'SecondaryStructureHard']:
        # DRAWING A SINGLE GRAPH
        for chain in unique_chains:
            final_name = name
            if len(chain.rstrip()):
                final_name+='-'+str(chain)
            #
            # Getting the X,Y,Z values for each entry
            #models, residues, Rs = grouped_data[chain]
            _cdf = structure_df[structure_df['chain']==chain]
            models, residues, Rs = list(_cdf['model']),list(_cdf['resid']),list(_cdf['R'])
            #
            # Finally, creating (but not showing) the graph 
            if True:
                assertion, ax = draw_xyz(X = models  ,      Y = residues  ,     Z = Rs
                        , xlabel =r'Frame #', ylabel =r"Residue #",zlabel =r'$\mathcal{R}$'
                        , title=r'Per-residue $\mathcal{R}$; CMAP: '+cmap+'\nPDB: ' + final_name
                        ,  cmap = cmap    ,  vmin=vmin, vmax=vmax)
                mark_figure(_cdf, ax)
                #
                #plt.show()
                # Now, we display the graph:
                FN = None
                if write:
                    FN = output_dir+'/pdb_%s_r_%s' %(final_name,cmap)
                figure_object = write_image(FN, write=write, show=show)#,new_fig)
                figures[f'chain_{chain}_val_cmap{cmap}'] = figure_object
                #
    #
    # Getting only those values for the particular chain 
    print(" 2.  \tHistogram (PDB: 1xqq)")
    for chain in unique_chains:
        final_name = name
        if len(chain.rstrip()):
            final_name+='-'+str(chain)
        
        # Getting the X,Y,Z values for each entry
        #models, residues, Rs = grouped_data[chain]
        _cdf = structure_df[structure_df['chain']==chain]
        models, residues, Rs = list(_cdf['model']),list(_cdf['resid']),list(_cdf['R'])
        
        X = []; Y=[]; Z=[]; # Will set X=model, Y=R, Z=P(R)
        # Bundling the three lists into one 2d array
        new_data =  np.array(list(zip(models,residues,Rs)))
        # Getting all R values, model by model
        for m in sorted(set(new_data[:,0])): # column 0 is the model column
            # Getting all Rs for that model #
            current_rs = new_data[np.where(new_data[:,0]==m)][:,2] # column 2 contains R
            # Getting the histogram
            a,b = np.histogram(current_rs,bins=np.arange(vmin,vmax+0.0001,0.01))
            max_count = float(np.max(a))
            for i in range(len(a)):
                X.append(m)
                Y.append((b[i]+b[i+1])/2.0)
                Z.append(a[i]/float(np.max(a)))
        
        if True:
            # Finally, creating (but not showing) the graph 
            plt.clf()
            is_true, ax = draw_xyz(X = X       ,      Y = Y  ,                Z = Z
            ,xlabel ='Frame #', ylabel =r"$\mathcal{R}$",zlabel =r"$P'(\mathcal{R})$:"
                ,cmap = 'Greys', ylim=[vmin,vmax],title=r'Per-model $\mathcal{R}$-histogram'+'\nPDB: '+final_name)
            plt.yticks(np.arange(vmin,vmax+0.00001,0.2))
            # Now, we display the graph:
            FN = None
            if write:
                FN = output_dir+'/pdb_%s_his'%(final_name)
            figure_object = write_image(FN, write=write, show=show)
            figures[f'chain_{chain}_his_cmap{cmap}'] = figure_object
            #print("\tSaved to:",FN)
    #
    #
    print(" 3.  \tRMSF (PDB: {})".format(name))
    for chain in unique_chains:
        final_name = name
        if len(chain.rstrip()):
            final_name+='-'+str(chain)
        
        # Getting the X,Y,Z values for each entry
        #models, residues, Rs = grouped_data[chain]
        _cdf = structure_df[structure_df['chain']==chain]
        models, residues, Rs = list(_cdf['model']),list(_cdf['resid']),list(_cdf['R'])
        
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
            if True:
                plt.clf()
                is_true, ax = draw_xyz(X = X       ,      Y = Y  ,                Z = Z
                    ,xlabel ='Frame #', ylabel =r"Residue #",zlabel ="$D_{-1}$"
                    ,cmap = 'Blues', title='Per-residue deviation $D_{-1} = |\\mathcal{R}_t - \\mathcal{R}_{t-1}|$\nPDB: '+ final_name)
                
                # Now, we display the graph:
                FN = None
                if write:
                    FN = output_dir+'/pdb_%s_rmsf'%(final_name)
                figure_object = write_image(FN, write=write, show=show)
                figures[f'chain_{chain}_rmsf_cmap{cmap}'] = figure_object
        else:
            print('\tChain "%s" has only one model. Not drawing this graph.' %(chain))
    #
    #
    print(f' 4.  \tRMSD (PDB: {name})')
    for chain in unique_chains:
        final_name = name
        if len(chain.rstrip()):
            final_name+='-'+str(chain)
        
        # Getting the X,Y,Z values for each entry
        #models, residues, Rs = grouped_data[chain]
        _cdf = structure_df[structure_df['chain']==chain]
        models, residues, Rs = list(_cdf['model']),list(_cdf['resid']),list(_cdf['R'])
        
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
            
            
            if True:
                # Finally, creating (but not showing) the graph 
                plt.clf()
                is_true, ax = draw_xyz(X = X, Y = Y, Z = Z,
                                xlabel = 'Frame #', ylabel = "Residue #", zlabel ="$D_{1}$",
                                cmap = 'Reds', title= r'Per-residue deviation $D_{1} = |\mathcal{R}_t - \mathcal{R}_{1}|$'+'\nPDB: '+final_name)
                #plt.yticks(np.arange(0,1.00001,0.2))
                # Now, we display the graph:
                FN = None
                if write:
                    FN = output_dir+'/pdb_%s_rmsd'%(final_name)
                figure_object = write_image(FN, write=write, show=show)
                figures[f'chain_{chain}_rmsd_cmap{cmap}'] = figure_object
        else:
            print('\tChain "%s" has only one model. Not drawing this graph.' %(chain))
    #
    plt.clf()
    #
    return True, figures
#

# from ._version import get_versions
# __version__ = get_versions()['version']
# del get_versions

if __name__ == "__main__":
    print('Please use "python -m backmap" for the standalone version of backmap.')
    
#from ._version import get_versions
#__version__ = get_versions()['version']
#del get_versions
