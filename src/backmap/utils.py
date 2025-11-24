import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import string,re
try:
    ascii_uppercase = string.ascii_uppercase
except:
    # A catch for Python 2.7 (will be deprecated in the next version)
    ascii_uppercase = string.uppercase

# Using the following flag to either use the Biopython PDB parser if True, vs
# the in house PDB parser if False
biopython_is_installed = False
try:
    # Possibly use: if 'Bio' in sys.modules:
    from Bio import PDB
    biopython_is_installed = True
except:
    biopython_is_installed = False

def median_filter(vals,nearest_neighbors=1):
    '''
    Simple smoothing function that returns the median value of 
    each element and its nearest neighbors (+/-).
    Inputs:
        vals [List of floats]
        nearest_neighbors [int]    The number of neighbors to use before and 
                                after the existing element.
                                Edge elements are not treated.
    Outouts:
        new_vals [List of floats]    Modified list with the median filter applied.
    '''
    new_vals = []
    len_vals = len(vals)
    for i in range(len_vals):
        val = vals[i]
        if i-nearest_neighbors >= 0 and i+nearest_neighbors < len_vals:
            val = np.median(vals[i-nearest_neighbors:i+nearest_neighbors+1])
        new_vals.append(val)
    return new_vals
#

# From: https://stackoverflow.com/questions/7965743/how-can-i-set-the-aspect-ratio-in-matplotlib
def forceAspect(aspect,ax=False):
    '''
    Setting the specific plot axis ratio. Works on the current plot/axis.
    Input:
        aspect [float] Desired aspect ratio  
                (E.g., if aspect == 2, then the x axis will have twicce the length of the y)
        ax [matplotlib ax object]    The existing plot to modify
    Output:
        None    The 
    '''
    if not ax: ax=plt.gca()
    extent = plt.axis()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
#
# Three-to-one amino acid conversion lookup
aa_three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','XXX':'X'}
#
# ===================================================================================
def calculate_dihedral_angle(p):
    '''
    Inputs:
        p [4-list of coordinates]    A set of four three dimensional coordinates 
                                    [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3],[x4,y4,z4]]
    Returns:
        d     Dihedral angle (range [-180,180])
    '''
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
    # Normalize vectors
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    b1 = b[1] / np.linalg.norm(b[1])
    x = np.dot(v[0], v[1])
    m = np.cross(v[0], b1)
    y = np.dot(m, v[1])
    d = np.degrees(np.arctan2( y, x ))
    return d 
    
def is_filehandle(x):
    return isinstance(x, io.IOBase)


from zipfile import ZipFile
from gzip import open as gzopen
import glob

def return_zip_file_handle(pdbfn, mode='r'):
    with ZipFile(pdbfn, 'r') as zip_archive:
        all_files_in_archive  = zip_archive.namelist()
        first_file_in_archive = all_files_in_archive[0]
        return zip_archive.open(first_file_in_archive, mode=mode)
#
def return_gz_filehandle(pdb_fn,mode='rt'):
    """Open a gzipped PDB file and return a file-like handle.
    
    Args:
        pdb_fn (str | os.PathLike): Path to the `.pdb.gz` file to open.
        mode (str, optional): Mode passed to ``gzip.open`` (e.g., ``'rt'`` for
            text or ``'rb'`` for binary). Defaults to ``'rt'``.

    Returns:
        io.BufferedReader | io.TextIOWrapper: Handle to the gzipped file opened
        with the requested mode.
    """
    return gzopen(pdb_fn,mode)
#
def return_filehandle_from_gz_zip_or_normal_file(pdbfn):
    """Return a readable handle for `.pdb`, `.pdb.gz`, or `.pdb.zip` inputs.

    Args:
        pdbfn (str): Path to a PDB file, optionally gzipped or zipped. For zip
            archives, the first file in the archive is opened.

    Returns:
        IOBase: File-like object opened with the appropriate handler for the
            given extension, or ``None`` if the extension is unsupported.
    """
    dict_extention_to_open_object = {
        '.pdb.zip':return_zip_file_handle,
        '.pdb.gz':return_gz_filehandle,
        '.pdb':open,
    }
    for file_extension in dict_extention_to_open_object.keys():
        if pdbfn[-len(file_extension):].lower() == file_extension:
            return dict_extention_to_open_object[file_extension](pdbfn)

def get_filename_and_filehandle(filename_or_filehandle):
    """Return a `(filename, filehandle)` tuple from a path or existing file object.

    Accepts either a file-like object (uses its `name` attribute) or a filepath
    string pointing to a PDB file, automatically opening `.pdb`, `.pdb.gz`, or
    `.pdb.zip` inputs with the appropriate handler. Raises `TypeError` for other
    input types.
    """
    fn = False
    f = False
    if is_filehandle(filename_or_filehandle):
        fn = f.name
        f = filename_or_filehandle
    elif isinstance(filename_or_filehandle, str) and not isinstance(filename_or_filehandle, list):
        fn = filename_or_filehandle
        f = return_filehandle_from_gz_zip_or_normal_file(fn)
    else:
        raise TypeError("fn_or_filehandle must be a filepath string or file-like object")
    return fn, f

# A PDB reader that uses Biopython
# in case the BipPython module was not installed, an in-house reader exists below
def read_pdb_biopython(fn_or_filehandle):
    """Parse a PDB file with BioPython and extract backbone atoms.

    Args:
        fn_or_filehandle: Path to a PDB file or readable file-like object
            supported by ``get_filename_and_filehandle`` (gz/zip/plain).

    Returns:
        pandas.DataFrame: Rows for backbone atoms (N, CA, C) across models,
        chains, and residues with columns
        ``['model','chain','resname','resid','atom','X','Y','Z']`` where
        coordinates are float32.
    """
    
    fn, f = get_filename_and_filehandle(filename_or_filehandle=fn_or_filehandle)

    # Parsing with PDB... note that "get_structure takes either a filename OR a filehandle"
    p=PDB.PDBParser() #(PERMISSIVE=1)
    structure=p.get_structure(fn[:-len(".pdb")], fn)
    #for model in structure:
    #    print [model.id]

    #model_to_chain_to_resno_atom_to_vals = {}
    # structure (models) -> model -> chain -> residue -> atom
    mol_rows = []
    for model in structure:
        model_number = model.id
        #
        # if not model_number in model_to_chain_to_resno_atom_to_vals:
        #     model_to_chain_to_resno_atom_to_vals[model_number] = {}
        #
        for chain in model:
            chain_id = chain.id
            # if not segname in model_to_chain_to_resno_atom_to_vals[model_number]:
            #     model_to_chain_to_resno_atom_to_vals[model_number][segname] = {}
            
            for residue in chain:
                #print(residue.__dict__)
                #continue
                resname = residue.resname
                resno   = residue.id[1]
                
                # Checking if this is an actual residue 
                # (e.g., ACE has a "C"-named atom, but no 'N' and 'CA')
                if(     'N' in residue.child_dict 
                    and 'CA' in residue.child_dict 
                    and 'C' in residue.child_dict):
                    pass
                else:
                    continue
                #
                for atom in residue:
                    if atom.name in ['N','CA','C']:
                        xyz = tuple(atom.coord)
                        mol_rows.append({'model':model_number,
                                        'chain':chain_id,
                                        'resname':resname,     
                                        'resid':resno, 
                                        'atom':atom.name, 
                                        'X':xyz[0],
                                        'Y':xyz[1],
                                        'Z':xyz[2]})
    df = pd.DataFrame(mol_rows)
    for c in ['X','Y','Z']:
        df[c] = df[c].astype('float32')
    return df
    #
    #return model_to_chain_to_resno_atom_to_vals

# OLD VERSION (IN HOUSE). IT IS FASTER THAN THE CURRENT "read_pdb", WHICH IS BIOPDB RUN, BUT IT IS NOT 
# AS WELL TESTED.
def read_pdb_inhouse(fn_or_filehandle):
    """Parse a PDB file without BioPython and extract backbone atoms.

    Args:
        fn: Path to a PDB file or a readable file-like object understood by
            ``get_filename_and_filehandle``.

    Returns:
        pandas.DataFrame: Backbone atoms (N, CA, C) for each model, chain, and
        residue with columns ``['model','chain','resname','resid','atom','X','Y','Z']``
        where coordinates are float32.
    """
    
    fn, f = get_filename_and_filehandle(filename_or_filehandle=fn_or_filehandle)
    # f = open(fn,"r")
    pdbblock = f.read()
    f.close()
    
    # Checking if the pdbblock is byte encoded, and switching to plain text
    if isinstance(pdbblock, (bytes, bytearray)):
        pdbblock = pdbblock.decode()

    # The regular expression that follow are expected to extract key values
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
    getlines       = re.compile(r"ATOM\s+(?P<atomno>\d+)\s+(?P<atomtype>\S+)\s+(?P<resname>...).(?P<chainname>.)\s+(?P<resno>\d+)\s+(?P<x>\-*\d+\.*\d*)\s+(?P<y>\-*\d+\.*\d*)\s+(?P<z>\-*\d+\.*\d*).{17}(?P<segname>.{5})",re.M)
    getlines_short = re.compile(r"ATOM\s+(?P<atomno>\d+)\s+(?P<atomtype>\S+)\s+(?P<resname>...).(?P<chainname>.)\s+(?P<resno>\d+)\s+(?P<x>\-*\d+\.*\d*)\s+(?P<y>\-*\d+\.*\d*)\s+(?P<z>\-*\d+\.*\d*)",re.M)
    
    resnos = []
    #models = pdbblock.split("\nEND\n")
    models = re.split(r"\nEND|\nMODEL|\nTER",pdbblock)
    
    mol_rows = []
    model_number = -1
    model_to_chain_to_resno_atom_to_vals = {}
    # structure (models) -> model -> chain -> residue -> atom
    
    #t0 = time.time()
    #print "#\treading...",
    for model_index in range(len(models)):
        model = models[model_index]
        if len(model.rstrip()) > 1:
            model_number = model_index
            
            segname_exists = 1
            currentlines = getlines.finditer(model)
            if not getlines.search(model):
                currentlines = getlines_short.finditer(model)
                segname_exists = 0
            
            for i in currentlines:
                vals = i.groupdict()
                atomtype = vals["atomtype"] #line[11:17].lstrip().rstrip()
                
                if atomtype in ["CA", "N", "C"]:
                    resname = vals["resname"]
                    resno = int(vals["resno"]) #int(resno) #int(line[22:26].lstrip().rstrip())
                    xyz = np.array([float(vals["x"]),float(vals["y"]),float(vals["z"])])
                    
                    segname = "A"
                    if segname_exists:
                        segname = vals["chainname"].lstrip().rstrip()
                    
                    mol_rows.append({'model':model_number,
                                     'chain':segname,
                                   'resname':resname,     
                                     'resid':resno, 
                                      'atom':atomtype, 
                                         'X':xyz[0],
                                         'Y':xyz[1],
                                         'Z':xyz[2]})
                    #
                    #model_to_chain_to_resno_atom_to_vals[model_number][segname][resno][atomtype] = xyz
                    #model_to_chain_to_resno_atom_to_vals[model_number][segname][resno]["resname"] = vals["resname"]
                #
            #
            #if not len(model_to_chain_to_resno_atom_to_vals[model_number]):
            #    del model_to_chain_to_resno_atom_to_vals[model_number]
            #    model_number-=1
        #
    #
    df = pd.DataFrame(mol_rows)
    for c in ['X','Y','Z']:
        df[c] = df[c].astype('float32')
    #
    return df
#

# Does a sniff test for whether this is a legitimate pdb file
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
    chainIDpossibilities+=ascii_uppercase # 'A' through 'Z'.
    for i in range(10):
        chainIDpossibilities+=str(i)
    chainIDpossibilities+=ascii_uppercase.lower() # 'a' through 'z'.
    lenchainIDpossibilities = len(chainIDpossibilities)
    largestchainIDindex = 0
    
    made_changes = 0
    f = open(fn,"r")
    lines = f.readlines()
    f.close()
    pdb_is_possibly_problematic = 0
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
#
def read_pdb(fn:str):
    '''
    Reads a PDB and outputs graphs representing Ramachandran number plots
    See: https://doi.org/10.7717/peerj.5745
    Inputs:
        fn:            PDB filename. 
                    Can be multiple structures separated by the MODEL term;
                    see "https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)")
        signed:        Whether the Ramachandran number plots should be signed or not.
    Outputs:
        pdb_df        A pandas dataframe with the following columns: 
                        ['model','chain','resid','R']
        pdb_matrix    Identical to pdb_df, but returned as a numpy matrix of dtype 'O'
    '''
    raw_pdb_data = False
    if biopython_is_installed:
        if check_pdb(fn):
            raw_pdb_data = read_pdb_biopython(fn)
        else:
            raw_pdb_data = read_pdb_inhouse(fn)
    else:
        raw_pdb_data = read_pdb_inhouse(fn)
    #
    return raw_pdb_data
    
#
