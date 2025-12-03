"""Utility functions important for loading and processing PDB files. 
Functions following ``calculate_dihedral_angle()`` are probably not very 
interesting to most."""

import zipfile
import gzip
import io, os, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from zipfile import ZipFile
from gzip import open as gzopen
import tarfile, io # to open tar.gz or tgz files
import string,re
from typing import Union
fn_or_filehandle:Union[str,os.PathLike,io.IOBase]

from typing import List, Tuple, Literal

try:
    ascii_uppercase = string.ascii_uppercase
except:
    # A catch for Python 2.7 (will be deprecated in the next version)
    ascii_uppercase = string.uppercase

# Using the following flag to either use the Biopython PDB parser if True, vs
# the in house PDB parser if False
biopython_is_installed = False
# Moving away from 
# try:
#     # Possibly use: if 'Bio' in sys.modules:
#     from Bio import PDB
#     biopython_is_installed = True
# except:
#     biopython_is_installed = False

# The regular expression that follow are expected to extract key values
# Pattern to extract PDB ATOM fields including atom id, type, residue info, coordinates, and segname
getlines       = re.compile(r"ATOM\s+(?P<atomno>\d+)\s+(?P<atomtype>\S+)\s+(?P<resname>...).(?P<chainname>.)\s+(?P<resno>\d+)\s+(?P<x>\-*\d+\.*\d*)\s+(?P<y>\-*\d+\.*\d*)\s+(?P<z>\-*\d+\.*\d*).{17}(?P<segname>.{5})",re.M)

# Fallback pattern to extract PDB ATOM fields including atom id, type, residue info, coordinates, and segname
getlines_short = re.compile(r"ATOM\s+(?P<atomno>\d+)\s+(?P<atomtype>\S+)\s+(?P<resname>...).(?P<chainname>.)\s+(?P<resno>\d+)\s+(?P<x>\-*\d+\.*\d*)\s+(?P<y>\-*\d+\.*\d*)\s+(?P<z>\-*\d+\.*\d*)",re.M)

#
# Three-to-one amino acid conversion lookup
aa_three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','XXX':'X'}
#
# ===================================================================================
def read_pdb(filename_or_filehandle:Union[str,os.PathLike]):
    """Reads a PDB and outputs graphs representing Ramachandran number plots.
    
    Load a PDB file, preferring BioPython parsing with an in-house fallback.

    The function first tries the BioPython-based parser when the library is
    installed and the file passes ``check_pdb``. Files that look problematic or
    environments without BioPython fall back to the in-house parser to extract
    backbone atoms.
    
    Args:
        filename_or_filehandle (str): Path to a PDB file to read, or file handle. 
        Can be multiple structures separated by the MODEL term.

    Returns:
        pandas.DataFrame: Backbone atoms (N, CA, C) with model, chain, residue
        identifiers and float32 coordinates.
    """
    raw_pdb_data = read_pdb_inhouse(filename_or_filehandle)
    #raw_pdb_data = False
    # if biopython_is_installed:
    #     raw_pdb_data = read_pdb_biopython(fn)
    # else:
    #     raw_pdb_data = read_pdb_inhouse(fn)
    #
    return raw_pdb_data

# --------------------------------------
# MOVING AWAY FROM BIOPYTHON ALLTOGETHER 
# --------------------------------------
# # A PDB reader that uses Biopython
# # in case the BipPython module was not installed, an in-house reader exists below
# def read_pdb_biopython(fn_or_filehandle:Union[str,os.PathLike,io.IOBase]):
#     """Parse a PDB file with BioPython and extract backbone atoms.
#
#     Args:
#         fn_or_filehandle: Path to a PDB file or readable file-like object
#             supported by ``get_filename_and_filehandle`` (gz/zip/plain).
#
#     Returns:
#         pandas.DataFrame: Rows for backbone atoms (N, CA, C) across models,
#         chains, and residues with columns
#         ``['model','chain','resname','resid','atom','X','Y','Z']`` where
#         coordinates are float32.
#     """
#   
#     fn, f = get_filename_and_filehandle(filename_or_filehandle=fn_or_filehandle)
#
#     # Parsing with PDB... note that "get_structure takes either a filename OR a filehandle"
#     p=PDB.PDBParser() #(PERMISSIVE=1)
#     structure=p.get_structure(fn[:-len(".pdb")], fn)
#     #for model in structure:
#     #    print [model.id]
#
#     #model_to_chain_to_resno_atom_to_vals = {}
#     # structure (models) -> model -> chain -> residue -> atom
#     mol_rows = []
#     for model in structure:
#         model_number = model.id
#         #
#         # if not model_number in model_to_chain_to_resno_atom_to_vals:
#         #     model_to_chain_to_resno_atom_to_vals[model_number] = {}
#         #
#         for chain in model:
#             chain_id = chain.id
#             # if not segname in model_to_chain_to_resno_atom_to_vals[model_number]:
#             #     model_to_chain_to_resno_atom_to_vals[model_number][segname] = {}
#            
#             for residue in chain:
#                 #print(residue.__dict__)
#                 #continue
#                 resname = residue.resname
#                 resno   = residue.id[1]
#                 #
#                 # Checking if this is an actual residue 
#                 # (e.g., ACE has a "C"-named atom, but no 'N' and 'CA')
#                 if(     'N' in residue.child_dict 
#                     and 'CA' in residue.child_dict 
#                     and 'C' in residue.child_dict):
#                     pass
#                 else:
#                     continue
#                 #
#                 for atom in residue:
#                     if atom.name in ['N','CA','C']:
#                         xyz = tuple(atom.coord)
#                         mol_rows.append({'model':model_number,
#                                         'chain':chain_id,
#                                         'resname':resname,     
#                                         'resid':resno, 
#                                         'atom':atom.name, 
#                                         'X':xyz[0],
#                                         'Y':xyz[1],
#                                         'Z':xyz[2]})
#     df = pd.DataFrame(mol_rows)
#     for c in ['X','Y','Z']:
#         df[c] = df[c].astype('float32')
#     return df
#     #
#     #return model_to_chain_to_resno_atom_to_vals

def  bytecheck(fileblock):
    """Return ``fileblock`` as text, decoding bytes-like input with UTF-8 (if 
    not UTF-8 already).

    Args:
        fileblock: Raw PDB content as ``str``, ``bytes``, or ``bytearray``.

    Returns:
        str: Decoded text when given bytes-like input; unchanged string otherwise.
    """
    if isinstance(fileblock, (bytes, bytearray)):
        fileblock = fileblock.decode()
    return fileblock

def parse_PDB_lines(pdb_block):
    """Convert a PDB text block into a list of atom dictionaries.
    File parsing based on: 
    https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    
    Args:
        pdb_block (str): Raw PDB contents containing ``ATOM`` records.

    Returns:
        list[dict]: One dict per atom with parsed fields ``atomno``, ``atomtype``,
        ``resname``, ``chainname``, ``resno``, ``x``, ``y``, ``z``, ``occupancy``,
        ``tf``, and ``segname``.
    """
    
    # Following the numbering from: 
    # 
    atom_dictionary_list = []
    for line in pdb_block.split('\n'):
        line = line.lstrip()
        len_line = len(line)
        if len_line and line[:4] == 'ATOM':
            d = {
                'atomno':int(line[4:11].lstrip().rstrip()),
                'atom':  line[11:16].lstrip().rstrip(),
                'resname':   line[17:21].lstrip().rstrip(),
                'chain': line[21],
                'resid': int(line[22:26].lstrip().rstrip()),
                'X':   float(line[28:38].lstrip().rstrip()),
                'Y':   float(line[38:46].lstrip().rstrip()),
                'Z':   float(line[46:54].lstrip().rstrip()),
            }
            occupancy = None
            if len_line >= 54:
                occupancy = line[55:60].lstrip().rstrip()
                if len(occupancy):
                    occupancy = float(occupancy)
                else:
                    occupancy = None
                            
            tf = None
            if len_line >= 60:
                tf = line[61:66].lstrip().rstrip()
                if len(tf):
                    tf = float(tf)
                else:
                    tf = None
            segname = None
            if len_line >= 72:
                segname = line[72:76].lstrip().rstrip()
                if not len(segname):
                    segname = None
            #
            d['occupancy'] = occupancy
            d['tf']        = tf
            d['segname']   = segname

            atom_dictionary_list.append(d)
    #
    return atom_dictionary_list

# OLD VERSION (IN HOUSE). IT IS FASTER THAN THE CURRENT "read_pdb", WHICH IS BIOPDB RUN, BUT IT IS NOT 
# AS WELL TESTED.
def read_pdb_inhouse(fn_or_filehandle:Union[str,os.PathLike]):
    """Parse a PDB file without BioPython and extract backbone atoms.

    Args:
        fn: Path to a PDB file or a readable file-like object understood by
            ``get_filename_and_filehandle``.

    Returns:
        pandas.DataFrame: Backbone atoms (N, CA, C) for each model, chain, and
        residue with columns ``['model','chain','resname','resid','atom','X','Y','Z']``
        where coordinates are float32.
    """
    # print([fn_or_filehandle]) 
    fn, f = get_filename_and_filehandle(filename_or_filehandle=fn_or_filehandle)
    # f = open(fn,"r")
    pdbblock = f.read()
    f.close()
    
    pdbblock = bytecheck(pdbblock)
        
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
        model_block = models[model_index]
        if len(model_block.rstrip().lstrip()) > 1:
            model_number = model_index
            
            """
            ATOM     10 1H   LYS A   1       0.763   3.548  -0.564
            ATOM    482  N   PRO A  61      27.194  -5.761  14.684  1.00  9.09           N  
            ATOM      2  CA  BLYSX   1     -77.937 -26.325   6.934  1.00  0.00      U1    
            ATOM      3  HT2 MET U   1      -1.052  -0.551 -12.281  0.00  0.00      UBIQ  
                      |   |   |  |   |        |       |       |                     |
                 atomno   |   |  |   |        x       y       z                 segname
                   atom type  |  |   |                                          (CHAIN)
                        restype  |   resno
                                chainID
            """
            
            all_atom_mol_rows = parse_PDB_lines(model_block)
            backbone_mol_rows = []
            for row_dict in all_atom_mol_rows:
                if row_dict['atom'] in ["CA", "N", "C"]:
                    row_dict['model'] = model_number
                    backbone_mol_rows.append(row_dict)
            mol_rows += backbone_mol_rows
            #
        #
    #
    df = pd.DataFrame(mol_rows)
    for c in ['X','Y','Z']:
        df[c] = df[c].astype('float32')
    #
    return df
#

# Does a sniff test for whether this is a legitimate pdb file
def check_pdb(fn:Union[str,os.PathLike]):
    """Heuristically validate that a PDB file is safe to parse with Biopython.

    The check catches common CHARMM/PDB quirks (four-letter residue names shifting
    into the chain column or missing/placeholder chain identifiers). Any
    potential issue marks the file as problematic.

    Args:
        fn (str): Path to the PDB file to examine.

    Returns:
        int: ``1`` when the file appears well-formed, otherwise ``0``.
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

def calculate_dihedral_angle(p:np.ndarray):
    """Compute the dihedral angle defined by four 3D points. 
    
    E.g., providing positions for four contiguous backbone `(C-)(N)(CA)(C)` 
    would give you the ``phi`` backbone dihedral angle.
    Similarly, providing positions for `(N)(CA)(C)(N+)` would give you the 
    ``psi`` dihedral angle.

    Args:
        p (np.ndarray): Array of four 3D coordinates with shape ``(4, 3)`` representing 
        four contigusous backbone atoms.

    Returns:
        float: Signed dihedral angle in degrees in the range ``[-180, 180]``.
    """
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
    """Determine whether an object behaves like a file handle.

    Args:
        x (Any): Object to inspect.

    Returns:
        bool: True if ``x`` is an instance of ``io.IOBase``, otherwise False.
    """
    return isinstance(x, (io.IOBase, tarfile.TarFile, gzip.GzipFile, zipfile.ZipFile))


def get_pdb_filenames(pdbfn_or_dir:str):
    """Collect `.pdb` filenames from a file or directory and sort them numerically.
    
    Args:
        pdbfn_or_dir (str): Path to a single PDB file or a directory containing
            PDB files.

    Returns:
        tuple[list[str], str]: Sorted list of absolute `.pdb` file paths and the
            directory path they reside in. If the path is neither a file nor a
            directory, an empty list and the derived directory path are
            returned.
    """
    
    pdbfn = os.path.abspath(pdbfn_or_dir)
    pdbdir = os.path.dirname(pdbfn)
    list_pdbfilenames = []
    #
    if os.path.isfile(pdbfn):
        # then this pathname leads to a FILE
        # ... so keep as is
        list_pdbfilenames = [pdbfn]
        name = re.split(r'[\/\.]',pdbfn)[-2]
    elif os.path.isdir(pdbfn):
        pdbdir = pdbfn
        list_pdbfilenames = sorted(glob.glob(pdbdir+"/*.pdb.*"))
        name = re.split(r'[\/\.]',pdbfn)[-1]
    else:
        #print(helpme)
        print("Either filename or directory expected. Exiting.")
        return list_pdbfilenames, pdbdir
    #
    # JUST "CLEVERLY" ARRANGING THE FILENAMES, IF WE HAVE A SET OF FILENAMES RATHER THAN ONE
    # (e.g., list_pdbfilenames = [something2part1,something1part2,something1part1,something10part1]
    # list_pdbfilenames.sort() this list to: [something1part1,something1part2,something2part1,something10part1]
    REXP = re.compile( r'\d+' )
    def key_function( s ): return list(map(int, re.findall(REXP, s )))
    list_pdbfilenames.sort(key=key_function)
    #
    return list_pdbfilenames, pdbdir


def return_tgz_filehandle(pdbfn):
    """Return file-like handles for every regular file in a tar.gz or tgz archive.

    Args:
        pdbfn (str | os.PathLike): Path to a ``.tar.gz``/``.tgz`` archive.

    Returns:
        io.BufferedReader: Handle to the first file found in the archive.
    """
    file_objects = []
    # Open the tar archive in read mode
    tar = tarfile.open(pdbfn, 'r:gz')
    # Iterate through the members of the archive
    for member in tar.getmembers():
        if member.isfile():  # Check if it's a regular file
            # Get a file-like object for the member
            fileobj = tar.extractfile(member)
            if fileobj:  # Ensure the file object is not None
                file_objects.append(fileobj)
                # # Read the content from the file-like object
                #content = fileobj.read()
    #first_filehandle_in_list = file_objects[0]
    return file_objects

def return_zip_file_handle(pdbfn, mode='r'):
    """Open the first file inside a zip archive and return a readable handle.
    
    Args:
        pdbfn (str | os.PathLike): Path to the `.zip` archive containing a PDB file.
        mode (str, optional): Mode passed to ``ZipFile.open`` (e.g., ``'r'`` or ``'rb'``).
            Defaults to ``'r'``.

    Returns:
        ZipExtFile: Handle to the first file found in the archive.
    """
    with ZipFile(pdbfn, 'r') as zip_archive:
        all_files_in_archive  = zip_archive.namelist()
        all_handles_in_archive = [ zip_archive.open(file, mode=mode) 
                                    for file in all_files_in_archive]
        first_file_in_archive = all_handles_in_archive[0]
        #
        return all_handles_in_archive
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
    return [gzopen(pdb_fn,mode)]
#
def return_normal_filehandle(pdb,mode='rt'):
    """Open an uncompressed PDB file and return it as a single-item list.

    Args:
        pdb (str | os.PathLike): Path to the PDB file to open.
        mode (str, optional): Mode passed to ``open`` (e.g., ``'rt'`` for text).
            Defaults to ``'rt'``.

    Returns:
        list[io.TextIOWrapper]: List containing the opened file handle.
    """
    return [open(pdb,mode)]

# A dictionary that returns the appropriate filhandle given the file extension
dict_extention_to_open_object = {
    '.pdb.tar.gz':return_tgz_filehandle,
    '.pdb.tgz':return_tgz_filehandle,
    '.pdb.zip':return_zip_file_handle,
    '.pdb.gz':return_gz_filehandle,
    '.pdb':return_normal_filehandle,
}

def return_filehandle_from_gz_zip_or_normal_file(pdbfn):
    """Return a readable handle for `.pdb`, `.pdb.gz`, or `.pdb.zip` inputs.

    Args:
        pdbfn (str): Path to a PDB file, optionally gzipped or zipped. For zip
            archives, the first file in the archive is opened.

    Returns:
        IOBase: File-like object opened with the appropriate handler for the
            given extension, or ``None`` if the extension is unsupported.
    """
    #
    for file_extension in dict_extention_to_open_object.keys():
        if pdbfn[-len(file_extension):].lower() == file_extension:
            return dict_extention_to_open_object[file_extension](pdbfn)

def get_filename_and_filehandle(filename_or_filehandle):
    """Return lists of filenames and corresponding file handles.
    
    Accepts either a path to a PDB file (plain or compressed) or an existing
    file-like object. Compressed archives (.gz, .zip, .tgz/.tar.gz) are expanded
    into multiple handles, so filenames and handles are always returned as
    parallel lists.
    
    Args:
        filename_or_filehandle (str | io.IOBase): Path to a PDB-like file or an
            already opened file-like object.

    Returns:
        tuple[str, io.IOBase]: Filenames and readable handles.

    Raises:
        TypeError: If input is neither a filepath string nor a file-like object.
    """
    #
    # In order to accomodate .zip and .tgz files, we must report lists of files and filehandles
    fn = None
    f  = None
    if is_filehandle(filename_or_filehandle):
        fn = filename_or_filehandle.name
        f  = filename_or_filehandle
    elif isinstance(filename_or_filehandle, str) and not isinstance(filename_or_filehandle, list):
        list_f = return_filehandle_from_gz_zip_or_normal_file(filename_or_filehandle)
        f  = list_f[0]
        # deriving the filenames directly from the filehandle
        fn = f.name
    else:
        raise TypeError("fn_or_filehandle must be a filepath string or file-like object")
    #
    return fn, f

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

def median_filter(vals,nearest_neighbors=1):
    """Apply a 1D median filter, leaving edge elements unchanged. 

    Args:
        vals (Sequence[float]): Input values to smooth.
        nearest_neighbors (int, optional): Number of neighbors to include on
            each side when computing the median. Defaults to 1.

    Returns:
        list[float]: Filtered values where interior points are replaced by the
            median of their neighborhood and edge points are left as-is.
    """
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
    """Set the axis data aspect ratio using current plot limits.

    Args:
        aspect (float): Desired width-to-height ratio for the data units.
                        (E.g., if ``aspect == 2``, then the x axis will have twice the length of the y)
        ax (matplotlib.axes.Axes, optional): Axis to update. Defaults to the
            current axes.
    """
    if not ax: ax=plt.gca()
    extent = plt.axis()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

def get_coord_by_key_value_match(key,value,list_of_dict):
    """Fetch the first matching record's coordinates for a given key/value pair.

    Args:
        key (str): Dictionary key to compare against ``value``.
        value: Target value to match within each dictionary.
        list_of_dict (list[dict]): Sequence of records expected to include ``'X'``, ``'Y'``, and ``'Z'`` entries.

    Returns:
        list[float]: Coordinates ``[X, Y, Z]`` of the first matching record, or an empty list if none match.
    """
    valid_records = [d for d in list_of_dict if d[key]==value]
    coordinates = []
    if len(valid_records):
        record_to_use = valid_records[0]
        coordinates = [record_to_use['X'],record_to_use['Y'],record_to_use['Z']]
    
    return coordinates
