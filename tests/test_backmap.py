#!/usr/bin/env python
import pytest, sys, glob
#from matplotlib.figure import Figure # to assert the returning of a proper figure
import matplotlib.pyplot as plt
import matplotlib as mpl # currently imported to compare returned axes object to mpl.axes.Axes
import matplotlib.figure

import backmap

"""Tests for `backmap` package."""

__author__ = "ranjanmannige"
__copyright__ = "Ranjan Mannige"
__license__ = "MIT"

# # @pytest.fixture
# # def response():
# #     """Sample pytest fixture.

# #     See more at: http://doc.pytest.org/en/latest/fixture.html
# #     """
# #     # import requests
# #     # return requests.get('https://github.com/audreyfeldroy/cookiecutter-pypackage')


# # def test_content(response):
# #     """Sample pytest test function with the pytest fixture as an argument."""
# #     # from bs4 import BeautifulSoup
# #     # assert 'GitHub' in BeautifulSoup(response.content).title.string

def test_file_opener():
    """ 
    Testing the custom file opening function that opens .gz, .zip, and plain PDBs
    """
    file_blocks = []
    #
    # Three filenames that contain different compression 
    # types (zip, gz, and no compression)
    # List obtained by: glob.glob('backmap/tests/pdbs/1mba.pdb*')
    same_pdb_content_but_different_compressions = [ 
        'tests/pdbs/1mba.pdb.zip',
        'tests/pdbs/1mba.pdb',
        'tests/pdbs/1mba.pdb.gz',
        'tests/pdbs/1mba.pdb.tar.gz',
        'tests/pdbs/1mba.pdb.tgz'
    ]
    # TEST: opening files of different types
    for pdbfn in same_pdb_content_but_different_compressions:
        fille_list = backmap.utils.return_filehandle_from_gz_zip_or_normal_file(pdbfn)
        file_block = fille_list[0].read()
        #
        if isinstance(file_block, (bytes, bytearray)):
            file_block = file_block.decode()
        #
        file_blocks.append(file_block)
    #
    assert len(set(file_blocks))==1

    # TEST: parsing files of different types
    df_collection = []
    for pdbfn in glob.glob('backmap/tests/pdbs/1mba.pdb*'):
        print(pdbfn)
        #f = backmap.utils.return_filehandle_from_gz_zip_or_normal_file(pdbfn)
        #file_block = f.read()
        #if isinstance(file_block, (bytes, bytearray)):
        #    file_block = file_block.decode()
        df = backmap.process_PDB(pdbfn)
        
        df_collection.append(df)

    #[for i in range(len(df_collection))]
    all_comparison_indices = [(i,j) for i in range(0, len(df_collection)) 
                            for j in range(i+1, len(df_collection))]
    assert all([df_collection[i].equals(df_collection[j]) for i,j in all_comparison_indices])

def test_ramachandran_number():
    # Testing the basic calculations
    assert backmap.R( 180, 180) == 1.0
    assert backmap.R(-180,-180) == 0.0
    assert backmap.R(   0,   0) == 0.5

def test_load_pdb():
    # Set pdb name 
    pdbfn = 'tests/pdbs/1mba.pdb'
    
    # TEST: Known features of the test file are queried here
    pdb_df = backmap.utils.read_pdb(filename_or_filehandle=pdbfn)
    assert set(pdb_df['model']) == {0}
    assert set(pdb_df['chain']) == {'A'}
    assert min(pdb_df['resid']) == 1
    assert max(pdb_df['resid']) == 146
#

def test_processing_data():
    pdbfn = 'tests/pdbs/1mba.pdb'
    structure_df = backmap.process_PDB(pdbfn=pdbfn, signed=False)
    assert structure_df.shape[0] == 146
    #
    should_be_true, figure_dict = backmap.draw_figures(structure_df=structure_df, 
                                   output_dir='', write=False, show=False)
    # Basic test
    assert should_be_true
    # Check if the first and last keys of figure_dict are of 
    # type matplotlib.figure.Figure
    fig_keys = list(figure_dict.keys())
    possible_fig = figure_dict[fig_keys[0]]
    assert isinstance(possible_fig, matplotlib.figure.Figure)
    possible_fig = figure_dict[fig_keys[-1]]
    assert isinstance(possible_fig, matplotlib.figure.Figure)
