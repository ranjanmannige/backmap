.. image:: https://github.com/ranjanmannige/backmap/raw/main/docs/images/banner.jpg
    :alt: BackMAP Banner
    :width: 75%

|pypi| |license| |pyversions| |status| |downloads|

.. |pypi| image:: https://img.shields.io/pypi/v/backmap.svg?label=version
    :target: http://pypi.org/project/backmap
    :alt: PyPi Release Version 

.. |format| image:: https://img.shields.io/pypi/format/Backmap.svg
    :alt: Format

.. |license| image:: https://img.shields.io/pypi/l/Backmap.svg
    :target: https://github.com/ranjanmannige/backmap/blob/master/LICENSE
    :alt: License

.. |downloads| image:: https://static.pepy.tech/personalized-badge/backmap?period=total&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=BLUE&left_text=downloads
    :alt: Github All Releases

.. |status| image:: https://img.shields.io/pypi/status/Backmap.svg
    :alt: Status

.. |pyversions| image:: https://img.shields.io/pypi/pyversions/Backmap.svg
    :alt: Allowed python environments_current_filenames

.. contents::


Full documentation available here: https://ranjanmannige.github.io/backmap/.


Introduction
============

BackMAP is a Python module (and stand-alone tool) that helps with the visualization of large amounts of structural (space-time) backbone data in a single graph. It utilizes a new per-residue backbone metric -- the Ramachandran number -- to provide easily readable "pictures" (multi-angle pictures or MAPS) of protein conformations, ensembles, and trajectories. Input structures can be either a combined protein databank (PDB) structure file, or a directory of such files, and produces graphs.


Installation
============

PIP Installation
-----------------

Running the following at a command prompt (terminal) would get the job done (the '-I' is not necessary, but ensures the latest sub-version is installed):

.. code-block:: bash

    $ pip install -I backmap


GIT Installation
----------------


.. code-block:: bash

    $ git clone https://github.com/ranjanmannige/backmap.git
    $ cd backmap
    $ pip install .
    # For testing
    $ pip install pytest
    $ pytest


Manual Installation
-------------------

Manually download the source code (`main.zip <https://github.com/ranjanmannige/backmap/archive/refs/heads/main.zip>`_) from the `git repository <https://github.com/ranjanmannige/backmap>`_. Then, same as above:

.. code-block:: bash
    
    # In stead of downloading, you can follow the next two commands (tested only on linux)
    $ wget https://github.com/ranjanmannige/backmap/archive/refs/heads/main.zip
    $ unzip main.zip # Should giv you a directory called "backmap-main"
    # The rest is the same as with installing using `git clone`
    $ cd backmap-main
    $ pip install .
    # For testin
    $ pip install pytest
    $ pytest


Usage
=====

In-script usage
---------------

.. code-block:: python

    import backmap
    print backmap.R(phi=0,psi=0) # Should print '0.5'

For more information about in-script module usage, refer to the `manuscript <https://raw.githubusercontent.com/ranjanmannige/backmap/master/manuscript/manuscript/backmap.pdf>`_ associated with this module.

Standalone usage
----------------

After installation, the following commands produce a variety of graphs (exampled below).

.. code-block:: bash

    $ python -m backmap -pdb ./pdbs/ProteinDatabankStructureFilename.pdb
    $ python -m backmap -pdb /directory/containing/pdbs/
    

Examples
========

Example 1: A stable protein (`1xqq <https://www.rcsb.org/structure/1XQQ>`_)
------------------------------------------------------------------------------

The Panels **(b)** through **(f)** were created by running the following command within thin the downloaded directory (Panel **(a)** was created using `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_).

.. code-block:: bash

    $ python -m backmap -pdb ./tests/pdbs/1xqq.pdb

As evident below, the graphs generated from the protein ensemble `1xqq <https://www.rcsb.org/structure/1XQQ>`_ describes a conformationally stable protein (each graph is detailed below). 

.. image:: https://github.com/ranjanmannige/backmap/raw/main/docs/images/fig-11-2x.jpg

Each column in Panel **(b)** describes the histogram in Ramachandran number (R) space for a single model/timeframe. These histograms show the presence of both helices (at R \~ 0.34) and sheets (at R \~ 0.52). Additionally, Panels **(c)** and **(d)** describe the per-residue conformational plots (colored by two different metrics or CMAPs), which show that most of the protein backbone remains relatively stable (e.g., few fluctuations in state or 'color' are evident over the frame \#). Finally, Panel **(e)** describes the extent towards which a single residue's state has deviated from the first frame, and Panel **(f)** describes the extent towards which a single residue's state has deviated from its state in the previous frame. Both these graphs, as expected from Panels **(c)** and **(d)**, show that this protein is relatively conformationally stable.


Example 2: An intrinsically disordered protein (`2fft <https://www.rcsb.org/structure/2FFT>`_)
----------------------------------------------------------------------------------------------

As compared to the conformationally stable protein above, an intrinsically disordered protein (`2fft <https://www.rcsb.org/structure/2FFT>`_)
is much more flexible

.. image:: https://github.com/ranjanmannige/backmap/raw/main/docs/images/fig-12-2x.jpg

Panel **(b)** shows that the states accessed per model are diverse and dramatically fluctuate over the entire range of R (this is especially true when compared to a stable protein, see above). 

The diverse states occupied by each residue (Panels **(c)** and **(d)**) confirm this conformational variation within most residues (Panels **(e)** and **(f)** similarly show how most of the residues fluctuate dramatically).

Yet, interestingly, Panels **(c)** through **(f)** also show an unusually stable region -- residues 15 through 25 -- which consistently display the same conformational (alpha-helical) state at R \~ 0.33 (interpreted as the color red in Panel **(c)**). This trend would be hard to recognize by simply looking at the structure (Panel **(a)**). 

.. _publications:

Publications
============
The Ramachandran number concept is discussed in the following manuscripts (this tool is discussed in the first reference):

1. Mannige (2018) "The Backmap Python Module: How a Simpler Ramachandran Number Can Simplify the Life of a Protein Simulator" PeerJ 6:e5745 [`PeerJ Journal Link <https://doi.org/10.7717/peerj.5745>`_].

2. Mannige, Kundu, Whitelam (2016) "The Ramachandran Number: An Order Parameter for Protein Geometry" PLoS ONE 11(8): e0160023 [`PLoS ONE Journal Link <https://doi.org/10.1371/journal.pone.0160023>`_].
