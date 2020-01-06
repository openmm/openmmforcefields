.. openmmforcefields documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


AMBER and CHARMM force fields for OpenMM
========================================

This repository provides support for AMBER and CHARMM force fields and small molecule parameterization with GAFF and the Open Force Field Toolkit for OpenMM.

Supported force fields
----------------------

* **AMBER:** All major AMBER force fields distributed with `AmberTools <https://ambermd.org/AmberTools.php>`_ 19.9 (except ff19SB), as well as all released `GAFF small molecule force fields <http://ambermd.org/antechamber/gaff.html>`_ through 1.81 (GAFF 1.x) and 2.11 (GAFF 2.x).
* **CHARMM:** Non-polarizable protein, nucleic acid, and pre-parameterized small molecule force fields available in in the [Aug 2015 CHARMM36 force field release from the Mackerell website](http://mackerell.umaryland.edu/charmm_ff.shtml). *Note that this conversion has not yet been fully validated.
* **Open Force Field Initiative force fields:** All distributed `Open Force Field Initiative <http://openforcefield.org>`_ force fields, including the `smirnoff99Frosst` series and `openforcefield-1.0.0 ("Parsley") <https://openforcefield.org/news/introducing-openforcefield-1.0/>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   generators


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
