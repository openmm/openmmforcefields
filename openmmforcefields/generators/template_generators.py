"""
Residue template generators for GAFF, SMIRNOFF, and espaloma small molecule force fields.
"""

################################################################################
# IMPORTS
################################################################################

import contextlib
import logging
import os
import warnings

from ..utils import classproperty

################################################################################
# LOGGER
################################################################################

_logger = logging.getLogger("openmmforcefields.generators.template_generators")

################################################################################
# Small molecule OpenMM ForceField template generation utilities
################################################################################


class ForceException(Exception):
    """Exception for forces"""


class SmallMoleculeTemplateGenerator:
    """
    Abstract base class for small molecule template generation for OpenMM ForceField.

    This class should not be used directly, but provides utility routines for
    subclasses that generate small molecule residue templates via external tools.

    Parameters
    ----------
    debug_ffxml_filename : str
        If not None, the generated ffxml file will be written to this filename.
        Default is None.
    """

    def __init__(self, molecules=None, cache=None):
        """
        Create a tempalte generator with some OpenFF toolkit molecules

        Requies the openff-toolkit toolkit: http://openforcefield.org

        Parameters
        ----------
        molecules : openff.toolkit.Molecule or list, optional, default=None
            Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can be used
            to construct a Molecule.
            Can also be a list of Molecule objects or objects that can be used to construct a Molecule.
            If specified, these molecules will be recognized and parameterized as needed.
            The parameters will be cached in case they are encountered again the future.
        cache : str, optional, default=None
            Filename for global caching of parameters.
            If specified, parameterized molecules will be stored in a TinyDB instance as a JSON file.
        """

        # Store specified molecules
        self._molecules = dict()
        self.add_molecules(molecules)

        # Set up cache
        self._cache = cache
        self._smiles_added_to_db = set()  # set of SMILES added to the database this session
        self._database_table_name = None  # this must be set by subclasses for cache to function

        # Name of the force field
        self._forcefield = None  # this must be set by subclasses

        # File to write ffxml to if requested
        self.debug_ffxml_filename = None

    @property
    def forcefield(self):
        """The current force field name in use"""
        return self._forcefield

    @contextlib.contextmanager
    def _open_db(self):
        """Open the cache database."""
        from tinydb import TinyDB

        tinydb_kwargs = {
            "sort_keys": True,
            "indent": 4,
            "separators": (",", ": "),
        }  # for pretty-printing
        db = TinyDB(self._cache, **tinydb_kwargs)
        try:
            yield db
        finally:
            db.close()

    def add_molecules(self, molecules=None):
        """
        Add specified list of Molecule objects to cached molecules that will be recognized.

        Parameters
        ----------
        molecules : openff.toolkit.Molecule or list of Molecules, optional, default=None
            Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can be used
            to construct a Molecule.
            Can also be a list of Molecule objects or objects that can be used to construct a Molecule.
            If specified, these molecules will be recognized and parameterized as needed.
            The parameters will be cached in case they are encountered again the future.

        Examples
        --------

        Add some Molecules later on after the generator has been registered:

        >>> from openff.toolkit import Molecule
        >>> smiles = ('c1ccccc1', 'O=Cc1ccc(O)c(OC)c1', 'CN1CCC[C@H]1c2cccnc2')
        >>> mol1, mol2, mol3 = [Molecule.from_smiles(s, allow_undefined_stereo=True) for s in smiles]
        >>> generator.add_molecules(mol1)  # doctest: +SKIP
        >>> generator.add_molecules([mol2, mol3])  # doctest: +SKIP
        """

        # Return if empty
        if not molecules:
            return

        # Ensure molecules is an iterable
        try:
            iter(molecules)
        except TypeError:
            molecules = [molecules]

        # Create copies
        # TODO: Do we need to try to construct molecules with other schemes, such as Molecule.from_smiles(), if needed?
        import copy

        molecules = [copy.deepcopy(molecule) for molecule in molecules]

        # Cache molecules
        self._molecules.update({molecule.to_smiles(): molecule for molecule in molecules})

    @staticmethod
    def _match_residue(residue, molecule_template):
        """
        Determine whether a residue matches a Molecule template and return a list of corresponding atoms.

        This implementation uses NetworkX for graph isomorphism determination.

        Parameters
        ----------
        residue : openmm.app.topology.Residue
            The residue to check
        molecule_template : openff.toolkit.Molecule
            The Molecule template to compare it to

        Returns
        -------
        matches : dict of int : int
            matches[residue_atom_index] is the corresponding Molecule template atom index
            or None if it does not match the template

        .. todo :: Can this be replaced by an isomorphism matching call to the OpenFF toolkit?
        """

        # TODO: Speed this up by rejecting molecules that do not have the same chemical formula

        # TODO: Can this NetworkX implementation be replaced by an isomorphism
        # matching call to the OpenFF toolkit?

        import networkx as nx

        # Build list of external bonds for residue
        number_of_external_bonds = {atom: 0 for atom in residue.atoms()}
        for bond in residue.external_bonds():
            if bond[0] in number_of_external_bonds:
                number_of_external_bonds[bond[0]] += 1
            if bond[1] in number_of_external_bonds:
                number_of_external_bonds[bond[1]] += 1

        # Residue graph
        residue_graph = nx.Graph()
        for atom in residue.atoms():
            residue_graph.add_node(
                atom,
                element=atom.element.atomic_number,
                number_of_external_bonds=number_of_external_bonds[atom],
            )
        for bond in residue.internal_bonds():
            residue_graph.add_edge(bond[0], bond[1])

        # Template graph
        # TODO: We can support templates with "external" bonds or atoms using attached string data in future
        # See https://docs.eyesopen.com/toolkits/python/oechemtk/OEChemClasses/OEAtomBase.html
        template_graph = nx.Graph()
        for atom_index, atom in enumerate(molecule_template.atoms):
            template_graph.add_node(atom_index, element=atom.atomic_number, number_of_external_bonds=0)
        for bond in molecule_template.bonds:
            template_graph.add_edge(bond.atom1_index, bond.atom2_index)

        # DEBUG
        # print(f'residue_graph: nodes {len(list(residue_graph.nodes))} edges {len(list(residue_graph.edges))}')
        # print(f'template_graph: nodes {len(list(template_graph.nodes))} edges {len(list(template_graph.edges))}')

        # Determine graph isomorphism
        from networkx.algorithms import isomorphism

        def node_match(n1, n2):
            """Return True if nodes match, and False if not"""
            return (n1["element"] == n2["element"]) and (
                n1["number_of_external_bonds"] == n2["number_of_external_bonds"]
            )

        graph_matcher = isomorphism.GraphMatcher(residue_graph, template_graph, node_match=node_match)
        if not graph_matcher.is_isomorphic():
            return None

        # Translate to local residue atom indices
        atom_index_within_residue = {atom: index for (index, atom) in enumerate(residue.atoms())}
        matches = {
            atom_index_within_residue[residue_atom]: template_atom
            for (residue_atom, template_atom) in graph_matcher.mapping.items()
        }

        return matches

    def _molecule_has_user_charges(self, molecule):
        """
        Checks whether or not a molecule has user charges (assigned and not all
        approximately equal to zero).  If so, checks to ensure that their sum
        matches the formal charge of the molecule, and issues a warning if not.

        Returns
        -------
        result : bool
            True if user charges are assigned and are not all approximately
            equal to zero; False otherwise.
        """

        import numpy as np
        from openff.units import unit

        if molecule.partial_charges is None:
            return False

        # We should be able to trust that charges from an OpenFF molecule will
        # be in OpenFF units, not OpenMM units or dimensionless values: see
        # https://github.com/openforcefield/openff-toolkit/pull/1619.

        partial_charges = molecule.partial_charges.m_as(unit.elementary_charge)
        if np.allclose(partial_charges, 0):
            return False

        actual_sum = partial_charges.sum()
        expected_sum = molecule.total_charge.m_as(unit.elementary_charge)
        if not np.isclose(actual_sum, expected_sum):
            warnings.warn(
                f"Sum of user-provided partial charges {actual_sum} does not match formal charge {expected_sum}"
            )

        return True

    def _generate_unique_atom_names(self, molecule):
        """
        Generate unique atom names

        Parameters
        ----------
        molecule : openff.toolkit.Molecule
            The molecule whose atom names are to be modified in-place
        """

        from collections import defaultdict

        element_counts = defaultdict(int)
        for atom in molecule.atoms:
            symbol = atom.symbol
            element_counts[symbol] += 1
            atom.name = symbol + str(element_counts[symbol])

    def generator(self, forcefield, residue):
        """
        Residue template generator method to register with openmm.app.ForceField

        Parameters
        ----------
        forcefield : openmm.app.ForceField
            The ForceField object to which residue templates and/or parameters are to be added.
        residue : openmm.app.Topology.Residue
            The residue topology for which a template is to be generated.

        Returns
        -------
        success : bool
            If the generator is able to successfully parameterize the residue, `True` is returned.
            If the generator cannot parameterize the residue, it should return `False` and not modify `forcefield`.
        """

        if self._database_table_name is None:
            raise NotImplementedError(
                "SmallMoleculeTemplateGenerator is an abstract base class and cannot be used directly."
            )

        from io import StringIO

        # TODO: Refactor to reduce code duplication

        _logger.info(f"Requested to generate parameters for residue {residue}")

        # If a database is specified, check against molecules in the database
        if self._cache is not None:
            with self._open_db() as db:
                table = db.table(self._database_table_name)
                for entry in table:
                    # Skip any molecules we've added to the database this session
                    if entry["smiles"] in self._smiles_added_to_db:
                        continue

                    # See if the template matches
                    from openff.toolkit import Molecule

                    molecule_template = Molecule.from_smiles(entry["smiles"], allow_undefined_stereo=True)
                    _logger.debug(f"Checking against {entry['smiles']}")
                    if self._match_residue(residue, molecule_template):
                        ffxml_contents = entry["ffxml"]

                        # Write to debug file if requested
                        if self.debug_ffxml_filename is not None:
                            with open(self.debug_ffxml_filename, "w") as outfile:
                                _logger.debug(f"writing ffxml to {self.debug_ffxml_filename}")
                                outfile.write(ffxml_contents)

                        # Add parameters and residue template for this residue
                        forcefield.loadFile(StringIO(ffxml_contents))
                        # Signal success
                        return True

        # Check against the molecules we know about
        for smiles, molecule in self._molecules.items():
            # See if the template matches
            if self._match_residue(residue, molecule):
                # Generate template and parameters.
                ffxml_contents = self.generate_residue_template(molecule)

                # Write to debug file if requested
                if self.debug_ffxml_filename is not None:
                    with open(self.debug_ffxml_filename, "w") as outfile:
                        _logger.debug(f"writing ffxml to {self.debug_ffxml_filename}")
                        outfile.write(ffxml_contents)

                # Add the parameters and residue definition
                forcefield.loadFile(StringIO(ffxml_contents))
                # If a cache is specified, add this molecule
                if self._cache is not None:
                    with self._open_db() as db:
                        table = db.table(self._database_table_name)
                        _logger.debug(f"Writing residue template for {smiles} to cache {self._cache}")
                        record = {"smiles": smiles, "ffxml": ffxml_contents}
                        # Add the IUPAC name for convenience if we can
                        try:
                            record["iupac"] = molecule.to_iupac()
                        except Exception:
                            pass
                        # Store the record
                        table.insert(record)
                        self._smiles_added_to_db.add(smiles)

                # Signal success
                return True

        # Report that we have failed to parameterize the residue
        _logger.warning(
            f"Did not recognize residue {residue.name}; did you forget to call .add_molecules() to add it?"
        )
        return False


################################################################################
# GAFF-specific OpenMM ForceField template generation utilities
################################################################################


class GAFFTemplateGenerator(SmallMoleculeTemplateGenerator):
    """
    OpenMM ForceField residue template generator for GAFF/AM1-BCC using pre-cached OpenFF toolkit molecules.

    One template generator can be registered in multiple OpenMM ForceField objects.

    Examples
    --------

    Create a template generator for GAFF for a single Molecule and register it with ForceField:

    >>> # Define a Molecule using the OpenFF Molecule object
    >>> from openff.toolkit import Molecule
    >>> molecule = Molecule.from_smiles('c1ccccc1')
    >>> # Create the GAFF template generator
    >>> from openmmforcefields.generators import GAFFTemplateGenerator
    >>> template_generator = GAFFTemplateGenerator(molecules=molecule)
    >>> # Create an OpenMM ForceField
    >>> from openmm.app import ForceField
    >>> amber_forcefields = ['amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml']
    >>> forcefield = ForceField(*amber_forcefields)
    >>> # Register the template generator
    >>> forcefield.registerTemplateGenerator(template_generator.generator)

    Create a template generator for a specific GAFF version and register multiple molecules:

    >>> molecule1 = Molecule.from_smiles('c1ccccc1')
    >>> molecule2 = Molecule.from_smiles('CCO')
    >>> template_generator = GAFFTemplateGenerator(molecules=[molecule1, molecule2], forcefield='gaff-2.2.20')

    You can also add some Molecules later on after the generator has been registered:

    >>> template_generator.add_molecules(molecule)
    >>> template_generator.add_molecules([molecule1, molecule2])

    You can optionally create or use a tiny database cache of pre-parameterized molecules:

    >>> template_generator = GAFFTemplateGenerator(cache='gaff-molecules.json')

    Newly parameterized molecules will be written to the cache, saving time next time!

    You can see which force fields are supported with

    >>> template_generator.INSTALLED_FORCEFIELDS
    ['gaff-1.4', 'gaff-1.8', 'gaff-1.81', 'gaff-2.1', 'gaff-2.11', 'gaff-2.2.20']
    """

    INSTALLED_FORCEFIELDS = [
        "gaff-1.4",
        "gaff-1.8",
        "gaff-1.81",
        "gaff-2.1",
        "gaff-2.11",
        "gaff-2.2.20",
    ]

    def __init__(self, molecules=None, forcefield=None, cache=None, template_generator_kwargs=None):
        """
        Create a GAFFTemplateGenerator with some OpenFF toolkit molecules

        Requies the OpenFF toolkit: http://openforcefield.org

        Parameters
        ----------
        molecules : openff.toolkit.Molecule or list, optional, default=None
            Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can be used
            to construct a Molecule.
            Can also be a list of Molecule objects or objects that can be used to construct a Molecule.
            If specified, these molecules will be recognized and parameterized with antechamber as needed.
            The parameters will be cached in case they are encountered again the future.
        forcefield : str, optional, default=None
            GAFF force field to use, one of
            ['gaff-1.4', 'gaff-1.8', 'gaff-1.81', 'gaff-2.1', 'gaff-2.11', 'gaff-2.2.20']
            If not specified, the latest GAFF supported version is used.
            GAFFTemplateGenerator.INSTALLED_FORCEFIELDS contains a complete up-to-date list of supported force fields.
        cache : str, optional, default=None
            Filename for global caching of parameters.
            If specified, parameterized molecules will be stored in a TinyDB instance as a JSON file.
        template_generator_kwargs : dict, optional, default=None
            Additional parameters for the template generator (ignored by GAFFTemplateGenerator).

        Examples
        --------

        Create a GAFF template generator for a single molecule (benzene, created from SMILES) and register it:

        >>> from openff.toolkit import Molecule
        >>> molecule = Molecule.from_smiles('c1ccccc1')
        >>> from openmmforcefields.generators import GAFFTemplateGenerator
        >>> gaff = GAFFTemplateGenerator(molecules=molecule)
        >>> from openmm.app import ForceField
        >>> amber_forcefields = [
        ...     'amber/protein.ff14SB.xml',
        ...     'amber/tip3p_standard.xml',
        ...     'amber/tip3p_HFE_multivalent.xml',
        ... ]
        >>> forcefield = ForceField(*amber_forcefields)
        >>> forcefield.registerTemplateGenerator(gaff)

        The latest GAFF version is used if none is specified.
        You can check which GAFF version is in use with

        >>> gaff.forcefield
        'gaff-2.2.20'

        Create a template generator for a specific GAFF version for multiple molecules read from an SDF file:

        >>> molecules = Molecule.from_file('molecules.sdf')  # doctest: +SKIP
        >>> gaff = GAFFTemplateGenerator(molecules=molecules, forcefield='gaff-2.2.20')  # doctest: +SKIP

        You can also add molecules later on after the generator has been registered:

        >>> molecule = Molecule.from_smiles("CCO")
        >>> gaff.add_molecules(molecule)

        >>> molecule1, molecule2 = [Molecule.from_smiles(smiles) for smiles in ["CC", "c1ccccc1"]]
        >>> gaff.add_molecules([molecule1, molecule2])

        To check which GAFF versions are supported, check the `INSTALLED_FORCEFIELDS` attribute:

        >>> print(GAFFTemplateGenerator.INSTALLED_FORCEFIELDS)
        ['gaff-1.4', 'gaff-1.8', 'gaff-1.81', 'gaff-2.1', 'gaff-2.11', 'gaff-2.2.20']

        You can optionally create or use a tiny database cache of pre-parameterized molecules:

        >>> gaff = GAFFTemplateGenerator(cache='gaff-molecules.json', forcefield='gaff-2.2.20')

        Newly parameterized molecules will be written to the cache, saving time next time!
        """

        # Initialize molecules and cache
        super().__init__(molecules=molecules, cache=cache)

        # Use latest supported GAFF version if none is specified
        if forcefield is None:
            forcefield = self.INSTALLED_FORCEFIELDS[-1]

        # Ensure a valid GAFF version is specified
        if forcefield not in self.INSTALLED_FORCEFIELDS:
            raise ValueError(f"Specified 'forcefield' ({forcefield}) must be one of {self.INSTALLED_FORCEFIELDS}")

        # Store user-specified GAFF version
        self._forcefield = forcefield
        import re

        result = re.match(r"^gaff-(?P<major_version>\d+)\.(?P<minor_version>[\d.]+)$", forcefield)
        if result is None:
            msg = "'forcefield' must be of form 'gaff-X.Y', where X and Y denote major and minor version\n"
            msg += f"Provided 'forcefield' argument was '{forcefield}'\n"
            msg += f"Supported values are: {self.INSTALLED_FORCEFIELDS}"
            raise ValueError(msg)
        self._gaff_major_version = result["major_version"]
        self._gaff_minor_version = result["minor_version"]
        self._gaff_version = f"{self._gaff_major_version}.{self._gaff_minor_version}"

        # Track parameters by GAFF version string
        # TODO: Use file hash instead of name?
        import os

        self._database_table_name = os.path.basename(forcefield)

    @property
    def gaff_version(self):
        """The current GAFF version in use"""
        return self._gaff_major_version + "." + self._gaff_minor_version

    @property
    def gaff_major_version(self):
        """The current GAFF major version in use"""
        return self._gaff_major_version

    @property
    def gaff_minor_version(self):
        """The current GAFF minor version in use"""
        return self._gaff_minor_version

    @property
    def gaff_dat_filename(self):
        """File path to the GAFF .dat AMBER force field file"""
        import importlib_resources

        filename = str(
            importlib_resources.files("openmmforcefields")
            / "ffxml"
            / "amber"
            / "gaff"
            / "dat"
            / f"{self._forcefield}.dat"
        )
        return filename

    @property
    def gaff_xml_filename(self):
        """File path to the GAFF .ffxml OpenMM force field file"""
        import importlib_resources

        filename = str(
            importlib_resources.files("openmmforcefields")
            / "ffxml"
            / "amber"
            / "gaff"
            / "ffxml"
            / f"{self._forcefield}.xml"
        )
        return filename

    def generator(self, forcefield, residue):
        """
        Residue template generator method to register with openmm.app.ForceField

        Parameters
        ----------
        forcefield : openmm.app.ForceField
            The ForceField object to which residue templates and/or parameters are to be added.
        residue : openmm.app.Topology.Residue
            The residue topology for which a template is to be generated.

        Returns
        -------
        success : bool
            If the generator is able to successfully parameterize the residue, `True` is returned.
            If the generator cannot parameterize the residue, it should return `False` and not modify `forcefield`.
        """

        return super().generator(forcefield, residue)

    def generate_residue_template(self, molecule):
        """
        Generate a residue template and additional parameters for the specified Molecule.

        Parameters
        ----------
        molecule : openff.toolkit.Molecule
            A molecule to generate a residue template for. The molecule will be
            parameterized with antechamber as needed. The parameters will be
            cached in case it is encountered again the future.

        Returns
        -------
        ffxml_contents : str
            Contents of ForceField `ffxml` file containing additional parameters and residue template.

        Notes
        -----

        * The residue template will be named after the SMILES of the molecule.
        * This method preserves stereochemistry during AM1-BCC charge parameterization.
        * Atom names in molecules will be assigned Tripos atom names if any are blank or not unique.
        """

        from openff.units import unit

        # Use the canonical isomeric SMILES to uniquely name the template
        smiles = molecule.to_smiles()
        _logger.info(f"Generating a residue template for {smiles} using {self._forcefield}")

        # Generate unique atom names
        self._generate_unique_atom_names(molecule)

        # Compute partial charges if required
        if self._molecule_has_user_charges(molecule):
            _logger.debug("Using user-provided charges because partial charges are nonzero...")
        else:
            _logger.debug("Computing AM1-BCC charges...")
            molecule.assign_partial_charges(partial_charge_method="am1bcc", normalize_partial_charges=True)

        # Geneate a single conformation
        _logger.debug("Generating a conformer...")
        molecule.generate_conformers(n_conformers=1)

        # Create temporary directory for running antechamber
        import os
        import tempfile

        tmpdir = tempfile.mkdtemp()
        prefix = "molecule"
        input_sdf_filename = os.path.join(tmpdir, prefix + ".sdf")
        gaff_mol2_filename = os.path.join(tmpdir, prefix + ".gaff.mol2")
        frcmod_filename = os.path.join(tmpdir, prefix + ".frcmod")

        # Write MDL SDF file for input into antechamber
        molecule.to_file(input_sdf_filename, file_format="sdf")

        # Parameterize the molecule with antechamber (without charging)
        _logger.debug("Running antechamber...")
        self._run_antechamber(
            molecule_filename=input_sdf_filename,
            input_format="mdl",
            gaff_mol2_filename=gaff_mol2_filename,
            frcmod_filename=frcmod_filename,
        )

        # Read the resulting GAFF mol2 file atom types
        _logger.debug("Reading GAFF atom types...")
        self._read_gaff_atom_types_from_mol2(gaff_mol2_filename, molecule)

        # Generate additional parameters if needed
        _logger.debug("Creating ffxml contents for additional parameters...")
        from inspect import (  # use introspection to support multiple parmed versions
            signature,
        )
        from io import StringIO

        leaprc = StringIO(f"parm = loadamberparams {frcmod_filename}")
        import parmed

        params = parmed.amber.AmberParameterSet.from_leaprc(leaprc)
        kwargs = {}
        if "remediate_residues" in signature(parmed.openmm.OpenMMParameterSet.from_parameterset).parameters:
            kwargs["remediate_residues"] = False
        params = parmed.openmm.OpenMMParameterSet.from_parameterset(params, **kwargs)
        ffxml = StringIO()
        kwargs = {}
        if "write_unused" in signature(params.write).parameters:
            kwargs["write_unused"] = True

        params.write(ffxml, **kwargs)
        ffxml_contents = ffxml.getvalue()

        # Create the residue template
        _logger.debug("Creating residue template...")
        from lxml import etree

        root = etree.fromstring(ffxml_contents)
        # Create residue definitions
        residues = etree.SubElement(root, "Residues")
        residue = etree.SubElement(residues, "Residue", name=smiles)
        for atom in molecule.atoms:
            charge_string = str(atom.partial_charge.m_as(unit.elementary_charge))

            etree.SubElement(
                residue,
                "Atom",
                name=atom.name,
                type=atom.gaff_type,
                charge=charge_string,
            )

        for bond in molecule.bonds:
            etree.SubElement(
                residue,
                "Bond",
                atomName1=bond.atom1.name,
                atomName2=bond.atom2.name,
            )
        # Render XML into string and append to parameters
        ffxml_contents = etree.tostring(root, pretty_print=True, encoding="unicode")
        _logger.debug("ffxml creation complete.")

        return ffxml_contents

    def _run_antechamber(
        self,
        molecule_filename,
        input_format="sdf",
        gaff_mol2_filename=None,
        frcmod_filename=None,
        verbosity=0,
    ):
        """
        Run AmberTools antechamber and parmchk2 to create GAFF mol2 and frcmod files.

        Parameters
        ----------
        molecule_filename : str
            The molecule to be parameterized.
        input_format : str, optional, default='sdf'
            Format specifier for input file to pass to antechamber.
        gaff_mol2_filename : str, optional, default=None
            Name of GAFF mol2 filename to output.  If None, uses local directory
            and molecule_name
        frcmod_filename : str, optional, default=None
            Name of GAFF frcmod filename to output.  If None, uses local directory
            and molecule_name
        verbosity : int, default=0
            Verbosity for antechamber

        Returns
        -------
        gaff_mol2_filename : str
            GAFF format mol2 filename produced by antechamber containing GAFF 1/2 atom types
        frcmod_filename : str
            Amber frcmod file containing additional parameters for the molecule not found in corresponding gaff.dat
        """

        if gaff_mol2_filename is None:
            gaff_mol2_filename = "molecule.gaff.mol2"
        if frcmod_filename is None:
            frcmod_filename = "molecule.frcmod"

        # Build absolute paths for input and output files
        import os

        molecule_filename = os.path.abspath(molecule_filename)
        gaff_mol2_filename = os.path.abspath(gaff_mol2_filename)
        frcmod_filename = os.path.abspath(frcmod_filename)

        def read_file_contents(filename):
            infile = open(filename)
            contents = infile.read()
            infile.close()
            return contents

        # Use temporary directory context to do this to avoid issues with spaces in filenames, etc.
        import subprocess
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)

            local_input_filename = "in." + input_format
            import shutil

            shutil.copy(molecule_filename, local_input_filename)

            # Determine whether antechamber supports -dr [yes/no] option
            cmd = "antechamber -h | grep dr"
            supports_acdoctor = False
            if "acdoctor" in subprocess.getoutput(cmd):
                supports_acdoctor = True

            if self._gaff_major_version == "1":
                atom_type = "gaff"
            elif self._gaff_major_version == "2":
                atom_type = "gaff2"
            else:
                raise ValueError(f"gaff major version {self._gaff_major_version} unknown")

            # Run antechamber without charging (which is done separately)
            cmd = (
                f"antechamber -i {local_input_filename} -fi {input_format} "
                f"-o out.mol2 -fo mol2 -s {verbosity} -at {atom_type}"
            )
            if supports_acdoctor:
                cmd += " -dr " + ("yes" if verbosity else "no")

            _logger.debug(cmd)
            output = subprocess.getoutput(cmd)
            import os

            if not os.path.exists("out.mol2"):
                msg = "antechamber failed to produce output mol2 file\n"
                msg += f"command: {cmd}\n"
                msg += "output:\n"
                msg += 8 * "----------" + "\n"
                msg += output
                msg += 8 * "----------" + "\n"
                msg += "input:\n"
                msg += 8 * "----------" + "\n"
                msg += read_file_contents(local_input_filename)
                msg += 8 * "----------" + "\n"
                # TODO: Run antechamber again with acdoctor mode on (-dr yes) to get more debug info, if supported
                os.chdir(cwd)
                raise Exception(msg)
            _logger.debug(output)

            # Run parmchk.
            shutil.copy(self.gaff_dat_filename, "gaff.dat")
            cmd = f"parmchk2 -i out.mol2 -f mol2 -p gaff.dat -o out.frcmod -s {self._gaff_major_version} -a Y"

            _logger.debug(cmd)
            output = subprocess.getoutput(cmd)
            if not os.path.exists("out.frcmod"):
                msg = "parmchk2 failed to produce output frcmod file\n"
                msg += f"command: {cmd}\n"
                msg += "output:\n"
                msg += 8 * "----------" + "\n"
                msg += output
                msg += 8 * "----------" + "\n"
                msg += "input mol2:\n"
                msg += 8 * "----------" + "\n"
                msg += read_file_contents("out.mol2")
                msg += 8 * "----------" + "\n"
                os.chdir(cwd)
                raise Exception(msg)
            _logger.debug(output)
            self._check_for_errors(output)

            # Copy back
            shutil.copy("out.mol2", gaff_mol2_filename)
            shutil.copy("out.frcmod", frcmod_filename)

            os.chdir(cwd)

        return gaff_mol2_filename, frcmod_filename

    def _read_gaff_atom_types_from_mol2(self, gaff_mol2_filename, molecule):
        """
        Read the GAFF atom types specified in an antechamber-generated mol2 file into atom.gaff_type in
        the specified Molecule

        Parameters
        ----------
        gaff_mol2_filename : str
            The antechamber-generated mol2 file containing GAFF/GAFF2 atom types
        molecule : Molecule
            The Molecule to receive atom types
            The Atom objects within the molecule will have a ``gaff_type`` field added containing the GAFF
            atom type as a string
        """

        # Read the resulting GAFF mol2 file atom types
        #       1 C1           1.8850    -1.0360    -0.1120 ca         1 MOL       0.000000
        # 012345678901234567890123456789012345678901234567890123456789012345678901234567890
        # 0         1         2         3         4         5         6         7         8
        with open(gaff_mol2_filename) as infile:
            line = infile.readline()
            # Seek to ATOM block
            while line:
                if line.strip() == "@<TRIPOS>ATOM":
                    break
                line = infile.readline()
            # Read GAFF atom types
            for index, atom in enumerate(molecule.atoms):
                line = infile.readline()
                atom.gaff_type = line[50:58].strip()

    def _check_for_errors(self, outputtext, other_errors=None, ignore_errors=None):
        """
        Check AMBER package output for the string 'ERROR' (upper or lowercase) and (optionally) specified other
        strings and raise an exception if it is found (to avoid silent failures which might be noted to log but
        otherwise ignored).

        Parameters
        ----------
        outputtext : str
            String listing output text from an (AMBER) command which should be checked for errors.
        other_errors : list(str), default None
            If specified, provide strings for other errors which will be checked for, such as
            "improper number of arguments", etc.
        ignore_errors: list(str), default None
            If specified, AMBER output lines containing errors but also containing any of the specified strings
            will be ignored (because, for example, AMBER issues an "ERROR" for non-integer charges in some cases
            when only a warning is needed).

        Notes
        -----
        If error(s) are found, raise a RuntimeError and attept to print the appropriate errors from the processed text.
        """

        lines = outputtext.split("\n")
        error_lines = []
        for line in lines:
            if "ERROR" in line.upper():
                error_lines.append(line)
            if other_errors is not None:
                for err in other_errors:
                    if err.upper() in line.upper():
                        error_lines.append(line)

        if ignore_errors is not None and len(error_lines) > 0:
            new_error_lines = []
            for ign in ignore_errors:
                ignore = False
                for err in error_lines:
                    if ign in err:
                        ignore = True
                if not ignore:
                    new_error_lines.append(err)
            error_lines = new_error_lines

        if len(error_lines) > 0:
            _logger.warning("Unexpected errors encountered running AMBER tool. Offending output:")
            for line in error_lines:
                _logger.warning(line)
            raise (RuntimeError("Error encountered running AMBER tool. Exiting."))


################################################################################
# Mixin for OpenMM System-generating template generators (used by SMIRNOFF and espaloma)
################################################################################


class OpenMMSystemMixin:
    """Mixin for force field template generators that produce OpenMM System objects"""

    def clear_system_cache(self):
        """Initialize the OpenMM System cache"""
        self._system_cache = dict()

    def cache_system(self, smiles, system):
        """
        Transiently cache a copy of the OpenMM System object

        Parameters
        ----------
        smiles : str
            The SMILES corresponding to the System object
        system : openmm.System
            The OpenMM System to cache
        """

        self._system_cache[smiles] = system

    def get_openmm_system(self, molecule):
        """
        Retrieve the OpenMM System object generated for a particular molecule for testing/validation.

        Parameters
        ----------
        molecule : openff.toolkit.Molecule
            The Molecule object

        Returns
        -------
        system : openmm.System or None
            If the Molecule object has already been parameterized by this instance, this molecule is returned.
            Otherwise, None is returned.
        """

        smiles = molecule.to_smiles()
        if smiles in self._system_cache:
            return self._system_cache[smiles]
        else:
            return None

    def convert_system_to_ffxml(self, molecule, system, improper_atom_ordering="smirnoff"):
        """
        Convert OpenMM System object to molecule-specific OpenMM ffxml

        Parameters
        ----------
        molecule : openff.toolkit.Molecule
            The Molecule to be converted
        system : openmm.System
            The System corresponding to molecule
        improper_atom_ordering : str, optional, default='smirnoff'
            OpenMM openmm.app.ForceField improper atom ordering scheme to use

        Returns
        -------
        ffxml_contents : str
            The OpenMM ffxml contents for the given molecule.
        """

        from openmm import CMMotionRemover
        from lxml import etree

        # Remove CMMotionRemover if present
        # See https://github.com/openmm/openmmforcefields/issues/365
        # and https://github.com/openmm/openmmforcefields/pull/367
        for f_idx in reversed(range(system.getNumForces())):
            force = system.getForce(f_idx)
            if isinstance(force, CMMotionRemover):
                system.removeForce(f_idx)

        # Generate OpenMM ffxml definition for this molecule
        root = etree.Element("ForceField")

        def as_attrib(quantity):
            """Format openff.units.Quantity or openmm.unit.Quantity as XML attribute."""
            import openff.units

            if isinstance(quantity, str):
                return quantity
            elif isinstance(quantity, (float, int)):
                return str(quantity)
            elif isinstance(quantity, openff.units.Quantity):
                # TODO: Match behavior of Quantity.value_in_unit_system
                return str(quantity.m)
            else:
                from openmm.unit import Quantity as OpenMMQuantity

                if isinstance(quantity, OpenMMQuantity):
                    from openmm import unit

                    return str(quantity.value_in_unit_system(unit.md_unit_system))
                else:
                    raise ValueError(f"Found unexpected type {type(quantity)}.")

        # Append unique type names to atoms
        smiles = molecule.to_smiles()
        for index, atom in enumerate(molecule.atoms):
            setattr(atom, "typename", f"{smiles}${atom.name}#{index}")

        # Generate atom types
        atom_types = etree.SubElement(root, "AtomTypes")
        for atom_index, atom in enumerate(molecule.atoms):
            # Create a new atom type for each atom in the molecule
            element_symbol = atom.symbol
            atom_type = etree.SubElement(
                atom_types, "Type", name=atom.typename, element=element_symbol, mass=as_attrib(atom.mass)
            )
            atom_type.set("class", atom.typename)  # 'class' is a reserved Python keyword, so use alternative API

        supported_forces = {
            "NonbondedForce",
            "HarmonicAngleForce",
            "HarmonicBondForce",
            "PeriodicTorsionForce",
        }

        # Compile forces into a dict
        forces = dict()
        for force in system.getForces():
            force_name = force.__class__.__name__

            if force_name in forces:
                raise ForceException(f"Two instances of force {force_name} appear in System")
            if force_name not in supported_forces:
                raise ForceException(f"Custom forces not supported. Found force of type {force_name}.")

            forces[force_name] = force

        def classes(atom_indices):
            """Build a dict of 'class#=typename' for use in creating XML tags for forces.

            Parameters
            ----------
            atom_indices : list of int
                Particle indices for molecule.atoms

            Returns
            -------
            classmap : dict of str : str
                Dict of format { 'class1' : typename1, ... }
            """
            return {
                f"class{class_index + 1}": molecule.atoms[atom_index].typename
                for class_index, atom_index in enumerate(atom_indices)
            }

        # Lennard-Jones
        # In case subclasses specifically set the 1-4 scaling factors, use those.
        nonbonded_types = etree.SubElement(
            root,
            "NonbondedForce",
            coulomb14scale=getattr(self, "_coulomb14scale", "0.833333"),
            lj14scale=getattr(self, "_lj14scale", "0.5"),
        )
        etree.SubElement(nonbonded_types, "UseAttributeFromResidue", name="charge")
        for atom_index in range(forces["NonbondedForce"].getNumParticles()):
            charge, sigma, epsilon = forces["NonbondedForce"].getParticleParameters(atom_index)
            nonbonded_type = etree.SubElement(
                nonbonded_types,
                "Atom",
                sigma=as_attrib(sigma),
                epsilon=as_attrib(epsilon),
            )
            nonbonded_type.set(
                "class", molecule.atoms[atom_index].typename
            )  # 'class' is a reserved Python keyword, so use alternative API

        # Bonds
        bond_types = etree.SubElement(root, "HarmonicBondForce")
        atom_indices = [-1] * 2
        for bond_index in range(forces["HarmonicBondForce"].getNumBonds()):
            atom_indices[0], atom_indices[1], length, k = forces["HarmonicBondForce"].getBondParameters(bond_index)

            etree.SubElement(
                bond_types,
                "Bond",
                **classes(atom_indices),
                length=as_attrib(length),
                k=as_attrib(k),
            )

        # Angles
        angle_types = etree.SubElement(root, "HarmonicAngleForce")
        atom_indices = [-1] * 3
        for angle_index in range(forces["HarmonicAngleForce"].getNumAngles()):
            atom_indices[0], atom_indices[1], atom_indices[2], angle, k = forces[
                "HarmonicAngleForce"
            ].getAngleParameters(angle_index)

            etree.SubElement(
                angle_types,
                "Angle",
                **classes(atom_indices),
                angle=as_attrib(angle),
                k=as_attrib(k),
            )

        # Torsions
        def torsion_tag(atom_indices):
            """Return 'Proper' or 'Improper' depending on torsion type"""
            atoms = [molecule.atoms[atom_index] for atom_index in atom_indices]
            # TODO: Check to make sure all atoms are in fact atoms and not virtual sites
            if atoms[0].is_bonded_to(atoms[1]) and atoms[1].is_bonded_to(atoms[2]) and atoms[2].is_bonded_to(atoms[3]):
                return "Proper"
            else:
                return "Improper"

        # Collect torsions
        torsions = dict()
        for torsion_index in range(forces["PeriodicTorsionForce"].getNumTorsions()):
            atom_indices = [-1] * 4
            (
                atom_indices[0],
                atom_indices[1],
                atom_indices[2],
                atom_indices[3],
                periodicity,
                phase,
                k,
            ) = forces["PeriodicTorsionForce"].getTorsionParameters(torsion_index)
            atom_indices = tuple(atom_indices)
            if atom_indices in torsions.keys():
                torsions[atom_indices].append((periodicity, phase, k))
            else:
                torsions[atom_indices] = [(periodicity, phase, k)]

        # Create torsion definitions
        torsion_types = etree.SubElement(root, "PeriodicTorsionForce", ordering=improper_atom_ordering)
        for atom_indices in torsions.keys():
            params = dict()  # build parameter dictionary
            nterms = len(torsions[atom_indices])
            for term in range(nterms):
                periodicity, phase, k = torsions[atom_indices][term]
                params[f"periodicity{term + 1}"] = as_attrib(periodicity)
                params[f"phase{term + 1}"] = as_attrib(phase)
                params[f"k{term + 1}"] = as_attrib(k)

            etree.SubElement(
                torsion_types,
                torsion_tag(atom_indices),
                **classes(atom_indices),
                **params,
            )

        # TODO: Handle virtual sites
        virtual_sites = [
            atom_index for atom_index in range(system.getNumParticles()) if system.isVirtualSite(atom_index)
        ]
        if len(virtual_sites) > 0:
            raise Exception("Virtual sites are not yet supported")

        # Create residue definitions
        # TODO: Handle non-Atom atoms too (virtual sites)
        residues = etree.SubElement(root, "Residues")
        residue = etree.SubElement(residues, "Residue", name=smiles)
        for atom_index, atom in enumerate(molecule.atoms):
            charge, sigma, epsilon = forces["NonbondedForce"].getParticleParameters(atom_index)
            etree.SubElement(
                residue,
                "Atom",
                name=atom.name,
                type=atom.typename,
                charge=as_attrib(charge),
            )
        for bond in molecule.bonds:
            etree.SubElement(residue, "Bond", atomName1=bond.atom1.name, atomName2=bond.atom2.name)

        # Render XML into string
        ffxml_contents = etree.tostring(root, pretty_print=True, encoding="unicode")

        # _logger.debug(f'{ffxml_contents}') # DEBUG

        return ffxml_contents


################################################################################
# Open Force Field Initiative SMIRNOFF specific OpenMM ForceField template generation utilities
################################################################################


class SMIRNOFFTemplateGenerator(SmallMoleculeTemplateGenerator, OpenMMSystemMixin):
    """
    OpenMM ForceField residue template generator for Open Force Field Initiative SMIRNOFF
    force fields using pre-cached OpenFF toolkit molecules.

    Open Force Field Initiative: http://openforcefield.org
    SMIRNOFF force field specification: https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html

    Examples
    --------

    Create a template generator for a single Molecule using the latest Open Force Field Initiative
    small molecule force field and register it with ForceField:

    >>> # Define a Molecule using the OpenFF Molecule object
    >>> from openff.toolkit import Molecule
    >>> molecule = Molecule.from_smiles('c1ccccc1')
    >>> # Create the SMIRNOFF template generator
    >>> from openmmforcefields.generators import SMIRNOFFTemplateGenerator
    >>> template_generator = SMIRNOFFTemplateGenerator(molecules=molecule)
    >>> # Create an OpenMM ForceField
    >>> from openmm.app import ForceField
    >>> amber_forcefields = ['amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml']
    >>> forcefield = ForceField(*amber_forcefields)
    >>> # Register the template generator
    >>> forcefield.registerTemplateGenerator(template_generator.generator)

    Create a template generator for a specific pre-installed SMIRNOFF version ('openff-2.0.0')
    and register multiple molecules:

    >>> molecule1 = Molecule.from_smiles('c1ccccc1')
    >>> molecule2 = Molecule.from_smiles('CCO')
    >>> template_generator = SMIRNOFFTemplateGenerator(molecules=[molecule1, molecule2], forcefield='openff-2.0.0')

    Alternatively, you can specify a local .offxml file in the SMIRNOFF specification:

    >>> template_generator = SMIRNOFFTemplateGenerator(molecules=[molecule1, molecule2], forcefield='mysmirnoff.offxml')  # doctest: +SKIP

    You can also add some Molecules later on after the generator has been registered:

    >>> template_generator.add_molecules(molecule)
    >>> template_generator.add_molecules([molecule1, molecule2])

    You can optionally create or use a tiny database cache of pre-parameterized molecules:

    >>> template_generator = SMIRNOFFTemplateGenerator(cache='smirnoff-molecules.json')

    Newly parameterized molecules will be written to the cache, saving time next time!
    """  # noqa

    def __init__(self, molecules=None, cache=None, forcefield=None, template_generator_kwargs=None):
        """
        Create a SMIRNOFFTemplateGenerator with some OpenFF toolkit molecules

        Requies the OpenFF toolkit: http://openforcefield.org

        Parameters
        ----------
        molecules : openff.toolkit .Molecule or list, optional, default=None
            Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can be used to construct a Molecule.
            Can also be a list of Molecule objects or objects that can be used to construct a Molecule.
            If specified, these molecules will be recognized and parameterized with SMIRNOFF as needed.
            The parameters will be cached in case they are encountered again the future.
        cache : str, optional, default=None
            Filename for global caching of parameters.
            If specified, parameterized molecules will be stored in a TinyDB instance as a JSON file.
        forcefield : str, optional, default=None
            Name of installed SMIRNOFF force field (without .offxml) or local .offxml filename (with extension).
            If not specified, the latest Open Force Field Initiative release is used.
        template_generator_kwargs : dict, optional, default=None
            Additional parameters for the template generator (ignored by SMIRNOFFTemplateGenerator).

        Examples
        --------

        Create a SMIRNOFF template generator for a single molecule (benzene, created from SMILES) and register it with ForceField:

        >>> from openff.toolkit import Molecule
        >>> molecule = Molecule.from_smiles('c1ccccc1')
        >>> from openmmforcefields.generators import SMIRNOFFTemplateGenerator
        >>> smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule)
        >>> from openmm.app import ForceField
        >>> amber_forcefields = ['amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml']
        >>> forcefield = ForceField(*amber_forcefields)

        The latest Open Force Field Initiative release is used if none is specified.

        >>> smirnoff.forcefield
        'openff-2.1.0'

        You can check which SMIRNOFF force field filename is in use with

        >>> smirnoff.smirnoff_filename  # doctest:+ELLIPSIS
        '/.../openff-2.1.0.offxml'

        Create a template generator for a specific SMIRNOFF force field for multiple
        molecules read from an SDF file:

        >>> molecules = Molecule.from_file('molecules.sdf')  # doctest: +SKIP
        >>> smirnoff = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield='openff-2.1.0')  # doctest: +SKIP

        You can also add molecules later on after the generator has been registered:

        >>> smirnoff.add_molecules(molecules)  # doctest: +SKIP

        To check which SMIRNOFF versions are supported, check the `INSTALLED_FORCEFIELDS` attribute:

        >>> print(SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS)  # doctest: +SKIP
        ['openff-1.0.1', 'openff-1.1.1', 'openff-1.0.0-RC1', 'openff-1.2.0', 'openff-1.1.0', 'openff-1.0.0', 'openff-1.0.0-RC2', 'smirnoff99Frosst-1.0.2', 'smirnoff99Frosst-1.0.0', 'smirnoff99Frosst-1.1.0', 'smirnoff99Frosst-1.0.4', 'smirnoff99Frosst-1.0.8', 'smirnoff99Frosst-1.0.6', 'smirnoff99Frosst-1.0.3', 'smirnoff99Frosst-1.0.1', 'smirnoff99Frosst-1.0.5', 'smirnoff99Frosst-1.0.9', 'smirnoff99Frosst-1.0.7']

        You can optionally create or use a cache of pre-parameterized molecules:

        >>> smirnoff = SMIRNOFFTemplateGenerator(cache='smirnoff.json', forcefield='openff-2.1.0')  # doctest: +SKIP

        Newly parameterized molecules will be written to the cache, saving time next time!
        """  # noqa

        self._lj14scale = None
        self._coulomb14scale = None

        # Initialize molecules and cache
        super().__init__(molecules=molecules, cache=cache)

        if forcefield is None:
            # Use latest supported Open Force Field Initiative release if none is specified
            forcefield = "openff-2.1.0"
            # TODO: After toolkit provides date-ranked force fields,
            # use latest dated version if we can sort by date, such as self.INSTALLED_FORCEFIELDS[-1]
        self._forcefield = forcefield

        # Track parameters by provided SMIRNOFF name
        # TODO: Can we instead use the force field hash, or some other unique identifier?
        # TODO: Use file hash instead of name?
        import os

        self._database_table_name = os.path.basename(forcefield)

        # Create ForceField object
        import openff.toolkit.typing.engines.smirnoff

        # check for an installed force field
        available_force_fields = openff.toolkit.typing.engines.smirnoff.get_available_force_fields()
        if (filename := forcefield + ".offxml") in available_force_fields or (
            filename := forcefield
        ) in available_force_fields:
            self._smirnoff_forcefield = openff.toolkit.typing.engines.smirnoff.ForceField(filename)

        # just try parsing the input and let openff handle the error
        else:
            try:
                self._smirnoff_forcefield = openff.toolkit.typing.engines.smirnoff.ForceField(forcefield)
            except Exception as e:
                _logger.error(e)
                raise ValueError(
                    f"Can't find specified SMIRNOFF force field ({forcefield}) in install paths "
                    "or parse the input as a string."
                ) from e

        self._coulomb14scale = str(self._smirnoff_forcefield.get_parameter_handler("Electrostatics").scale14)
        self._lj14scale = str(self._smirnoff_forcefield.get_parameter_handler("vdW").scale14)

        # Delete constraints, if present
        if "Constraints" in self._smirnoff_forcefield._parameter_handlers:
            del self._smirnoff_forcefield._parameter_handlers["Constraints"]

        # Find SMIRNOFF filename
        smirnoff_filename = self._search_paths(filename)
        self._smirnoff_filename = smirnoff_filename

        # Cache a copy of the OpenMM System generated for each molecule for testing purposes
        self.clear_system_cache()

    @classproperty
    def INSTALLED_FORCEFIELDS(cls):
        """
        Return a list of the offxml files shipped with the openff-forcefields package.

        Returns
        -------
        file_names : str
           The file names of available force fields

        .. todo ::

           Replace this with an API call once this issue is addressed:
           https://github.com/openforcefield/openff-toolkit/issues/477
        """

        from openff.toolkit.typing.engines.smirnoff import get_available_force_fields

        file_names = list()
        for filename in get_available_force_fields(full_paths=False):
            root, ext = os.path.splitext(filename)
            # Only add variants without '_unconstrained'
            if "_unconstrained" in root:
                continue
            # The OpenFF Toolkit ships two versions of its ff14SB port, one with SMIRNOFF-style
            # impropers and one with Amber-style impropers. The latter requires a special handler
            # (`AmberImproperTorsionHandler`) that is not shipped with the toolkit. See
            # https://github.com/openforcefield/amber-ff-porting/tree/0.0.3
            if root.startswith("ff14sb") and "off_impropers" not in root:
                continue
            file_names.append(root)

        return file_names

    def _search_paths(self, filename):
        """
        Search registered OpenFF plugin directories

        Parameters
        ----------
        filename : str
            The filename to find the full path for

        Returns
        -------
        fullpath : str
            Full path to identified file, or None if no file found

        .. todo ::

           Replace this with an API call once this issue is addressed:
           https://github.com/openforcefield/openff-toolkit/issues/477
        """

        # TODO: Replace this method once there is a public API in the OpenFF toolkit for doing this

        from openff.toolkit.typing.engines.smirnoff.forcefield import (
            _get_installed_offxml_dir_paths,
        )

        # Check whether this could be a file path
        if isinstance(filename, str):
            # Try first the simple path.
            searched_dirs_paths = [""]
            # Then try a relative file path w.r.t. an installed directory.
            searched_dirs_paths.extend(_get_installed_offxml_dir_paths())

            # Determine the actual path of the file.
            # TODO: What is desired toolkit behavior if two files with the desired name are available?
            for dir_path in searched_dirs_paths:
                file_path = os.path.join(dir_path, filename)
                if os.path.isfile(file_path):
                    return file_path
        # No file found
        return None

    @property
    def smirnoff_filename(self):
        """Full path to the SMIRNOFF force field file"""
        return self._smirnoff_filename

    def generate_residue_template(self, molecule):
        """
        Generate a residue template and additional parameters for the specified Molecule.

        Parameters
        ----------
        molecule : openff.toolkit.Molecule
            A molecule to generate a residue template for. The molecule will be
            parameterized with SMIRNOFF as needed. The parameters will be
            cached in case it is encountered again the future.

        Returns
        -------
        ffxml_contents : str
            Contents of ForceField `ffxml` file containing additional parameters and residue template.

        Notes
        -----

        * The residue template will be named after the SMILES of the molecule.
        * This method preserves stereochemistry during AM1-BCC charge parameterization.
        * Atom names in molecules will be assigned Tripos atom names if any are blank or not unique.
        """

        # Use the canonical isomeric SMILES to uniquely name the template
        smiles = molecule.to_smiles()
        _logger.info(f"Generating a residue template for {smiles} using {self._forcefield}")

        # Generate unique atom names
        self._generate_unique_atom_names(molecule)

        has_user_charges = self._molecule_has_user_charges(molecule)

        # Determine which molecules (if any) contain user-specified partial charges that should be used
        charge_from_molecules = list()
        if has_user_charges:
            charge_from_molecules = [molecule]
            _logger.debug("Using user-provided charges because partial charges are nonzero...")

        # Parameterize molecule
        _logger.debug("Generating parameters...")
        system = self._smirnoff_forcefield.create_openmm_system(
            molecule.to_topology(),
            charge_from_molecules=charge_from_molecules,
            # "allow_nonintegral_charges" is a misnomer since the actual check
            # that OpenFF does will raise an error even if the user charges sum
            # to an integer but do not match the formal charge.  Since we have
            # already warned about this if it is the case, allow it.
            allow_nonintegral_charges=has_user_charges,
        )

        self.cache_system(smiles, system)

        # Convert to ffxml
        ffxml_contents = self.convert_system_to_ffxml(molecule, system)
        return ffxml_contents


################################################################################
# Espaloma template generation utilities
################################################################################


class EspalomaTemplateGenerator(SmallMoleculeTemplateGenerator, OpenMMSystemMixin):
    """
    OpenMM ForceField residue template generator for espaloma force fields using pre-cached OpenFF toolkit molecules.

    Espaloma uses a graph net approach to chemical perception to assign parameters and charges.

    * Espaloma docs and papers: https://docs.espaloma.org/
    * Espaloma code and models: https://github.com/choderalab/espaloma
    * Open Force Field Initiative: http://openforcefield.org

    .. warning :: This API is experimental and subject to change.

    Examples
    --------

    Create a template generator for a single Molecule using the latest Open Force Field Initiative
    small molecule force field and register it with ForceField:

    >>> # Define a Molecule using the OpenFF Molecule object
    >>> from openff.toolkit import Molecule
    >>> molecule = Molecule.from_smiles('c1ccccc1')
    >>> # Create the Espaloma template generator
    >>> from openmmforcefields.generators import EspalomaTemplateGenerator
    >>> template_generator = EspalomaTemplateGenerator(molecules=molecule)
    >>> # Create an OpenMM ForceField
    >>> from openmm.app import ForceField
    >>> amber_forcefields = ['amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml']
    >>> forcefield = ForceField(*amber_forcefields)
    >>> # Register the template generator
    >>> forcefield.registerTemplateGenerator(template_generator.generator)

    Create a template generator for a specific Espaloma release ('espaloma-0.3.2')
    and register multiple molecules:

    >>> molecule1 = Molecule.from_smiles('c1ccccc1')
    >>> molecule2 = Molecule.from_smiles('CCO')
    >>> template_generator = EspalomaTemplateGenerator(
    ...     molecules=[molecule1, molecule2],
    ...     forcefield='espaloma-0.3.2',
    ... )

    Alternatively, you can specify a local .pt parameter file for Espaloma:

    >>> template_generator = EspalomaTemplateGenerator(
    ...     molecules=[molecule1, molecule2],
    ...     forcefield='espaloma-0.3.2.pt',
    ... )

    You can also add some Molecules later on after the generator has been registered:

    >>> template_generator.add_molecules(molecule)
    >>> template_generator.add_molecules([molecule1, molecule2])

    You can optionally create or use a tiny database cache of pre-parameterized molecules:

    >>> template_generator = EspalomaTemplateGenerator(cache='espaloma-molecules.json')

    Newly parameterized molecules will be written to the cache, saving time next time!
    """

    CHARGE_METHODS = ("nn", "am1-bcc", "gasteiger", "from-molecule")

    def __init__(
        self,
        molecules=None,
        cache=None,
        forcefield=None,
        model_cache_path=None,
        template_generator_kwargs=None,
    ):
        """
        Create an EspalomaTemplateGenerator with some OpenFF toolkit molecules

        Requies the OpenFF toolkit: http://openforcefield.org
        and espaloma: http://github.com/choderalab/espaloma

        Parameters
        ----------
        molecules : openff.toolkit.Molecule or list, optional, default=None
            Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can
            be used to construct a Molecule.
            Can also be a list of Molecule objects or objects that can be used to construct a Molecule.
            If specified, these molecules will be recognized and parameterized with espaloma as needed.
            The parameters will be cached in case they are encountered again the future.
        cache : str, optional, default=None
            Filename for global caching of parameters.
            If specified, parameterized molecules will be stored in a TinyDB instance as a JSON file.
        forcefield : str, optional, default=None
            Name of installed Espaloma force field version (e.g. 'espaloma-0.3.2') to retrieve remotely,
            a local Espaloma .pt parmaeters filename (with extension),
            or a URL to an online espaloma force field.
        model_cache_path : str, optional, default=None
            If specified, use this directory to cache espaloma models
            default: ~/.espaloma/
        template_generator_kwargs : dict, optional, default=None
            An optional dictionary containing additional optional arguments:
            {"reference_forcefield": str, Openff force field supported by https://github.com/openforcefield/openff-forcefields
            without .offxml extension}
            {"charge_method": str, Charge method supported by espaloma ['nn', 'am1-bcc', 'gasteiger', 'from-molecule']}

            Default behavior is to use ``openff_unconstrained-2.0.0`` for ``reference_forcefield`` and
            `nn` for `charge_method`.
            User defined charges can be assigned by setting the ``charge_method`` to ``from_molecule``
            if charges are assigned to openff.toolkit.Molecule.

        Examples
        --------

        Create an Espaloma template generator for a single molecule (benzene, created from SMILES) and register it with ForceField:

        >>> from openff.toolkit import Molecule
        >>> molecule = Molecule.from_smiles('c1ccccc1')
        >>> from openmmforcefields.generators import EspalomaTemplateGenerator
        >>> espaloma_generator = EspalomaTemplateGenerator(molecules=molecule)
        >>> from openmm.app import ForceField
        >>> amber_forcefields = ['amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml']
        >>> forcefield = ForceField(*amber_forcefields)

        >>> espaloma_generator.forcefield
        'espaloma-0.3.2'

        You can check which espaloma parameter set force field filename is in use with

        >>> espaloma_generator.espaloma_filename
        '/.../espaloma-0.3.2.pt'

        Create a template generator for a specific SMIRNOFF force field for multiple
        molecules read from an SDF file or list of SMILES strings:

        >>> molecules = Molecule.from_file('molecules.sdf')  # doctest: +SKIP
        >>> molecules = [Molecule.from_smiles(smiles) for smiles in ["CCO", "c1ccccc1"]]
        >>> espaloma_generator = EspalomaTemplateGenerator(molecules=molecules, forcefield='espaloma-0.3.2')

        You can also add molecules later on after the generator has been registered:

        >>> espaloma_generator.add_molecules(molecules[-1])

        You can optionally create or use a cache of pre-parameterized molecules:

        >>> espaloma_generator = EspalomaTemplateGenerator(cache='smirnoff.json', forcefield='espaloma-0.3.2')

        Newly parameterized molecules will be written to the cache, saving time next time!

        You can also pass a template_generator_kwargs to specify the reference_forcefield and/or charge_method in EspalomaTemplateGenerator:

        >>> template_generator_kwargs = {"reference_forcefield": "openff_unconstrained-2.0.0", "charge_method": "nn"}
        >>> espaloma_generator = EspalomaTemplateGenerator(cache='smirnoff.json', forcefield='espaloma-0.3.2', template_generator_kwargs=template_generator_kwargs)
        """  # noqa

        # Initialize molecules and cache
        super().__init__(molecules=molecules, cache=cache)

        # Espaloma model cache path
        if model_cache_path is None:
            import os

            self.ESPALOMA_MODEL_CACHE_PATH = f"{os.getenv('HOME')}/.espaloma"
        else:
            self.ESPALOMA_MODEL_CACHE_PATH = model_cache_path

        if forcefield is None:
            # Use latest supported Espaloma force field release if none is specified
            forcefield = "espaloma-0.3.2"
            # TODO: After toolkit provides date-ranked force fields,
            # use latest dated version if we can sort by date, such as self.INSTALLED_FORCEFIELDS[-1]
        self._forcefield = forcefield

        # Check that espaloma model parameters can be found or locally cached
        self.espaloma_model_filepath = self._get_model_filepath(forcefield)

        # Check reference forcefield and charge method
        if template_generator_kwargs is not None:
            self._reference_forcefield = template_generator_kwargs.get(
                "reference_forcefield", "openff_unconstrained-2.0.0"
            )
            self._charge_method = template_generator_kwargs.get("charge_method", "nn")
        else:
            # Consider upgrading to 2.1.0, the recommended small moleucle force field for general use
            self._reference_forcefield = "openff_unconstrained-2.0.0"
            self._charge_method = "from-molecule"

        # Check to make sure dependencies are installed
        try:
            import espaloma  # noqa
        except ImportError:
            raise ValueError("The EspalomaResidueTemplateGenerator requires espaloma to be installed")

        # Check force field can be found

        # Track parameters by provided force field name
        # TODO: Can we instead use the force field hash, or some other unique identifier?
        # TODO: Use file hash instead of name?
        import os

        self._database_table_name = os.path.basename(forcefield)

        # Load torch model
        import torch

        self.espaloma_model = torch.load(
            self.espaloma_model_filepath, map_location=torch.device("cpu"), weights_only=False
        )
        self.espaloma_model.eval()

        # Cache a copy of the OpenMM System generated for each molecule for testing purposes
        self.clear_system_cache()

    @classproperty
    def INSTALLED_FORCEFIELDS(self):
        """Return list of available force field versions."""
        # TODO: Does this belong here? Is there a better way to do this?
        # TODO: Update this
        # TODO: Can we list force fields installed locally?
        # TODO: Maybe we can check ~/.espaloma and ESPALOMA_PATH?
        return ["espaloma-0.3.2"]

    def _get_model_filepath(self, forcefield):
        """
        Retrieve local file path to cached espaloma model parameters, or retrieve remote model if needed.

        Parameters
        ----------
        forcefield : str
            Version to locate in local cache (or retrieve if needed)

        Returns
        -------
        cached_filename : str
            Path to local cache of espaloma .pt model parameters
        """

        import os

        if os.path.exists(forcefield):
            # A specific file path has been specified
            _logger.info(f"Using espaloma model found at {forcefield}")
            return forcefield
        # TODO: This isn't quite right---we should be checking this in the previous branch?
        elif os.path.exists(os.path.join(self.ESPALOMA_MODEL_CACHE_PATH, forcefield)):
            # A specific file path has been specified
            filepath = os.path.join(self.ESPALOMA_MODEL_CACHE_PATH, forcefield)
            _logger.info(f"Using espaloma model found at {filepath}")
            return filepath
        else:
            import validators

            if validators.url(forcefield):
                # URL has been provided
                url = forcefield
                filename = os.path.basename(url)  # local filename for caching
            else:
                # Identify version number
                import re

                m = re.match(r"espaloma-(\d+\.\d+\.\d+)", forcefield)
                if m is None:
                    raise ValueError(
                        f'Espaloma model must be filepath or formatted like "espaloma-0.3.2" (found: "{forcefield}")'
                    )
                version = m.group(1)
                # Construct URL
                url = f"https://github.com/choderalab/espaloma/releases/download/{version}/espaloma-{version}.pt"
                filename = f"espaloma-{version}.pt"  # local filename for caching

            # Check cache
            cached_filename = os.path.join(self.ESPALOMA_MODEL_CACHE_PATH, filename)
            if os.path.exists(cached_filename):
                _logger.info(f"Using espaloma model cached at {cached_filename}")
                return cached_filename
            else:
                # Create the cache directory
                os.makedirs(self.ESPALOMA_MODEL_CACHE_PATH, exist_ok=True)

                # Attempt to retrieve from URL
                _logger.info(f"Attempting to retrieve espaloma model from {url}")
                import tempfile
                import urllib
                import urllib.error
                import urllib.request

                with tempfile.TemporaryDirectory(dir=self.ESPALOMA_MODEL_CACHE_PATH) as temp_dir:
                    temp_filename = os.path.join(temp_dir, filename)
                    try:
                        urllib.request.urlretrieve(url, filename=temp_filename)
                    except urllib.error.URLError:
                        raise ValueError(f"No espaloma model found at expected URL: {url}")
                    except urllib.error.HTTPError as e:
                        raise ValueError(f"An error occurred while retrieving espaloma model from {url} : {e}")
                    os.replace(temp_filename, cached_filename)
                return cached_filename

    @property
    def espaloma_filename(self):
        """Full path to the SMIRNOFF force field file"""
        return self.espaloma_model_filepath

    def generate_residue_template(self, molecule):
        """
        Generate a residue template and additional parameters for the specified Molecule.

        Parameters
        ----------
        molecule : openff.toolkit.Molecule
            A molecule to generate a residue template for. The molecule will be
            parameterized with espaloma as needed. The parameters will be
            cached in case it is encountered again the future.

        Returns
        -------
        ffxml_contents : str
            Contents of ForceField `ffxml` file containing additional parameters and residue template.

        Notes
        -----

        * The residue template will be named after the SMILES of the molecule.
        """

        from openff.units import unit

        # Use the canonical isomeric SMILES to uniquely name the template
        smiles = molecule.to_smiles()
        _logger.info(f"Generating a residue template for {smiles} using {self._forcefield}")

        # Generate unique atom names
        self._generate_unique_atom_names(molecule)

        # Parameterize molecule
        _logger.debug("Generating espaloma parameters...")

        # create an Espaloma Graph object to represent the molecule of interest
        import espaloma as esp

        molecule_graph = esp.Graph(molecule)

        # Regenerate SMIRNOFF impropers
        from espaloma.graphs.utils.regenerate_impropers import regenerate_impropers

        regenerate_impropers(molecule_graph)

        # Book keep partial charges if molecule has user charges
        # NOTE: Charges will be overwritten when the Espaloma Graph object is loaded into an espaloma model
        if self._molecule_has_user_charges(molecule):
            user_charges = molecule.partial_charges.m_as(unit.elementary_charge)
        else:
            user_charges = None

        # Assign parameters
        self.espaloma_model(molecule_graph.heterograph)

        # Create an OpenMM System
        # Update partial charges if charge_method is "from_molecule"
        if self._charge_method == "from-molecule":
            if user_charges is not None:
                import numpy as np
                import torch

                # Handle ValueError:
                # "ValueError: given numpy array has byte order different from the native byte order.
                # Conversion between byte orders is currently not supported."
                molecule_graph.nodes["n1"].data["q"] = (
                    torch.from_numpy(user_charges.astype(np.float32)).unsqueeze(-1).float()
                )
            else:
                # No charges were found in molecule -- defaulting to nn charge method
                warnings.warn("No charges found in molecule. Defaulting to 'nn' charge method.")
                self._charge_method = "nn"

        system = esp.graphs.deploy.openmm_system_from_graph(
            molecule_graph,
            charge_method=self._charge_method,
            forcefield=self._reference_forcefield,
        )
        _logger.info(
            f"Generating a system with charge method {self._charge_method} and "
            f"{self._reference_forcefield} to assign nonbonded parameters"
        )
        self.cache_system(smiles, system)

        # Convert to ffxml
        ffxml_contents = self.convert_system_to_ffxml(molecule, system)
        return ffxml_contents
