"""
Residue template generator for the AMBER GAFF1/2 small molecule force fields.

"""

################################################################################
# IMPORTS
################################################################################

################################################################################
# LOGGER
################################################################################

import os
import logging
import contextlib
_logger = logging.getLogger("openmmforcefields.generators.template_generators")

################################################################################
# Small molecule OpenMM ForceField template generation utilities
################################################################################

class SmallMoleculeTemplateGenerator(object):
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
        Create a tempalte generator with some openforcefield toolkit molecules

        Requies the openforcefield toolkit: http://openforcefield.org

        Parameters
        ----------
        molecules : openforcefield.topology.Molecule or list, optional, default=None
            Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can be used to construct a Molecule.
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
        self._smiles_added_to_db = set() # set of SMILES added to the database this session
        self._database_table_name = None # this must be set by subclasses for cache to function

        # Name of the force field
        self._forcefield = None # this must be set by subclasses

        # File to write ffxml to if requested
        self.debug_ffxml_filename = None

    @property
    def forcefield(self):
        """The current force field name in use"""
        return self._forcefield

    @contextlib.contextmanager
    def _open_db(self):
        """Open the cache database.
        """
        from tinydb import TinyDB
        tinydb_kwargs = { 'sort_keys' : True, 'indent' : 4, 'separators' : (',', ': ') } # for pretty-printing
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
        molecules : openforcefield.topology.Molecule or list of Molecules, optional, default=None
            Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can be used to construct a Molecule.
            Can also be a list of Molecule objects or objects that can be used to construct a Molecule.
            If specified, these molecules will be recognized and parameterized with antechamber as needed.
            The parameters will be cached in case they are encountered again the future.

        Examples
        --------

        Add some Molecules later on after the generator has been registered:

        >>> from openforcefield.topology import Molecule
        >>> mol1, mol2, mol3 = [ Molecule.from_smiles(smiles) for smiles in ('c1ccccc1', 'O=Cc1ccc(O)c(OC)c1', 'CN1CCC[C@H]1c2cccnc2') ]
        >>> generator.add_oemols(mol1)
        >>> generator.add_oemols([mol2, mol3])

        """
        # Return if empty
        if not molecules:
            return

        # Ensure molecules is an iterable
        try:
            iterator = iter(molecules)
        except TypeError as te:
            molecules = [ molecules ]

        # Create copies
        # TODO: Do we need to try to construct molecules with other schemes, such as Molecule.from_smiles(), if needed?
        import copy
        molecules = [ copy.deepcopy(molecule) for molecule in molecules ]

        # Cache molecules
        self._molecules.update( { molecule.to_smiles() : molecule for molecule in molecules } )

    @staticmethod
    def _match_residue(residue, molecule_template):
        """Determine whether a residue matches a Molecule template and return a list of corresponding atoms.

        This implementation uses NetworkX for graph isomorphism determination.

        Parameters
        ----------
        residue : simtk.openmm.app.topology.Residue
            The residue to check
        molecule_template : openforcefield.topology.Molecule
            The Molecule template to compare it to

        Returns
        -------
        matches : dict of int : int
            matches[residue_atom_index] is the corresponding Molecule template atom index
            or None if it does not match the template

        .. todo :: Can this be replaced by an isomorphism matching call to the openforcefield toolkit?

        """
        # TODO: Speed this up by rejecting molecules that do not have the same chemical formula

        # TODO: Can this NetworkX implementation be replaced by an isomorphism
        # matching call to the openforcefield toolkit?

        import networkx as nx

        # Build list of external bonds for residue
        number_of_external_bonds = { atom : 0 for atom in residue.atoms() }
        for bond in residue.external_bonds():
            if bond[0] in number_of_external_bonds: number_of_external_bonds[bond[0]] += 1
            if bond[1] in number_of_external_bonds: number_of_external_bonds[bond[1]] += 1

        # Residue graph
        residue_graph = nx.Graph()
        for atom in residue.atoms():
            residue_graph.add_node(atom, element=atom.element.atomic_number, number_of_external_bonds=number_of_external_bonds[atom])
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
        #print(f'residue_graph: nodes {len(list(residue_graph.nodes))} edges {len(list(residue_graph.edges))}')
        #print(f'template_graph: nodes {len(list(template_graph.nodes))} edges {len(list(template_graph.edges))}')

        # Determine graph isomorphism
        from networkx.algorithms import isomorphism
        def node_match(n1, n2):
            """Return True of nodes match, and False if not"""
            return (n1['element']==n2['element']) and (n1['number_of_external_bonds']==n2['number_of_external_bonds'])
        graph_matcher = isomorphism.GraphMatcher(residue_graph, template_graph, node_match=node_match)
        if graph_matcher.is_isomorphic() == False:
            return None

        # Translate to local residue atom indices
        # TODO: This can be simplified because molecule_template uses atom index as key
        atom_index_within_residue = { atom : index for (index, atom) in enumerate(residue.atoms()) }
        atom_index_within_template = { index : index for (index, atom) in enumerate(molecule_template.atoms) }
        matches = { atom_index_within_residue[residue_atom] : atom_index_within_template[template_atom] for (residue_atom, template_atom) in graph_matcher.mapping.items() }

        return matches

    def _generate_unique_atom_names(self, molecule):
        """
        Generate unique atom names

        Parameters
        ----------
        molecule : openforcefield.topology.Molecule
            The molecule whose atom names are to be modified in-place
        """
        from collections import defaultdict
        element_counts = defaultdict(int)
        for atom in molecule.atoms:
            symbol = atom.element.symbol
            element_counts[symbol] += 1
            atom.name = symbol + str(element_counts[symbol])

    def generator(self, forcefield, residue):
        """
        Residue template generator method to register with simtk.openmm.app.ForceField

        Parameters
        ----------
        forcefield : simtk.openmm.app.ForceField
            The ForceField object to which residue templates and/or parameters are to be added.
        residue : simtk.openmm.app.Topology.Residue
            The residue topology for which a template is to be generated.

        Returns
        -------
        success : bool
            If the generator is able to successfully parameterize the residue, `True` is returned.
            If the generator cannot parameterize the residue, it should return `False` and not modify `forcefield`.

        """
        if self._database_table_name is None:
            raise NotImplementedError('SmallMoleculeTemplateGenerator is an abstract base class and cannot be used directly.')

        from io import StringIO

        # TODO: Refactor to reduce code duplication

        _logger.info(f'Requested to generate parameters for residue {residue}')

        # If a database is specified, check against molecules in the database
        if self._cache is not None:
            with self._open_db() as db:
                table = db.table(self._database_table_name)
                for entry in table:
                    # Skip any molecules we've added to the database this session
                    if entry['smiles'] in self._smiles_added_to_db:
                        continue

                    # See if the template matches
                    from openforcefield.topology import Molecule
                    molecule_template = Molecule.from_smiles(entry['smiles'])
                    print(f"Checking against {entry['smiles']}")
                    if self._match_residue(residue, molecule_template):
                        ffxml_contents = entry['ffxml']

                        # Write to debug file if requested
                        if self.debug_ffxml_filename is not None:
                            with open(self.debug_ffxml_filename, 'w') as outfile:
                                _logger.info(f'writing ffxml to {self.debug_ffxml_filename}')
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
                    with open(self.debug_ffxml_filename, 'w') as outfile:
                        _logger.info(f'writing ffxml to {self.debug_ffxml_filename}')
                        outfile.write(ffxml_contents)

                # Add the parameters and residue definition
                forcefield.loadFile(StringIO(ffxml_contents))
                # If a cache is specified, add this molecule
                if self._cache is not None:
                    with self._open_db() as db:
                        table = db.table(self._database_table_name)
                        _logger.info(f'Writing residue template for {smiles} to cache {self._cache}')
                        record = {'smiles' : smiles, 'ffxml' : ffxml_contents}
                        # Add the IUPAC name for convenience if we can
                        try:
                            record['iupac'] = molecule.to_iupac()
                        except Exception as e:
                            pass
                        # Store the record
                        table.insert(record)
                        self._smiles_added_to_db.add(smiles)

                # Signal success
                return True

        # Report that we have failed to parameterize the residue
        _logger.warning(f'Did not recognize residue {residue.name}; did you forget to call .add_molecules() to add it?')
        return False

################################################################################
# GAFF-specific OpenMM ForceField template generation utilities
################################################################################

class GAFFTemplateGenerator(SmallMoleculeTemplateGenerator):
    """
    OpenMM ForceField residue template generator for GAFF/AM1-BCC using pre-cached openforcefield toolkit molecules.

    One template generator can be registered in multiple OpenMM ForceField objects.

    Examples
    --------

    Create a template generator for GAFF for a single Molecule and register it with ForceField:

    >>> # Define a Molecule using the openforcefield Molecule object
    >>> from openforcefield.topology import Molecule
    >>> molecule = Molecule.from_smiles('c1ccccc1')
    >>> # Create the GAFF template generator
    >>> from openmoltools.forcefield_generators import GAFFTemplateGenerator
    >>> template_generator = GAFFTemplateGenerator(molecules=molecule)
    >>> # Create an OpenMM ForceField
    >>> from simtk.openmm.app import ForceField
    >>> amber_forcefields = ['amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml']
    >>> forcefield = ForceField(*amber_forcefields)
    >>> # Register the template generator
    >>> forcefield.registerTemplateGenerator(template_generator.generator)

    Create a template generator for a specific GAFF version and register multiple molecules:

    >>> molecule1 = Molecule.from_smiles('c1ccccc1')
    >>> molecule2 = Molecule.from_smiles('CCO')
    >>> template_generator = GAFFTemplateGenerator(molecules=[molecule1, molecuel2], forcefield='gaff-2.11')

    You can also add some Molecules later on after the generator has been registered:

    >>> template_generator.add_molecule(molecule)
    >>> template_generator.add_molecules([molecule1, molecule2])

    You can optionally create or use a tiny database cache of pre-parameterized molecules:

    >>> template_generator = GAFFTemplateGenerator(cache='gaff-molecules.json')

    Newly parameterized molecules will be written to the cache, saving time next time!

    You can see which force fields are supported with

    >>> template_generator.INSTALLED_FORCEFIELDS
    ['gaff-1.4', 'gaff-1.8', 'gaff-1.81', 'gaff-2.1', 'gaff-2.11']

    """
    INSTALLED_FORCEFIELDS = ['gaff-1.4', 'gaff-1.8', 'gaff-1.81', 'gaff-2.1', 'gaff-2.11']

    def __init__(self, molecules=None, forcefield=None, cache=None):
        """
        Create a GAFFTemplateGenerator with some openforcefield toolkit molecules

        Requies the openforcefield toolkit: http://openforcefield.org

        Parameters
        ----------
        molecules : openforcefield.topology.Molecule or list, optional, default=None
            Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can be used to construct a Molecule.
            Can also be a list of Molecule objects or objects that can be used to construct a Molecule.
            If specified, these molecules will be recognized and parameterized with antechamber as needed.
            The parameters will be cached in case they are encountered again the future.
        cache : str, optional, default=None
            Filename for global caching of parameters.
            If specified, parameterized molecules will be stored in a TinyDB instance as a JSON file.
        forcefield : str, optional, default=None
            GAFF force field to use, one of ['gaff-1.4', 'gaff-1.8', 'gaff-1.81', 'gaff-2.1', 'gaff-2.11']
            If not specified, the latest GAFF supported version is used.
            GAFFTemplateGenerator.INSTALLED_FORCEFIELDS contains a complete up-to-date list of supported force fields.

        Examples
        --------

        Create a GAFF template generator for a single molecule (benzene, created from SMILES) and register it with ForceField:

        >>> from openforcefield.topology import Molecule
        >>> molecule = Molecule.from_smiles('c1ccccc1')
        >>> from openmoltools.forcefield_generators import GAFFTemplateGenerator
        >>> gaff = GAFFTemplateGenerator(molecules=molecule)
        >>> from simtk.openmm.app import ForceField
        >>> amber_forcefields = ['amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml']
        >>> forcefield = ForceField(*amber_forcefields)
        >>> forcefield.registerTemplateGenerator(gaff)

        The latest GAFF version is used if none is specified.
        You can check which GAFF version is in use with

        >>> gaff.forcefield
        'gaff-2.11'

        Create a template generator for a specific GAFF version for multiple molecules read from an SDF file:

        >>> molecules = Molecule.from_file('molecules.sdf')
        >>> gaff = GAFFTemplateGenerator(molecules=molecules, forcefield='gaff-2.11')

        You can also add molecules later on after the generator has been registered:

        >>> gaff.add_molecules(molecule)
        >>> gaff.add_molecules([molecule1, molecule2])

        To check which GAFF versions are supported, check the `INSTALLED_FORCEFIELDS` attribute:

        >>> print(GAFFTemplateGenerator.INSTALLED_FORCEFIELDS)
        ['gaff-1.4', 'gaff-1.8', 'gaff-1.81', 'gaff-2.1', 'gaff-2.11']

        You can optionally create or use a tiny database cache of pre-parameterized molecules:

        >>> gaff = GAFFTemplateGenerator(cache='gaff-molecules.json', forcefield='gaff-1.80')

        Newly parameterized molecules will be written to the cache, saving time next time!
        """
        # Initialize molecules and cache
        super().__init__(molecules=molecules, cache=cache)

        # Use latest supported GAFF version if none is specified
        if forcefield is None:
            forcefield = self.INSTALLED_FORCEFIELDS[-1]

        # Ensure a valid GAFF version is specified
        if not forcefield in self.INSTALLED_FORCEFIELDS:
            raise ValueError(f"Specified 'forcefield' ({forcefield}) must be one of {self.INSTALLED_FORCEFIELDS}")

        # Store user-specified GAFF version
        self._forcefield = forcefield
        import re
        result = re.match('^gaff-(?P<major_version>\d+)\.(?P<minor_version>\d+)$', forcefield)
        if result is None:
            msg = "'forcefield' must be of form 'gaff-X.Y', where X and Y denote major and minor version\n"
            msg += f"Provided 'forcefield' argument was '{forcefield}'\n"
            msg += f"Supported values are: {self.INSTALLED_FORCEFIELDS}"
            raise ValueError(msg)
        self._gaff_major_version = result['major_version']
        self._gaff_minor_version = result['minor_version']
        self._gaff_version = f'{self._gaff_major_version}.{self._gaff_minor_version}'

        # Track parameters by GAFF version string
        self._database_table_name = forcefield

        # Track which OpenMM ForceField objects have loaded the relevant GAFF parameters
        self._gaff_parameters_loaded = dict()

    @property
    def gaff_version(self):
        """The current GAFF version in use"""
        return self._gaff_major_version + '.' + self._gaff_minor_version

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
        from pkg_resources import resource_filename
        filename = resource_filename('openmmforcefields', os.path.join('ffxml', 'amber', 'gaff', 'dat', f'{self._forcefield}.dat'))
        return filename

    @property
    def gaff_xml_filename(self):
        """File path to the GAFF .ffxml OpenMM force field file"""
        from pkg_resources import resource_filename
        filename = resource_filename('openmmforcefields', os.path.join('ffxml', 'amber', 'gaff', 'ffxml', f'{self._forcefield}.xml'))
        return filename

    def generator(self, forcefield, residue):
        """
        Residue template generator method to register with simtk.openmm.app.ForceField

        Parameters
        ----------
        forcefield : simtk.openmm.app.ForceField
            The ForceField object to which residue templates and/or parameters are to be added.
        residue : simtk.openmm.app.Topology.Residue
            The residue topology for which a template is to be generated.

        Returns
        -------
        success : bool
            If the generator is able to successfully parameterize the residue, `True` is returned.
            If the generator cannot parameterize the residue, it should return `False` and not modify `forcefield`.

        """
        # Load the GAFF parameters if we haven't done so already for this force field
        if not forcefield in self._gaff_parameters_loaded:
            # Instruct the ForceField to load the GAFF parameters
            forcefield.loadFile(self.gaff_xml_filename)
            # Note that we've loaded the GAFF parameters
            self._gaff_parameters_loaded[forcefield] = True

        return super().generator(forcefield, residue)

    def generate_residue_template(self, molecule, residue_atoms=None):
        """
        Generate a residue template and additional parameters for the specified Molecule.

        Parameters
        ----------
        molecules : openforcefield.topology.Molecule or list of Molecules, optional, default=None
            Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can be used to construct a Molecule.
            Can also be a list of Molecule objects or objects that can be used to construct a Molecule.
            If specified, these molecules will be recognized and parameterized with antechamber as needed.
            The parameters will be cached in case they are encountered again the future.
        residue_atoms : list of openforcefield.topology.Atom, optional, default=None
            If specified, the subset of atoms to use in constructing a residue template

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
        _logger.info(f'Generating a residue template for {smiles}')

        # Generate unique atom names
        self._generate_unique_atom_names(molecule)

        # Compute net formal charge
        net_charge = molecule.total_charge
        _logger.debug(f'Total charge is {net_charge}')

        # Compute partial charges if required
        # TODO: Replace this with a call to Molecule API once we have a way to check for user-specified charges
        import numpy as np
        from simtk import unit
        if np.all(molecule.partial_charges / unit.elementary_charge == 0):
            _logger.debug(f'Computing AM1-BCC charges...')
            molecule.compute_partial_charges_am1bcc()
        else:
            _logger.debug(f'Using user-provided charges because partial charges are nonzero...')

        # Geneate a single conformation
        _logger.debug(f'Generating a conformer...')
        molecule.generate_conformers(n_conformers=1)

        # Create temporary directory for running antechamber
        import tempfile
        import os
        tmpdir = tempfile.mkdtemp()
        prefix = 'molecule'
        input_sdf_filename = os.path.join(tmpdir, prefix + '.sdf')
        gaff_mol2_filename = os.path.join(tmpdir, prefix + '.gaff.mol2')
        frcmod_filename = os.path.join(tmpdir, prefix + '.frcmod')

        # Write MDL SDF file for input into antechamber
        molecule.to_file(input_sdf_filename, file_format='sdf')

        # Parameterize the molecule with antechamber (without charging)
        _logger.debug(f'Running antechamber...')
        self._run_antechamber(molecule_filename=input_sdf_filename, input_format='mdl',
                gaff_mol2_filename=gaff_mol2_filename, frcmod_filename=frcmod_filename)

        # Read the resulting GAFF mol2 file atom types
        _logger.debug(f'Reading GAFF atom types...')
        self._read_gaff_atom_types_from_mol2(gaff_mol2_filename, molecule)

        # If residue_atoms = None, add all atoms to the residues
        if residue_atoms == None:
            residue_atoms = [ atom for atom in molecule.atoms ]

        # Modify partial charges so that charge on residue atoms is integral
        # TODO: This may require some modification to correctly handle API changes
        # when openforcefield makes charge quantities consistently unit-bearing or
        # pure numbers.
        _logger.debug(f'Fixing partial charges...')
        from simtk import unit
        residue_charge = 0.0 * unit.elementary_charge
        total_charge = unit.sum(molecule.partial_charges)
        sum_of_absolute_charge = unit.sum(abs(molecule.partial_charges))
        charge_deficit = net_charge * unit.elementary_charge - total_charge
        if sum_of_absolute_charge / unit.elementary_charge > 0.0:
            # Redistribute excess charge proportionally to absolute charge
            molecule.partial_charges += charge_deficit * abs(molecule.partial_charges) / sum_of_absolute_charge

        # Generate additional parameters if needed
        # TODO: Do we have to make sure that we don't duplicate existing parameters already loaded in the forcefield?
        _logger.debug(f'Creating ffxml contents for additional parameters...')
        from inspect import signature # use introspection to support multiple parmed versions
        from io import StringIO
        leaprc = StringIO('parm = loadamberparams %s' % frcmod_filename)
        import parmed
        params = parmed.amber.AmberParameterSet.from_leaprc(leaprc)
        kwargs = {}
        if 'remediate_residues' in signature(parmed.openmm.OpenMMParameterSet.from_parameterset).parameters:
            kwargs['remediate_residues'] = False
        params = parmed.openmm.OpenMMParameterSet.from_parameterset(params, **kwargs)
        ffxml = StringIO()
        kwargs = {}
        if 'write_unused' in signature(params.write).parameters:
            kwargs['write_unused'] = True
        params.write(ffxml, **kwargs)
        ffxml_contents = ffxml.getvalue()

        # Create the residue template
        _logger.debug(f'Creating residue template...')
        from lxml import etree
        root = etree.fromstring(ffxml_contents)
        # Create residue definitions
        residues = etree.SubElement(root, "Residues")
        residue = etree.SubElement(residues, "Residue", name=smiles)
        for atom in molecule.atoms:
            atom = etree.SubElement(residue, "Atom", name=atom.name, type=atom.gaff_type, charge=str(atom.partial_charge / unit.elementary_charge))
        for bond in molecule.bonds:
            if (bond.atom1 in residue_atoms) and (bond.atom2 in residue_atoms):
                bond = etree.SubElement(residue, "Bond", atomName1=bond.atom1.name, atomName2=bond.atom2.name)
            elif (bond.atom1 in residue_atoms) and (bond.atom2 not in residue_atoms):
                bond = etree.SubElement(residue, "ExternalBond", atomName=bond.atom1.name)
            elif (bond.atom1 not in residue_atoms) and (bond.atom2 in residue_atoms):
                bond = etree.SubElement(residue, "ExternalBond", atomName=bond.atom2.name)
        # Render XML into string and append to parameters
        ffxml_contents = etree.tostring(root, pretty_print=True, encoding='unicode')
        _logger.debug(f'ffxml creation complete.')

        return ffxml_contents

    def _run_antechamber(self, molecule_filename, input_format='sdf',
                         gaff_mol2_filename=None, frcmod_filename=None, verbosity=0):
        """Run AmberTools antechamber and parmchk2 to create GAFF mol2 and frcmod files.

        Parameters
        ----------
        molecule_filename : str
            The molecule to be parameterized.
        input_format : str
            antechamber input format for molecule_filename
        gaff_mol2_filename : str, optional, default=None
            Name of GAFF mol2 filename to output.  If None, uses local directory
            and molecule_name
        frcmod_filename : str, optional, default=None
            Name of GAFF frcmod filename to output.  If None, uses local directory
            and molecule_name
        input_format : str, optional, default='mol2'
            Format specifier for input file to pass to antechamber.
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
            gaff_mol2_filename = 'molecule.gaff.mol2'
        if frcmod_filename is None:
            frcmod_filename = 'molecule.frcmod'

        # Build absolute paths for input and output files
        import os
        molecule_filename = os.path.abspath( molecule_filename )
        gaff_mol2_filename = os.path.abspath( gaff_mol2_filename )
        frcmod_filename = os.path.abspath( frcmod_filename )

        def read_file_contents(filename):
            infile = open(filename, 'r')
            contents = infile.read()
            infile.close()
            return contents

        # Use temporary directory context to do this to avoid issues with spaces in filenames, etc.
        import tempfile, subprocess
        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)

            local_input_filename = 'in.' + input_format
            import shutil
            shutil.copy(molecule_filename, local_input_filename)

            # Determine whether antechamber supports -dr [yes/no] option
            cmd = f'antechamber -h | grep dr'
            supports_acdoctor = False
            if ('acdoctor' in subprocess.getoutput(cmd)):
                supports_acdoctor = True

            # Run antechamber without charging (which is done separately)
            cmd = f'antechamber -i {local_input_filename} -fi {input_format} -o out.mol2 -fo mol2 -s {verbosity} -at {self._gaff_major_version}'
            if supports_acdoctor:
                cmd += ' -dr ' + ('yes' if verbosity else 'no')

            _logger.debug(cmd)
            output = subprocess.getoutput(cmd)
            import os
            if not os.path.exists('out.mol2'):
                msg  = "antechamber failed to produce output mol2 file\n"
                msg += "command: %s\n" % cmd
                msg += "output:\n"
                msg += 8 * "----------" + '\n'
                msg += output
                msg += 8 * "----------" + '\n'
                msg += "input:\n"
                msg += 8 * "----------" + '\n'
                msg += read_file_contents(local_input_filename)
                msg += 8 * "----------" + '\n'
                # TODO: Run antechamber again with acdoctor mode on (-dr yes) to get more debug info, if supported
                raise Exception(msg)
            _logger.debug(output)

            # Run parmchk.
            cmd = f"parmchk2 -i out.mol2 -f mol2 -p {self.gaff_dat_filename} -o out.frcmod -s %{self._gaff_major_version}"
            _logger.debug(cmd)
            output = subprocess.getoutput(cmd)
            if not os.path.exists('out.frcmod'):
                msg  = "parmchk2 failed to produce output frcmod file\n"
                msg += "command: %s\n" % cmd
                msg += "output:\n"
                msg += 8 * "----------" + '\n'
                msg += output
                msg += 8 * "----------" + '\n'
                msg += "input mol2:\n"
                msg += 8 * "----------" + '\n'
                msg += read_file_contents('out.mol2')
                msg += 8 * "----------" + '\n'
                raise Exception(msg)
            _logger.debug(output)
            self._check_for_errors(output)

            # Copy back
            shutil.copy( 'out.mol2', gaff_mol2_filename )
            shutil.copy( 'out.frcmod', frcmod_filename )

            os.chdir(cwd)

        return gaff_mol2_filename, frcmod_filename

    def _read_gaff_atom_types_from_mol2(self, gaff_mol2_filename, molecule):
        """
        Read the GAFF atom types specified in an antechamber-generated mol2 file into atom.gaff_type in the specified Molecule

        Parameters
        ----------
        gaff_mol2_filename : str
            The antechamber-generated mol2 file containing GAFF/GAFF2 atom types
        molecule : Molecule
            The Molecule to receive atom types
            The Atom objects within the molecule will have a ``gaff_type`` field added containing the GAFF atom type as a string
        """
        # Read the resulting GAFF mol2 file atom types
        #       1 C1           1.8850    -1.0360    -0.1120 ca         1 MOL       0.000000
        # 012345678901234567890123456789012345678901234567890123456789012345678901234567890
        # 0         1         2         3         4         5         6         7         8
        with open(gaff_mol2_filename, 'r') as infile:
            line = infile.readline()
            # Seek to ATOM block
            while line:
                if line.strip() == '@<TRIPOS>ATOM':
                    break
                line = infile.readline()
            # Read GAFF atom types
            for index, atom in enumerate(molecule.atoms):
                line = infile.readline()
                atom.gaff_type = line[50:58].strip()

        return

    def _check_for_errors(self, outputtext, other_errors=None, ignore_errors=None):
        """Check AMBER package output for the string 'ERROR' (upper or lowercase) and (optionally) specified other strings and raise an exception if it is found (to avoid silent failures which might be noted to log but otherwise ignored).

        Parameters
        ----------
        outputtext : str
            String listing output text from an (AMBER) command which should be checked for errors.
        other_errors : list(str), default None
            If specified, provide strings for other errors which will be chcked for, such as "improper number of arguments", etc.
        ignore_errors: list(str), default None
            If specified, AMBER output lines containing errors but also containing any of the specified strings will be ignored (because, for example, AMBER issues an "ERROR" for non-integer charges in some cases when only a warning is needed).

        Notes
        -----
        If error(s) are found, raise a RuntimeError and attept to print the appropriate errors from the processed text.

        """

        lines = outputtext.split('\n')
        error_lines = []
        for line in lines:
            if 'ERROR' in line.upper():
                error_lines.append( line )
            if not other_errors == None:
                for err in other_errors:
                    if err.upper() in line.upper():
                        error_lines.append( line )

        if not ignore_errors == None and len(error_lines)>0:
            new_error_lines = []
            for ign in ignore_errors:
                ignore = False
                for err in error_lines:
                    if ign in err:
                        ignore = True
                if not ignore:
                    new_error_lines.append( err )
            error_lines = new_error_lines

        if len(error_lines) > 0:
            _logger.info("Unexpected errors encountered running AMBER tool. Offending output:")
            for line in error_lines:
                _logger.info(line)
            raise(RuntimeError("Error encountered running AMBER tool. Exiting."))

        return

################################################################################
# Open Force Field Initiative SMIRNOFF specific OpenMM ForceField template generation utilities
################################################################################

class SMIRNOFFTemplateGenerator(SmallMoleculeTemplateGenerator):
    """
    OpenMM ForceField residue template generator for Open Force Field Initiative SMIRNOFF
    force fields using pre-cached openforcefield toolkit molecules.

    Open Force Field Initiative: http://openforcefield.org
    SMIRNOFF force field specification: https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html

    Examples
    --------

    Create a template generator for a single Molecule using the latest Open Force Field Initiative
    small molecule force field and register it with ForceField:

    >>> # Define a Molecule using the openforcefield Molecule object
    >>> from openforcefield.topology import Molecule
    >>> molecule = Molecule.from_smiles('c1ccccc1')
    >>> # Create the SMIRNOFF template generator
    >>> from openmoltools.forcefield_generators import SMIRNOFFTemplateGenerator
    >>> template_generator = SMIRNOFFTemplateGenerator(molecules=molecule)
    >>> # Create an OpenMM ForceField
    >>> from simtk.openmm.app import ForceField
    >>> amber_forcefields = ['amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml']
    >>> forcefield = ForceField(*amber_forcefields)
    >>> # Register the template generator
    >>> forcefield.registerTemplateGenerator(template_generator.generator)

    Create a template generator for a specific pre-installed SMIRNOFF version ('openff-1.0.0')
    and register multiple molecules:

    >>> molecule1 = Molecule.from_smiles('c1ccccc1')
    >>> molecule2 = Molecule.from_smiles('CCO')
    >>> template_generator = SMIRNOFFTemplateGenerator(molecules=[molecule1, molecule2], forcefield='openff-1.0.0')

    Alternatively, you can specify a local .offxml file in the SMIRNOFF specification:

    >>> template_generator = SMIRNOFFTemplateGenerator(molecules=[molecule1, molecule2], forcefield='mysmirnoff.offxml')

    You can also add some Molecules later on after the generator has been registered:

    >>> template_generator.add_molecule(molecule)
    >>> template_generator.add_molecules([molecule1, molecule2])

    You can optionally create or use a tiny database cache of pre-parameterized molecules:

    >>> template_generator = GAFFTemplateGenerator(cache='gaff-molecules.json')

    Newly parameterized molecules will be written to the cache, saving time next time!

    """
    # TODO: Automatically populate this at import by examining plugin directories in order of semantic version
    INSTALLED_FORCEFIELDS = ['smirnoff99Frosst-1.1.0', 'openff-1.0.0']

    def __init__(self, molecules=None, cache=None, forcefield=None):
        """
        Create a SMIRNOFFTemplateGenerator with some openforcefield toolkit molecules

        Requies the openforcefield toolkit: http://openforcefield.org

        Parameters
        ----------
        molecules : openforcefield.topology.Molecule or list, optional, default=None
            Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can be used to construct a Molecule.
            Can also be a list of Molecule objects or objects that can be used to construct a Molecule.
            If specified, these molecules will be recognized and parameterized with antechamber as needed.
            The parameters will be cached in case they are encountered again the future.
        cache : str, optional, default=None
            Filename for global caching of parameters.
            If specified, parameterized molecules will be stored in a TinyDB instance as a JSON file.
        forcefield : str, optional, default=None
            Name of installed SMIRNOFF force field (without .offxml) or local .offxml filename (with extension).
            If not specified, the latest Open Force Field Initiative release is used.

        Examples
        --------

        Create a GAFF template generator for a single molecule (benzene, created from SMILES) and register it with ForceField:

        >>> from openforcefield.topology import Molecule
        >>> molecule = Molecule.from_smiles('c1ccccc1')
        >>> from openmoltools.forcefield_generators import SMIRNOFFTemplateGenerator
        >>> smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule)
        >>> from simtk.openmm.app import ForceField
        >>> amber_forcefields = ['amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml']
        >>> forcefield = ForceField(*amber_forcefields)

        The latest Open Force Field Initiative release is used if none is specified.

        >>> smirnof.forcefield
        'openff-1.0.0'

        You can check which SMIRNOFF force field filename is in use with

        >>> smirnoff.smirnoff_filename
        '/full/path/to/openff-1.0.0.offxml'

        Create a template generator for a specific SMIRNOFF force field for multiple
        molecules read from an SDF file:

        >>> molecules = Molecule.from_file('molecules.sdf')
        >>> smirnoff = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield='smirnoff99Frosst-1.1.0')

        You can also add molecules later on after the generator has been registered:

        >>> smirnoff.add_molecules(molecules)

        To check which SMIRNOFF versions are supported, check the `INSTALLED_FORCEFIELDS` attribute:

        >>> print(SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS)
        ['smirnoff99Frosst-1.1.0', 'openff-1.0.0']

        You can optionally create or use a cache of pre-parameterized molecules:

        >>> smirnoff = SMIRNOFFTemplateGenerator(cache='smirnoff.json', forcefield='openff-1.0.0')

        Newly parameterized molecules will be written to the cache, saving time next time!
        """
        # Initialize molecules and cache
        super().__init__(molecules=molecules, cache=cache)

        if forcefield is None:
            # Use latest supported Open Force Field Initiative release if none is specified
            forcefield = self.INSTALLED_FORCEFIELDS[-1]
        self._forcefield = forcefield

        # Track parameters by provided SMIRNOFF name
        # TODO: Can we instead use the force field hash, or some other unique identifier?
        self._database_table_name = forcefield

        # Create ForceField object
        import openforcefield.typing.engines.smirnoff
        try:
            filename = forcefield
            if not filename.endswith('.offxml'):
                filename += '.offxml'
            self._smirnoff_forcefield = openforcefield.typing.engines.smirnoff.ForceField(filename)
        except Exception as e:
            print(e)
            raise ValueError(f"Can't find specified SMIRNOFF force field ({forcefield}) in install paths")

        # Delete constraints, if present
        if 'Constraints' in self._smirnoff_forcefield._parameter_handlers:
            del self._smirnoff_forcefield._parameter_handlers['Constraints']

        # Find SMIRNOFF filename
        smirnoff_filename = self._search_paths(filename)
        self._smirnoff_filename = smirnoff_filename

        # Cache a copy of the OpenMM System generated for each molecule for testing purposes
        self._system_cache = dict()

    def _search_paths(self, filename):
        """Search registered openforcefield plugin directories

        Parameters
        ----------
        filename : str
            The filename to find the full path for

        Returns
        -------
        fullpath : str
            Full path to identified file, or None if no file found
        """
        # TODO: Replace this method once there is a public API in the openforcefield toolkit for doing this

        from openforcefield.utils import get_data_file_path
        from openforcefield.typing.engines.smirnoff.forcefield import _get_installed_offxml_dir_paths

        # Check whether this could be a file path
        if isinstance(filename, str):
            # Try first the simple path.
            searched_dirs_paths = ['']
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

    def get_openmm_system(self, molecule):
        """Retrieve the OpenMM System object generated for a particular molecule for testing/validation.

        Parameters
        ----------
        molecule : openforcefield.topology.Molecule
            The Molecule object

        Returns
        -------
        system : simtk.openmm.System or None
            If the Molecule object has already been parameterized by this instance, this molecule is returned.
            Otherwise, None is returned.
        """
        smiles = molecule.to_smiles()
        if smiles in self._system_cache:
            return self._system_cache[smiles]
        else:
            return None

    def generate_residue_template(self, molecule, residue_atoms=None):
        """
        Generate a residue template and additional parameters for the specified Molecule.

        Parameters
        ----------
        molecules : openforcefield.topology.Molecule or list of Molecules, optional, default=None
            Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can be used to construct a Molecule.
            Can also be a list of Molecule objects or objects that can be used to construct a Molecule.
            If specified, these molecules will be recognized and parameterized with antechamber as needed.
            The parameters will be cached in case they are encountered again the future.
        residue_atoms : list of openforcefield.topology.Atom, optional, default=None
            If specified, the subset of atoms to use in constructing a residue template

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
        _logger.info(f'Generating a residue template for {smiles}')

        # Generate unique atom names
        self._generate_unique_atom_names(molecule)

        # Determine which molecules (if any) contain user-specified partial charges that should be used
        # TODO: Replace this with a call to Molecule API once we have a way to check for user-specified charges
        import numpy as np
        from simtk import unit
        charge_from_molecules = None
        if np.all(molecule.partial_charges / unit.elementary_charge == 0):
            charge_from_molecules = [molecule]
            _logger.debug(f'Using user-provided charges because partial charges are nonzero...')

        # Parameterize molecule
        _logger.debug(f'Generating parameters...')
        system = self._smirnoff_forcefield.create_openmm_system(molecule.to_topology(), charge_from_molecules=charge_from_molecules)

        # Transiently cache a copy of the OpenMM System object generated for testing/verification purposes
        self._system_cache[smiles] = system

        # Generate OpenMM ffxml definition for this molecule
        from lxml import etree
        root = etree.Element("ForceField")

        def as_attrib(quantity):
            """Format simtk.unit.Quantity as XML attribute."""
            if isinstance(quantity, str):
                return quantity
            elif isinstance(quantity, float) or isinstance(quantity, int):
                return str(quantity)
            else:
                from simtk import unit
                return str(quantity.value_in_unit_system(unit.md_unit_system))

        # Append unique type names to atoms
        for index, particle in enumerate(molecule.particles):
            setattr(particle, 'typename', f'{smiles}${particle.name}#{index}')

        # Generate atom types
        atom_types = etree.SubElement(root, "AtomTypes")
        for particle_index, particle in enumerate(molecule.particles):
            # Create a new atom type for each atom in the molecule
            paricle_indices = [particle_index]
            atom_type = etree.SubElement(atom_types, "Type", name=particle.typename,
                element=particle.element.symbol, mass=as_attrib(particle.element.mass))
            atom_type.set('class', particle.typename) # 'class' is a reserved Python keyword, so use alternative API

        # Compile forces into a dict
        forces = dict()
        for force in system.getForces():
            force_name = force.__class__.__name__
            if force_name in forces:
                raise Exception("Two instances of force {force_name} appear in System")
            forces[force_name] = force

        def classes(particle_indices):
            """Build a dict of 'class#=typename' for use in creating XML tags for forces.

            Parameters
            ----------
            particle_indices : list of int
                Particle indices for molecule.particles

            Returns
            -------
            classmap : dict of str : str
                Dict of format { 'class1' : typename1, ... }
            """
            return { f'class{class_index+1}' : molecule.particles[particle_index].typename for class_index,particle_index in enumerate(particle_indices) }

        # Round parameters using strings for ease of comparison
        # DEBUG
        #from simtk import unit
        #def round_quantity(quantity):
        #    NDECIMALS = 3
        #    value = quantity.value_in_unit_system(unit.md_unit_system)
        #    value = round(value, NDECIMALS)
        #    return value
        #for particle_index in range(forces['NonbondedForce'].getNumParticles()):
        #    charge, sigma, epsilon = forces['NonbondedForce'].getParticleParameters(particle_index)
        #    forces['NonbondedForce'].setParticleParameters(particle_index, round_quantity(charge), round_quantity(sigma), round_quantity(epsilon))
        #for exception_index in range(forces['NonbondedForce'].getNumExceptions()):
        #    i, j, chargeProd, sigma, epsilon = forces['NonbondedForce'].getExceptionParameters(exception_index)
        #    forces['NonbondedForce'].setExceptionParameters(exception_index, i, j, round_quantity(chargeProd), round_quantity(sigma), round_quantity(epsilon))

        # Lennard-Jones
        # TODO: Get coulomb14scale and lj14scale from SMIRNOFF ForceField object,
        # though this must match the original AMBER values
        nonbonded_types = etree.SubElement(root, "NonbondedForce", coulomb14scale="0.833333", lj14scale="0.5")
        etree.SubElement(nonbonded_types, "UseAttributeFromResidue", name="charge")
        for particle_index in range(forces['NonbondedForce'].getNumParticles()):
            charge, sigma, epsilon = forces['NonbondedForce'].getParticleParameters(particle_index)
            nonbonded_type = etree.SubElement(nonbonded_types, "Atom",
                sigma=as_attrib(sigma), epsilon=as_attrib(epsilon))
            nonbonded_type.set('class', molecule.particles[particle_index].typename) # 'class' is a reserved Python keyword, so use alternative API

        # Bonds
        bond_types = etree.SubElement(root, "HarmonicBondForce")
        particle_indices = [-1]*2
        for bond_index in range(forces['HarmonicBondForce'].getNumBonds()):
            particle_indices[0], particle_indices[1], length, k = forces['HarmonicBondForce'].getBondParameters(bond_index)
            bond_type = etree.SubElement(bond_types, "Bond", **classes(particle_indices),
                length=as_attrib(length), k=as_attrib(k))

        # Angles
        angle_types = etree.SubElement(root, "HarmonicAngleForce")
        particle_indices = [-1]*3
        for angle_index in range(forces['HarmonicAngleForce'].getNumAngles()):
            particle_indices[0], particle_indices[1], particle_indices[2], angle, k = forces['HarmonicAngleForce'].getAngleParameters(angle_index)
            angle_type = etree.SubElement(angle_types, "Angle", **classes(particle_indices),
                angle=as_attrib(angle), k=as_attrib(k))

        # Torsions
        def torsion_tag(particle_indices):
            """Return 'Proper' or 'Improper' depending on torsion type"""
            atoms = [ molecule.particles[particle_index] for particle_index in particle_indices ]
            # TODO: Check to make sure all particles are in fact atoms and not virtual sites
            if atoms[0].is_bonded_to(atoms[1]) and atoms[1].is_bonded_to(atoms[2]) and atoms[2].is_bonded_to(atoms[3]):
                return "Proper"
            else:
                return "Improper"

        # Collect torsions
        torsions = dict()
        for torsion_index in range(forces['PeriodicTorsionForce'].getNumTorsions()):
            particle_indices = [-1]*4
            particle_indices[0], particle_indices[1], particle_indices[2], particle_indices[3], periodicity, phase, k = forces['PeriodicTorsionForce'].getTorsionParameters(torsion_index)
            particle_indices = tuple(particle_indices)
            if particle_indices in torsions.keys():
                torsions[particle_indices].append( (periodicity, phase, k) )
            else:
                torsions[particle_indices] = [ (periodicity, phase, k) ]

        # Create torsion definitions
        torsion_types = etree.SubElement(root, "PeriodicTorsionForce", ordering='smirnoff')
        for particle_indices in torsions.keys():
            params = dict() # build parameter dictionary
            nterms = len(torsions[particle_indices])
            for term in range(nterms):
                periodicity, phase, k = torsions[particle_indices][term]
                params[f'periodicity{term+1}'] = as_attrib(periodicity)
                params[f'phase{term+1}'] = as_attrib(phase)
                params[f'k{term+1}'] = as_attrib(k)
            torsion_type = etree.SubElement(torsion_types, torsion_tag(particle_indices), **classes(particle_indices), **params)

        # TODO: Handle virtual sites
        virtual_sites = [ particle_index for particle_index in range(system.getNumParticles()) if system.isVirtualSite(particle_index) ]
        if len(virtual_sites) > 0:
            raise Exception('Virtual sites are not yet supported')

        # Create residue definitions
        # TODO: Handle non-Atom particles too (virtual sites)
        from simtk import unit
        residues = etree.SubElement(root, "Residues")
        residue = etree.SubElement(residues, "Residue", name=smiles)
        for particle_index, particle in enumerate(molecule.particles):
            charge, sigma, epsilon = forces['NonbondedForce'].getParticleParameters(particle_index)
            atom = etree.SubElement(residue, "Atom", name=particle.name, type=particle.typename, charge=as_attrib(charge))
        for bond in molecule.bonds:
            bond = etree.SubElement(residue, "Bond", atomName1=bond.atom1.name, atomName2=bond.atom2.name)

        # Render XML into string
        ffxml_contents = etree.tostring(root, pretty_print=True, encoding='unicode')

        return ffxml_contents
