"""
Small molecule residue template generator abstract base class.

"""

################################################################################
# IMPORTS
################################################################################

import json
from abc import ABC, abstractmethod

################################################################################
# LOGGER
################################################################################

import logging
_logger = logging.getLogger("openmmforcefields.generators.molecule")

################################################################################
# Force field generators
################################################################################

class SmallMoleculeTemplateGenerator(object):
    """
    Abstract base class for OpenMM ForceField residue template generators for
    small molecules using pre-cached openforcefield toolkit molecules.

    Note that this base class cannot generate parameters on its own.
    It only contains useful methods for aiding in molecule parameterization.

    """
    def __init__(self, molecules=None, cache=None):
        """
        Create a SmallMoleculeTemplateGenerator capable of generating parameters
        for the atom type centric simtk.openmm.app.ForceField class given
        openforcefield toolkit Molecule objects.

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
            If specified, parameterized molecules will be stored in a TinyDB instance.
            Note that no checking is done to determine this cache was created with the same GAFF version.

        """
        # Add oemols to the dictionary
        self._molecules = dict()
        self.add_molecules(molecules)

        self._cache = cache
        self._smiles_added_to_db = set() # set of SMILES added to the database this session

    # TODO: Replace this encoder/decoder logic when openmm objects are properly serializable
    class _JSONEncoder(json.JSONEncoder):
        def default(self, o):
            from simtk.openmm.app import ForceField, Element
            if isinstance(o, ForceField._TemplateData):
                s = {'_type' : '_TemplateData'}
                s.update(o.__dict__)
                return s
            elif isinstance(o, ForceField._TemplateAtomData):
                s = {'_type' : '_TemplateAtomData'}
                s.update(o.__dict__)
                return s
            elif isinstance(o, Element):
                return {'_type' : 'Element', 'atomic_number' : o.atomic_number}
            else:
                return super(SmallMoleculeTemplateGenerator._JSONEncoder, self).default(o)

    class _JSONDecoder(json.JSONDecoder):
        def __init__(self, *args, **kwargs):
            json.JSONDecoder.__init__(self, object_hook=self.object_hook, *args, **kwargs)

        def object_hook(self, obj):
            from simtk.openmm.app import ForceField, Element
            if '_type' not in obj:
                return obj
            type = obj['_type']
            if type == '_TemplateData':
                template = ForceField._TemplateData.__new__(ForceField._TemplateData)
                del obj['_type']
                template.__dict__ = obj
                return template
            if type == '_TemplateAtomData':
                atom = ForceField._TemplateAtomData.__new__(ForceField._TemplateAtomData)
                del obj['_type']
                atom.__dict__ = obj
                return atom
            elif type == 'Element':
                return Element.getByAtomicNumber(obj['atomic_number'])
            return obj

    def add_molecules(self, molecules=None):
        """
        Add specified list of Molecule objects to cached molecules that will be
        recognized when parameterizing residues with missing parameters.

        Parameters
        ----------
        molecules : openforcefield.topology.Molecule or list, optional, default=None
            Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can be used to construct a Molecule.
            Can also be a list of Molecule objects or objects that can be used to construct a Molecule.
            If specified, these molecules will be recognized and parameterized with antechamber as needed.
            The parameters will be cached in case they are encountered again the future.

        Examples
        --------
        Add some Molecules later on after the generator has been registered:

        >>> forcefield.add_molecules(mol)
        >>> forcefield.add_molecules([mol1, mol2])

        """
        # Return if empty
        if not molecules:
            return

        # Ensure oemols is iterable
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

        # Create graph for residue being queried against template
        residue_graph = nx.Graph()
        for atom in residue.atoms():
            residue_graph.add_node(atom, element=atom.element.atomic_number, number_of_external_bonds=number_of_external_bonds[atom])
        for bond in residue.internal_bonds():
            residue_graph.add_edge(bond[0], bond[1])

        # Create graph for template molecule
        template_graph = nx.Graph()
        for atom_index, atom in enumerate(molecule_template.atoms):
            template_graph.add_node(atom_index, element=atom.atomic_number, number_of_external_bonds=0)
        for (atom1, atom2) in molecule_template.bonds:
            template_graph.add_edge(bond.atom1_index, bond.atom2_index)

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

    def generateResidueTemplate(molecule):
        """
        Abstract base method to generate an residue template for simtk.openmm.app.ForceField.

        This method must be implemented by all subclasses.

        Parameters
        ----------
        molecule : openforcefield.topology.Molecule
            Molecule for which new residue template (and optionally, new parameters) are to be generated.

        Returns
        -------
        template : simtk.openmm.app.forcefield._TemplateData
            Residue template for ForceField using atom types and parameters from `gaff.xml` or `gaff2.xml`.
        additional_parameters_ffxml : str
            Contents of ForceField `ffxml` file defining additional parameters.

        Notes
        -----
        New atom types cannot be defined.

        """
        raise NotImplementedError()

    def _open_db():
        """
        """
        from tinydb import TinyDB
        tinydb_kwargs = { 'sort_keys' : True, 'indent' : 4, 'separators' : (',', ': ') } # for pretty-printing
        db = TinyDB(self._cache, **tinydb_kwargs)
        return db

    def generator(self, forcefield, residue, structure=None):
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
        from io import StringIO

        # If a database is specified, check against molecules in the database
        if self._cache is not None:
            db = self._open_db()
            for entry in db:
                # Skip any molecules we've added to the database this session
                if entry['smiles'] in self._smiles_added_to_db:
                    continue

                # See if the template matches
                from openforcefield.topology import Molecule
                molecule_template = Molecule.from_smiles(entry['smiles'])
                if self._match_residue(residue, molecule_template):
                    # Register the template
                    template = self._JSONDecoder().decode(entry['template'])
                    forcefield.registerResidueTemplate(template)
                    # Add the parameters
                    # TODO: Do we have to worry about parameter collisions?
                    forcefield.loadFile(StringIO(entry['ffxml']))
                    # Signal success
                    return True

        # Check against the molecules we know about
        for smiles, molecule in self._molecules.items():
            # See if the template matches
            if self._match_residue(residue, molecule):
                # Generate template and parameters.
                [template, ffxml] = self.generateResidueTemplate(molecule)
                # Register the template
                forcefield.registerResidueTemplate(template)
                # Add the parameters
                # TODO: Do we have to worry about parameter collisions?
                #       What happens if two residues contain the same additional parameter?
                forcefield.loadFile(StringIO(ffxml))
                # If a cache is specified, add this molecule
                if self._cache is not None:
                    print('Writing {} to cache'.format(smiles))
                    record = {'smiles' : smiles, 'template' : self._JSONEncoder().encode(template), 'ffxml' : ffxml}
                    # Add the IUPAC name for convenience if we can
                    try:
                        record['iupac'] = molecule.to_iupac()
                    except Exception as e:
                        pass
                    # Store the record
                    db.insert(record)
                    self._smiles_added_to_db.add(smiles)
                    db.close()

                # Signal success
                return True

        # Make sure to close database
        if self._cache is not None:
            db.close()

        # Report that we have failed to parameterize the residue
        logger.warn("Didn't know how to parameterize residue {}".format(residue.name))
        return False
