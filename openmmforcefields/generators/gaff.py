"""
Residue template generator for the AMBER GAFF1/2 small molecule force fields.

"""

################################################################################
# IMPORTS
################################################################################

################################################################################
# LOGGER
################################################################################

import logging
_logger = logging.getLogger("openmmforcefields.generators.gaff")

################################################################################
# GAFF-specific force field generation utilities
################################################################################

class IncompatibleGAFFVersion(RuntimeError):
    """
    A cache file has been opened with an incompatible GAFF version.
    """
    pass

class GAFFTemplateGenerator(object):
    """
    OpenMM ForceField residue template generator for GAFF/AM1-BCC using pre-cached openforcefield toolkit molecules.

    Examples
    --------

    Create a template generator for GAFF for a single Molecule and register it with ForceField:

    >>> from openmoltools.forcefield_generators import GAFFTemplateGenerator
    >>> template_generator = GAFFTemplateGenerator(molecules=mol)
    >>> from simtk.openmm.app import ForceField
    >>> forcefield = ForceField('amber14-all.xml', 'tip3p.xml')
    >>> forcefield.registerTemplateGenerator(template_generator.generator)
    >>> forcefield.loadFile(template_generator.gaff_xml_filename)

    Create a template generator for a specific GAFF version and register multiple molecules:

    >>> template_generator = GAFFTemplateGenerator(molecules=[mol1, mol2], gaff_version='2.11')

    You can also add some Molecules later on after the generator has been registered:

    >>> forcefield.add_molecule(mol)
    >>> forcefield.add_molecules([mol1, mol2])

    You can optionally create or use a tiny database cache of pre-parameterized molecules:

    >>> template_generator = GAFFTemplateGenerator(cache='gaff-molecules.json')

    Newly parameterized molecules will be written to the cache, saving time next time!

    """
    SUPPORTED_GAFF_VERSIONS = ['1.4', '1.8', '1.81', '2.1', '2.11']

    def __init__(self, molecules=None, cache=None, gaff_version='2.11'):
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
            If specified, parameterized molecules will be stored in a TinyDB instance.
            Note that no checking is done to determine this cache was created with the same GAFF version.
        gaff_version : str, default = '2.11'
            GAFF version to use. One of ('1.4', '1.8', '1.81', '2.1', '2.11')

        Examples
        --------

        Create a template generator for GAFF for a single OEMol and register it with ForceField:

        >>> from openmoltools.forcefield_generators import GAFFTemplateGenerator
        >>> template_generator = GAFFTemplateGenerator(molecules=mol)

        Create a template generator for GAFF 1.81 for multiple molecules:

        >>> template_generator = GAFFTemplateGenerator(molecules=[mol1, mol2], gaff_version='1.81')

        You can optionally create or use a tiny database cache of pre-parameterized molecules:

        >>> template_generator = GAFFTemplateGenerator(cache='gaff-molecules.json')

        """
        if not gaff_version in SUPPORTED_GAFF_VERSIONS:
            raise Exception(f"Error: gaff_version must be one of {self.SUPPORTED_GAFF_VERSIONS}")
        self._gaff_version = gaff_version
        self._gaff_major_version, self._gaff_minor_version = gaff_version.split('.')

        # Add oemols to the dictionary
        self._molecules = dict()
        self.add_molecules(molecules)

        self._cache = cache
        self._smiles_added_to_db = set() # set of SMILES added to the database this session

        # Check GAFF version compatibility with cache
        # TODO: Add overwrite_cache optional argument to overwrite in this case?
        if self._cache is not None:
            db = self._open_db()
            if 'gaff_version' in db:
                if db['gaff_version'] != gaff_version:
                    raise IncompatibleGAFFVersion(f"Specified cache uses GAFF {db['gaff_version']} but this generator constructed for GAFF {gaff_version}")
            else:
                # Record which GAFF version we're using
                db['gaff_version'] = gaff_version
            db.close()

    @property
    def gaff_dat_filename(self):
        """Return the path to the GAFF dat file corresponding to the requested GAFF version.

        Returns
        -------
        filename : str
            The full path to the corresponding GAFF ffxml file
        """
        from openmmforcefields.utils import get_data_filename
        filename = get_data_filename(f'gaff-{self._gaff_version}.dat')
        return filename

    @property
    def gaff_xml_filename(self):
        """Return the path to the GAFF ffxml file corresponding to the requested GAFF version.

        Returns
        -------
        filename : str
            The full path to the corresponding GAFF ffxml file
        """
        # TODO: Should we return a StringIO instead?
        from openmmforcefields.utils import get_data_filename
        filename = get_data_filename(f'gaff-{self._gaff_version}.xml')
        return filename

    def _open_db():
        """
        """
        from tinydb import TinyDB
        tinydb_kwargs = { 'sort_keys' : True, 'indent' : 4, 'separators' : (',', ': ') } # for pretty-printing
        db = TinyDB(self._cache, **tinydb_kwargs)
        return db

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
        >>> forcefield.add_oemols(mol1)
        >>> forcefield.add_oemols([mol2, mol3])

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
        for (atom1, atom2) in molecule_template.bonds:
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

    def _getoutput(cmd):
        """Compatibility function to substitute deprecated commands.getoutput in Python2.7"""
        import subprocess
        try:
            out = subprocess.getoutput(cmd)
        except AttributeError:
            out = subprocess.Popen(cmd, shell=True, stderr=subprocess.STDOUT,
                                   stdout=subprocess.PIPE).communicate()[0]
        try:
            return str(out.decode())
        except:
            return str(out)

    def generate_residue_template(molecule, residue_atoms=None):
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
        _logger.info(f'Generating a residue template for {smiles}')

        # TODO: Generates unique atom names if they do not exist (or have Molecule do this)

        # Compute net formal charge
        net_charge = molecule.total_charge

        # Compute partial charge
        molecule.compute_partial_charges_am1bcc()

        # Geneate a single conformation
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
        molecule.to_file(input_sdf_filename)

        # Parameterize the molecule with antechamber (without charging)
        self._run_antechamber(molecule_sdf_filename=input_sdf_filename, input_format='mdl',
                gaff_mol2_filename=gaff_mol2_filename, frcmod_filename=frcmod_filename)

        # Read the resulting GAFF mol2 file atom types
        self._read_gaff_atom_types_from_mol2(gaff_mol2_filename, molecule)

        # If residue_atoms = None, add all atoms to the residues
        if residue_atoms == None:
            residue_atoms = [ atom for atom in molecule.atoms ]

        # Modify partial charges so that charge on residue atoms is integral
        residue_charge = 0.0
        sum_of_absolute_charge = 0.0
        for atom in residue_atoms:
            charge = atom.partial_charge
            residue_charge += charge
            sum_of_absolute_charge += abs(charge)
        excess_charge = residue_charge - molecule.net_charge
        # Redistribute excess charge proportionally to absolute charge
        if sum_of_absolute_charge == 0.0:
            sum_of_absolute_charge = 1.0
        for atom in residue_atoms:
            charge = atom.partial_charge
            atom.partial_charge = charge + excess_charge * (abs(charge) / sum_of_absolute_charge)

        # Generate additional parameters if needed
        # TODO: Do we have to make sure that we don't duplicate existing parameters already loaded in the forcefield?
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
        from simtk.openmm.app import ForceField, Element
        ffxml_contents += '<Residues>\n'
        ffxml_contents += f'  <Residue name="{molecule.to_smiles()}">'
        for (index, atom) in enumerate(molecule.atoms):
            element = Element.getByAtomicNumber(atom.atomic_number)
            ffxml_contents += f'    <Atom name="{atom.name}" type="{atom.gaff_type}" charge="{atom.partial_charge}"/>\n'
        for bond in molecule.bonds:
            if (bond.atom1 in residue_atoms) and (bond.atom2 in residue_atoms):
                ffxml_contents += f'    <Bond atomName1="{bond.atom1.name}" atomName1="{bond.atom2.name}"/>\n'
            elif (bond.atom1 in residue_atoms) and (bond.atom2 not in residue_atoms):
                ffxml_contents += f'    <ExternalBond atomName="{bond.atom1.name}"/>\n'
            elif (bond.atom1 not in residue_atoms) and (bond.atom2 in residue_atoms):
                ffxml_contents += f'    <ExternalBond atomName="{bond.atom2.name}"/>\n'
        ffxml_contents += f'  </Residue>'
        ffxml_contents += "</Residues>\n"

        return ffxml_contents

    def _run_antechamber(molecule_filename, input_format='sdf',
        gaff_mol2_filename=None, frcmod_filename=None,
        input_format='mol2', resname=False, log_debug_output=False, verbosity=0):
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
        resname : bool, optional, default=False
            Set the residue name used within output files to molecule_name
        log_debug_output : bool, optional, default=False
            If true, will send output of tleap to_logger.
        gaff_version : str, default = '2.11'
            GAFF version to use. One of ('1.4', '1.8', '1.81', '2.1', '2.11')
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
        gaff_mol2_filename = os.path.abspath( gaff_mol2_filename )
        frcmod_filename = os.path.abspath( frcmod_filename )
        input_filename = os.path.abspath( input_filename )

        def read_file_contents(filename):
            infile = open(filename, 'r')
            contents = infile.read()
            infile.close()
            return contents

        # Use temporary directory context to do this to avoid issues with spaces in filenames, etc.
        with enter_temp_directory():
            local_input_filename = 'in.' + input_format
            import shutil
            shutil.copy(molecule_filename, local_input_filename)

            # Run antechamber without charging (which is done separately)
            cmd = f'antechamber -i {local_input_filename} -fi {input_format} -o {gaff_mol2_filename} -fo mol2 -s {verbosity} -at {self._gaff_major_version}'

            if log_debug_output: _logger.debug(cmd)
            output = self._getoutput(cmd)
            import os
            if not os.path.exists('out.mol2'):
                msg  = "antechamber failed to produce output mol2 file\n"
                msg += "command: %s\n" % cmd
                msg += "output:\n"
                msg += 8 * "----------" + '\n'
                msg += output
                msg += 8 * "----------" + '\n'
                msg += "input mol2:\n"
                msg += 8 * "----------" + '\n'
                msg += read_file_contents(local_input_filename)
                msg += 8 * "----------" + '\n'
                raise Exception(msg)
            if log_debug_output: _logger.debug(output)

            # Run parmchk.
            cmd = f"parmchk2 -i out.mol2 -f mol2 -p {self.gaff_dat_filename} -o out.frcmod -s %{self._gaff_major_version}"
            if log_debug_output: _logger.debug(cmd)
            output = self._getoutput(cmd)
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
            if log_debug_output: _logger.debug(output)
            self._check_for_errors(output)

            # Copy back
            shutil.copy( 'out.mol2', gaff_mol2_filename )
            shutil.copy( 'out.frcmod', frcmod_filename )

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
            # Read GAFF atom types
            for index, atom in enumerate(molecule.n_atoms):
                line = infile.readline()
                atom.gaff_type = line[50:58].strip()

    def _check_for_errors(outputtext, other_errors=None, ignore_errors=None ):
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
            print("Unexpected errors encountered running AMBER tool. Offending output:")
            for line in error_lines: print(line)
            raise(RuntimeError("Error encountered running AMBER tool. Exiting."))

        return

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
                    # Add parameters and residue template for this residue
                    # TODO: What can we do to avoid parameter or template collisions?
                    forcefield.loadFile(StringIO(entry['ffxml']))
                    # Signal success
                    return True

        # Check against the molecules we know about
        for smiles, molecule in self._molecules.items():
            # See if the template matches
            if self._match_residue(residue, molecule):
                # Generate template and parameters.
                ffxml_contents = self.generate_residue_template(molecule)
                # Add the parameters and residue definition
                # TODO: Do we have to worry about parameter collisions?
                #       What happens if two residues contain the same additional parameter?
                forcefield.loadFile(StringIO(ffxml_contents))
                # If a cache is specified, add this molecule
                if self._cache is not None:
                    print('Writing {} to cache'.format(smiles))
                    record = {'smiles' : smiles, 'ffxml' : ffxml_contents}
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
