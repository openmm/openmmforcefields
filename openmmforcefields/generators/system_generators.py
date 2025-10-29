"""
System generators that build an OpenMM System object from a Topology object.
"""

from openmmforcefields.utils import classproperty

################################################################################
# LOGGER
################################################################################

import logging

_logger = logging.getLogger("openmmforcefields.generators.system_generators")

################################################################################
# System generator base class
################################################################################


class SystemGenerator:
    """
    Common interface for generating OpenMM Systems from OpenMM Topology objects
    that may contain both biopolymers (with parameters provided by OpenMM) and small molecules
    (with parameters provided by residue template generators).

    Currently, this class supports

    * GAFF, via ``GAFFTemplateGenerator``: see ``GAFFTemplateGenerator.INSTALLED_FORCEFIELDS``
    * SMIRNOFF, via ``SMIRNOFFTemplateGenerator``: see ``SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS``

    .. todo ::

       Once Open Force Field Topology objects support residue definitions, we will also be able
       to support Open Force Field Topology objects (which carry their own Molecule definitions).

    Parameters
    ----------
    forcefield : openmm.app.ForceField
        The ForceField object used to create new System objects.
        New ffxml files can be read in at any time.
    forcefield_kwargs : dict
        Keyword arguments fed to ``openmm.app.ForceField.createSystem()`` during System generation.
        These keyword arguments can be modified at any time.
    periodic_forcefield_kwargs : dict
        Keyword arguments fed to ``openmm.app.ForceField.createSystem()`` during System generation for
        periodic systems.
        These keyword arguments can be modified at any time.
    nonperiodic_forcefield_kwargs : dict
        Keyword arguments fed to ``openmm.app.ForceField.createSystem()`` during System generation for non-periodic
        systems.
        These keyword arguments can be modified at any time.
    barostat : openmm.MonteCarloBarostat
        If not None, this container holds the barostat parameters to use for newly created System objects.
    molecules : openff.toolkit.topology.Molecule or list, optional, default=None
        Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can be used to
        construct a Molecule.
        Can also be a list of Molecule objects or objects that can be used to construct a Molecule.
        If specified, these molecules will be recognized and parameterized as needed.
        The parameters will be cached in case they are encountered again the future.
    cache : filename, optional, default=None
        If not None, filename for caching small molecule residue templates.
    postprocess_system : method
        If not None, this method will be called as ``system = postprocess_system(system)`` to post-process the System
        object for create_system(topology) before it is returned.
    """

    def __init__(
        self,
        forcefields=None,
        small_molecule_forcefield="openff-2.2.0",
        forcefield_kwargs=None,
        nonperiodic_forcefield_kwargs=None,
        periodic_forcefield_kwargs=None,
        template_generator_kwargs=None,
        barostat=None,
        molecules=None,
        cache=None,
        postprocess_system=None,
    ):
        """
        This is a utility class to generate OpenMM Systems from Open Force Field Topology objects using AMBER
        protein force fields and GAFF small molecule force fields.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        forcefields : list of str, optional, default=None
            List of the names of ffxml files that will be used in System creation for the biopolymer.
        small_molecule_forcefield : str, optional, default='openff-2.2.0'
            Small molecule force field to use.
            Must be supported by one of the registered template generators:
                [GAFFTemplateGenerator, SMIRNOFFTemplateGenerator]
            Supported GAFF force fields include 'gaff-2.2.20', 'gaff-2.11', and others.
            (See ``GAFFTemplateGenerator.INSTALLED_FORCEFIELDS`` for a complete list.)
            Supported SMIRNOFF force fields include all installed force fields from Parsley and Sage lines, such as 'openff-1.0.0' and 'openff-2.2.1'.
            (See ``SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS`` for a complete list.)
        forcefield_kwargs : dict, optional, default=None
            Keyword arguments to be passed to ``openmm.app.ForceField.createSystem()`` during ``System``
            object creation.
        nonperiodic_forcefield_kwargs : dict, optional, default={'nonbondedMethod' : NoCutoff}
            Keyword arguments added to forcefield_kwargs when the Topology is non-periodic.
        periodic_forcefield_kwargs : NonbondedMethod, optional, default={'nonbondedMethod' : PME}
            Keyword arguments added to forcefield_kwargs when the Topology is periodic.
        template_generator_kwargs : dict, optional, default=None
            The value of an extra keyword argument called `template_generator_kwargs` to be passed to subclasses of
            ``openmmforcefields.generators.template_generators.SmallMoleculeTemplateGenerator``.
            Currently only used for ``openmmforcefields.generators.template_generators.EspalomaTemplateGenerator``.
        barostat : openmm.MonteCarloBarostat, optional, default=None
            If not None, a new ``MonteCarloBarostat`` with matching parameters (but a different random number seed)
            will be created and
            added to each newly created ``System``.
        molecules : openff.toolkit.topology.Molecule or list, optional, default=None
            Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can be used
            to construct a Molecule.
            Can also be a list of Molecule objects or objects that can be used to construct a Molecule.
            If specified, these molecules will be recognized and parameterized as needed.
            The parameters will be cached in case they are encountered again the future.
        cache : filename, optional, default=None
            If not None, filename for caching small molecule residue templates.
        postprocess_system : method, optiona, default=None
            If not None, this method will be called as ``system = postprocess_system(system)`` to post-process the
            System object for create_system(topology) before it is returned.

        Examples
        --------

        Here's an example that uses GAFF 2.2.20 along with the new ``ff14SB`` generation of AMBER force fields
        (and compatible solvent models) to generate an OpenMM ``System`` object from an
        `Open Force Field Topology <https://open-forcefield-toolkit.readthedocs.io/en/latest/api/generated/openff.toolkit.topology.Topology.html#openff.toolkit.topology.Topology>`_ object:

        >>> # Define the keyword arguments to feed to ForceField
        >>> from openmm import unit
        >>> from openmm import app
        >>> # Define standard OpenMM biopolymer and solvent force fields to use

        To initialize the ``SystemGenerator``, we specify the OpenMM force fields, the small molecule force field, and
        any ``kwargs`` to be fed to the OpenMM ``openmm.app.ForceField.createSystem()`` method:

        >>> from openmmforcefields.generators import SystemGenerator
        >>> amber_forcefields = [
        ...     "amber/protein.ff14SB.xml",
        ...     "amber/tip3p_standard.xml",
        ...     "amber/tip3p_HFE_multivalent.xml",
        ... ]
        >>> small_molecule_forcefield = "gaff-2.2.20"
        >>> forcefield_kwargs = {
        ...     "constraints": app.HBonds,
        ...     "rigidWater": True,
        ...     "removeCMMotion": False,
        ...     "hydrogenMass": 4 * unit.amu,
        ... }
        >>> system_generator = SystemGenerator(
        ...     forcefields=amber_forcefields,
        ...     small_molecule_forcefield=small_molecule_forcefield,
        ...     forcefield_kwargs=forcefield_kwargs,
        ... )

        If the ``cache`` argument is specified, parameterized molecules are cached in the corresponding file.

        >>> cache = 'db.json'
        >>> system_generator = SystemGenerator(forcefields=amber_forcefields, small_molecule_forcefield='gaff-2.2.20', forcefield_kwargs=forcefield_kwargs, cache=cache)  # doctest: +SKIP

        To use a barostat, you need to define a barostat whose parameters will be copied into each system (with a different random number seed):

        >>> import openmm
        >>> from openmm import unit
        >>> pressure = 1.0 * unit.atmospheres
        >>> temperature = 298.0 * unit.kelvin
        >>> frequency = 25 # steps
        >>> system_generator.barostat = openmm.MonteCarloBarostat(pressure, temperature, frequency)

        Now, you can create an OpenMM ``System`` object from an OpenMM ``Topology`` object and a list of OpenFF ``Molecule`` objects

        >>> from openff.toolkit import Molecule, Topology
        >>> molecules = [Molecule.from_smiles(smiles) for smiles in ["CCO", "c1ccccc1"]]
        >>> openmm_topology = Topology.from_molecules(molecules).to_openmm()
        >>> system = system_generator.create_system(openmm_topology, molecules=molecules)

        Parameters for multiple force fields can be held in the same cache file.

        To use the `Open Force Field 'openff-1.0.0' ("Parsley") force field <https://openforcefield.org/news/introducing-openforcefield-1.0/>`_ instead,
        simply change the ``small_molecule_forcefield`` parameter to one of the supported ``GAFFTemplateGenerator.INSTALLED_FORCEFIELDS``:

        >>> small_molecule_forcefield = 'openff-1.0.0'
        >>> system_generator = SystemGenerator(forcefields=amber_forcefields, small_molecule_forcefield=small_molecule_forcefield, forcefield_kwargs=forcefield_kwargs)

        For debugging convenience, you can also turn _off_ specific interactions during system creation, such as particle charges:

        >>> system_generator.particle_charges = False # will cause particle charges to be set to zero
        >>> system_generator.exception_charges = False # will zero out all 1-4 charge interactions
        >>> system_generator.particle_epsilons = False # will zero out Lennard-Jones particle-particle interactions
        >>> system_generator.particle_exceptions = False # will zero out all 1-4 Lennard-Jones interactions
        >>> system_generator.torsions = False # will zero out all torsion terms
        """  # noqa

        # Initialize
        self.barostat = barostat  # barostat to copy, or None if no barostat is to be added

        # Post-creation system transformations
        self.particle_charges = True  # include particle charges
        self.exception_charges = True  # include electrostatics nonzero exceptions
        self.particle_epsilons = True  # include LJ particles
        self.exception_epsilons = True  # include LJ nonzero exceptions
        self.torsions = True  # include torsions

        # Method to use for postprocessing system
        self.postprocess_system = postprocess_system

        # Create OpenMM ForceField object
        forcefields = forcefields if (forcefields is not None) else list()
        from openmm import app

        self.forcefield = app.ForceField(*forcefields)

        # Cache force fields and settings to use
        self.forcefield_kwargs = forcefield_kwargs if forcefield_kwargs is not None else dict()
        self.nonperiodic_forcefield_kwargs = (
            nonperiodic_forcefield_kwargs
            if nonperiodic_forcefield_kwargs is not None
            else {"nonbondedMethod": app.NoCutoff}
        )
        self.periodic_forcefield_kwargs = (
            periodic_forcefield_kwargs if periodic_forcefield_kwargs is not None else {"nonbondedMethod": app.PME}
        )
        self.template_generator_kwargs = template_generator_kwargs

        # Raise an exception if nonbondedForce is specified in forcefield_kwargs
        if "nonbondedMethod" in self.forcefield_kwargs:
            raise ValueError(
                """nonbondedMethod cannot be specified in forcefield_kwargs;
                 must be specified in either periodic_forcefield_kwargs (if it should be applied to periodic systems)
                 or nonperiodic_forcefield_kwargs (if it should be applied to non-periodic systems)"""
            )

        # Create and cache a residue template generator
        from openmmforcefields.generators.template_generators import (
            SmallMoleculeTemplateGenerator,
        )

        self.template_generator = None
        if small_molecule_forcefield is not None:
            for template_generator_cls in SmallMoleculeTemplateGenerator.__subclasses__():
                try:
                    _logger.debug(f"Trying {template_generator_cls.__name__} to load {small_molecule_forcefield}")
                    self.template_generator = template_generator_cls(
                        forcefield=small_molecule_forcefield,
                        cache=cache,
                        template_generator_kwargs=self.template_generator_kwargs,
                    )
                    break
                except (ValueError,) as e:
                    _logger.debug(f"  {template_generator_cls.__name__} cannot load {small_molecule_forcefield}")
                    _logger.debug(e)
            if self.template_generator is None:
                msg = (
                    "No registered small molecule template generators could load force field "
                    f"'{small_molecule_forcefield}'\n"
                    "Available installed force fields are:\n"
                )
                for template_generator_cls in SmallMoleculeTemplateGenerator.__subclasses__():
                    msg += f"  {template_generator_cls.__name__}: {template_generator_cls.INSTALLED_FORCEFIELDS}\n"
                raise ValueError(msg)
            self.forcefield.registerTemplateGenerator(self.template_generator.generator)

        # Inform the template generator about any specified molecules
        self.add_molecules(molecules)

    @classproperty
    def SMALL_MOLECULE_FORCEFIELDS(cls):
        """Return a listof available small molecule force fields"""
        forcefields = list()
        from openmmforcefields.generators.template_generators import (
            SmallMoleculeTemplateGenerator,
        )

        for template_generator_cls in SmallMoleculeTemplateGenerator.__subclasses__():
            forcefields += template_generator_cls.INSTALLED_FORCEFIELDS
        return forcefields

    def add_molecules(self, molecules):
        """
        Add molecules to registered template generator

        Parameters
        ----------
        molecules : openff.toolkit.topology.Molecule or list, optional, default=None
            Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can
            be used to construct a Molecule.
            Can also be a list of Molecule objects or objects that can be used to construct a Molecule.
            If specified, these molecules will be recognized and parameterized as needed.
            The parameters will be cached in case they are encountered again the future.
        """

        if self.template_generator is None:
            raise ValueError(
                "You must have a small molecule residue template generator registered to add small molecules"
            )

        self.template_generator.add_molecules(molecules)

    def _modify_forces(self, system):
        """
        Add barostat and modify forces if requested.
        """
        # Add barostat if requested and the system uses periodic boundary conditions
        if (self.barostat is not None) and system.usesPeriodicBoundaryConditions():
            import numpy as np
            import copy

            MAXINT = np.iinfo(np.int32).max

            # Get the barostat
            barostat = copy.deepcopy(self.barostat)
            seed = np.random.randint(MAXINT)
            barostat.setRandomNumberSeed(seed)
            system.addForce(barostat)

        # Modify forces if requested
        for force in system.getForces():
            if force.__class__.__name__ == "NonbondedForce":
                for index in range(force.getNumParticles()):
                    charge, sigma, epsilon = force.getParticleParameters(index)
                    if not self.particle_charges:
                        charge *= 0
                    if not self.particle_epsilons:
                        epsilon *= 0
                    force.setParticleParameters(index, charge, sigma, epsilon)
                for index in range(force.getNumExceptions()):
                    p1, p2, chargeProd, sigma, epsilon = force.getExceptionParameters(index)
                    if not self.exception_charges:
                        chargeProd *= 0
                    if not self.exception_epsilons:
                        epsilon *= 0
                    force.setExceptionParameters(index, p1, p2, chargeProd, sigma, epsilon)
            elif force.__class__.__name__ == "PeriodicTorsionForce":
                for index in range(force.getNumTorsions()):
                    p1, p2, p3, p4, periodicity, phase, K = force.getTorsionParameters(index)
                    if not self.torsions:
                        K *= 0
                    force.setTorsionParameters(index, p1, p2, p3, p4, periodicity, phase, K)

    def create_system(self, topology, molecules=None):
        """
        Create a system from the specified topology.

        .. todo :: Add support for OpenFF Topology objects once they can be converted to OpenMM Topology objects.

        Parameters
        ----------
        topology : openmm.app.Topology object
            The topology describing the system to be created
        molecules : openff.toolkit.topology.Molecule or list of Molecules, optional, default=None
            Can alternatively be an object (such as an OpenEye OEMol or RDKit Mol or SMILES string) that can be used
            to construct a Molecule.
            Can also be a list of Molecule objects or objects that can be used to construct a Molecule.
            If specified, these molecules will be recognized and parameterized as needed.
            The parameters will be cached in case they are encountered again the future.

        Returns
        -------
        system : openmm.System
            A system object generated from the topology
        """

        # Inform the template generator about any specified molecules
        if (self.template_generator is not None) and (molecules is not None):
            self.template_generator.add_molecules(molecules)

        # Build the kwargs to use
        import copy

        forcefield_kwargs = copy.deepcopy(self.forcefield_kwargs)
        if topology.getPeriodicBoxVectors() is None:
            forcefield_kwargs.update(self.nonperiodic_forcefield_kwargs)
        else:
            forcefield_kwargs.update(self.periodic_forcefield_kwargs)

        # Build the System
        system = self.forcefield.createSystem(topology, **forcefield_kwargs)

        # Modify other forces as requested
        self._modify_forces(system)

        # Post-process the System if requested
        if self.postprocess_system is not None:
            system = self.postprocess_system(system)

        return system
