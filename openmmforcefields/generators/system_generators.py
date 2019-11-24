"""
System generators that build an OpenMM System object from a Topology object.

"""

################################################################################
# LOGGER
################################################################################

import logging
_logger = logging.getLogger("perses.forcefields.system_generators")

################################################################################
# System generator base class
################################################################################

class SystemGenerator(object):
    """
    A utility class to generate OpenMM Systems from Open Force Field Topology objects
    that may contain both biopolymers and small molecules.

    """
    def __init__(self, particle_charges=True, exception_charges=True, particle_epsilons=True, exception_epsilons=True, torsions=True):
        self._particle_charges = particle_charges
        self._exception_charges = exception_charges
        self._particle_epsilons = particle_epsilons
        self._exception_epsilons = exception_epsilons
        self._torsions = torsions

    def postprocess_system(self):

        # Add barostat if requested.
        if self._barostat is not None:
            import numpy as np
            from simtk import openmm
            MAXINT = np.iinfo(np.int32).max
            barostat = openmm.MonteCarloBarostat(*self._barostat)
            seed = np.random.randint(MAXINT)
            barostat.setRandomNumberSeed(seed)
            system.addForce(barostat)

        # TODO: Remove CM Motion?

        for force in system.getForces():
            if force.__class__.__name__ == 'NonbondedForce':
                for index in range(force.getNumParticles()):
                    charge, sigma, epsilon = force.getParticleParameters(index)
                    if not self._particle_charges:
                        charge *= 0
                    if not self._particle_epsilons:
                        epsilon *= 0
                    force.setParticleParameters(index, charge, sigma, epsilon)
                for index in range(force.getNumExceptions()):
                    p1, p2, chargeProd, sigma, epsilon = force.getExceptionParameters(index)
                    if not self._exception_charges:
                        chargeProd *= 0
                    if not self._exception_epsilons:
                        epsilon *= 0
                    force.setExceptionParameters(index, p1, p2, chargeProd, sigma, epsilon)
            elif force.__class__.__name__ == 'PeriodicTorsionForce':
                for index in range(force.getNumTorsions()):
                    p1, p2, p3, p4, periodicity, phase, K = force.getTorsionParameters(index)
                    if not self._torsions:
                        K *= 0
                    force.setTorsionParameters(index, p1, p2, p3, p4, periodicity, phase, K)

    def create_system(self, topology):
        """
        Create a system from the specified topology.

        Parameters
        ----------
        topology : openmmtools.topology.Topology object
            The topology describing the system to be created

        Returns
        -------
        system : simtk.openmm.System
            A system object generated from the topology
        """
        raise NotImplementedError('Tried to call create_system() for abstract base class')

################################################################################
# GAFF system generator
################################################################################

class GAFFSystemGenerator(SystemGenerator):
    """
    This is a utility class to generate OpenMM Systems from Open Force Field Topology objects using AMBER
    protein force fields and GAFF small molecule force fields.

    .. warning :: This API is experimental and subject to change.

    >>> # Initialize a SystemGenerator
    >>> from openmmforcefields.generators import GAFFSystemGenerator
    >>> system_generator = GAFFSystemGenerator(forcefields=['amber/protein.ff14SB.xml', 'amber/tip3p.xml'], gaff_version='2.11', cache='gaf-2.11.json')
    >>> # Create an OpenMM System from an Open Force Field toolkit Topology object
    >>> system = system_generator.create_system(topology)

    """

    def __init__(self, forcefields=None, gaff_version='2.11', cache=None, **kwargs):
        """
        Create a new SystemGenerator for Amber systems that contain GAFF-parameterized small molecules.

        Parameters
        ----------
        forcefields : list of str, optional, default=None
            List of the names of ffxml files that will be used in System creation for the biopolymer.
        gaff_version : str, optional, default='2.11'
            GAFF version to use. One of ('1.4', '1.8', '1.81', '2.1', '2.11')
        cache : filename or TinyDB instance
            JSON filename or TinyDB instance that can be used to cache parameterized small molecules by OEGAFFTemplateGenerator
        """
        # Call base class
        super(GAFFSystemGenerator).__init__(**kwargs)

        # Cache force fields and settings to use
        self._forcefield_xmls = forcefields if (forcefields is not None) else list()
        self._forcefield_kwargs = forcefield_kwargs if forcefield_kwargs is not None else {}

        # Create and cache a ForceField object
        from simtk.openmm import app
        self._forcefield = app.ForceField(*self._forcefield_xmls)

        # Create and cache a residue template generator
        from .gaff import GAFFTemplateGenerator
        self._generator = GAFFTemplateGenerator(molecules=molecules, cache=cache)
        # TODO: Do we need to worry about parameter collisions with frcmod-generated additional parameters?
        self._forcefield.registerTemplateGenerator(self._generator.generator)

        # Ensure that center-of-mass motion removal is not added
        if 'removeCMMotion' not in self._forcefield_kwargs:
            self._forcefield_kwargs['removeCMMotion'] = False

        # Cache barostat if needed
        self._barostat = None
        if barostat is not None:
            pressure = barostat.getDefaultPressure()
            if hasattr(barostat, 'getDefaultTemperature'):
                temperature = barostat.getDefaultTemperature()
            else:
                temperature = barostat.getTemperature()
            frequency = barostat.getFrequency()
            self._barostat = (pressure, temperature, frequency)

    def get_forcefield(self):
        """
        Return the associated ForceField object.

        Returns
        -------
        forcefield : simtk.openmm.app.ForceField
            The current ForceField object.
        """
        return self._forcefield

    def create_system(self, topology):
        """
        Create a system from the specified topology.

        Parameters
        ----------
        topology : openmmtools.topology.Topology object
            The topology describing the system to be created

        Returns
        -------
        system : simtk.openmm.System
            A system object generated from the topology
        """
        try:
            system = self._forcefield.createSystem(new_topology, **self._forcefield_kwargs)
        except Exception as e:
            # TODO: Write topology that failed to parameterize to a PDB file
            from simtk import unit
            import numpy as np
            nparticles = sum([1 for atom in new_topology.atoms()])
            positions = unit.Quantity(np.zeros([nparticles,3], np.float32), unit.angstroms)
            # Write PDB file of failed topology
            from simtk.openmm.app import PDBFile
            outfile = open('BuildSystem-failure.pdb', 'w')
            pdbfile = PDBFile.writeFile(topology.to_openmm(), positions, outfile)
            outfile.close()
            msg = str(e)
            import traceback
            msg += traceback.format_exc()
            msg += "\n"
            msg += "PDB file written as 'BuildSystem-failure.pdb'"
            raise Exception(msg)

        # Postprocess system
        self.postprocess_system(system)

        return system

    @property
    def ffxmls(self):
        return self._forcefield_xmls

    @property
    def forcefield(self):
        return self._forcefield

    @property
    def generator(self):
        return self._generator

################################################################################
# SMIRNOFF system generator
################################################################################

class SMIRNOFFSystemGenerator(SystemGenerator):
    """
    This is a utility class to generate OpenMM Systems from Open Force Field Topology objects using AMBER
    biopolymer force fields and SMIRNOF small molecule force fields.

    .. warning :: This API is experimental and subject to change.

    >>> # Initialize a SystemGenerator
    >>> from openmmforcefields.generators import GAFFSystemGenerator
    >>> system_generator = GAFFSystemGenerator(forcefields=['amber/protein.ff14SB.xml', 'amber/tip3p.xml'], gaff_version='2.11', cache='gaf-2.11.json')
    >>> # Create an OpenMM System from an Open Force Field toolkit Topology object
    >>> system = system_generator.create_system(topology)

    """

    def __init__(self, forcefields=None, smirnoff_forcefield=None, metadata=None, cache=None, **kwargs):
        """
        Create a new SystemGenerator for AMBER systems that contain SMIRNOFF-parameterized small molecules.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        forcefields : list of str, optional, default=None
            List of the names of ffxml files that will be used in System creation for the biopolymer.
        sminroff_forcefield : str, optional, default=None
            SMIRNOFF force field filename
        cache : filename or TinyDB instance
            JSON filename or TinyDB instance that can be used to cache parameterized small molecules
        """
        # Call base class
        super(SMIRNOFFSystemGenerator).__init__(**kwargs)

        # TODO: Create openforcefield ForceField object
        self._openforcefield_forcefield = openforcefield.typing.engines.smirnoff.ForceField(smirnoff_forcefield)

        # Cache force fields and settings to use
        self._forcefield_xmls = forcefields_to_use
        self._forcefield_kwargs = forcefield_kwargs if forcefield_kwargs is not None else {}

        # Create and cache a ForceField object
        from simtk.openmm import app
        self._forcefield = app.ForceField(*self._forcefield_xmls)

        # Create and cache a residue template generator
        from .gaff import GAFFTemplateGenerator
        self._generator = GAFFTemplateGenerator(molecules=molecules, cache=cache)
        # TODO: Do we need to worry about parameter collisions with frcmod-generated additional parameters?
        self._forcefield.registerTemplateGenerator(self._generator.generator)

        # Ensure that center-of-mass motion removal is not added
        if 'removeCMMotion' not in self._forcefield_kwargs:
            self._forcefield_kwargs['removeCMMotion'] = False

        # Cache barostat if needed
        self._barostat = None
        if barostat is not None:
            pressure = barostat.getDefaultPressure()
            if hasattr(barostat, 'getDefaultTemperature'):
                temperature = barostat.getDefaultTemperature()
            else:
                temperature = barostat.getTemperature()
            frequency = barostat.getFrequency()
            self._barostat = (pressure, temperature, frequency)

    def create_system(self, topology):
        """
        Create a system from the specified topology.

        Parameters
        ----------
        topology : openmmtools.topology.Topology object
            The topology describing the system to be created

        Returns
        -------
        system : simtk.openmm.System
            A system object generated from the topology
        """
        # Create the receptor system
        receptor_system = self._openmm_forcefield.createSystem(new_topology, **self._forcefield_kwargs)

        # Parameterize any small molecules
        smirnoff_system = self._smirnoff_forcefield.create_system(small_molecule)

        # TODO: Combine with ParmEd
        system = None

        # Postprocess system
        self.postprocess_system(system)

        return system

################################################################################
# Dummy system generator
################################################################################

class DummyForceField(object):
    """
    Dummy force field that can add basic parameters to any system for testing purposes.

    * All particles are assigned carbon mass
    * All particles interact with a repulsive potential
    * All bonds have equilibrium length 1 A
    * All angles have equilibrium angle dependent on number of substituents of central atom
        * 2, 3 bonds: 120 degrees
        * 4 bonds: 109.8 degrees
        * 5 or more bonds: 90 degrees
    * Torsions are added with periodicity 3, but no barrier height

    """
    def create_system(self, topology, **kwargs):
        """
        Create a System object with simple parameters from the provided Topology

        Any kwargs are ignored.

        Parameters
        ----------
        topology : openforcefield.topology.Topology
            The Topology to be parameterized

        Returns
        -------
        system : simtk.openmm.System
            The System object

        """
        # TODO: Allow periodicity to be determined from topology

        from openmmtools.constants import kB
        kT = kB * 300*unit.kelvin # hard-coded temperature for setting energy scales

        # Create a System
        system = openmm.System()

        # Add particles
        mass = 12.0 * unit.amu
        for atom in topology.atoms:
            system.addParticle(mass)

        # Add simple repulsive interactions
        # TODO: Use softcore repulsive interaction; Gaussian times switch?
        nonbonded = openmm.CustomNonbondedForce('100/(r/0.1)^4')
        nonbonded.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic);
        nonbonded.setCutoffDistance(1*unit.nanometer)
        system.addForce(nonbonded)
        for atom in topology.atoms:
            nonbonded.addParticle([])

        # Build a list of which atom indices are bonded to each atom
        bondedToAtom = []
        for atom in topology.atoms():
            bondedToAtom.append(set())
        for (atom1, atom2) in topology.bonds():
            bondedToAtom[atom1.index].add(atom2.index)
            bondedToAtom[atom2.index].add(atom1.index)
        return bondedToAtom

        # Add bonds
        bond_force = openmm.HarmonicBondForce()
        r0 = 1.0 * unit.angstroms
        sigma_r = 0.1 * unit.angstroms
        Kr = kT / sigma_r**2
        for atom1, atom2 in topology.bonds():
            bond_force.addBond(atom1.index, atom2.index, r0, Kr)
        system.addForce(bond_force)

        # Add angles
        uniqueAngles = set()
        for bond in topology.bonds():
            for atom in bondedToAtom[bond.atom1]:
                if atom != bond.atom2:
                    if atom < bond.atom2:
                        uniqueAngles.add((atom, bond.atom1, bond.atom2))
                    else:
                        uniqueAngles.add((bond.atom2, bond.atom1, atom))
            for atom in bondedToAtom[bond.atom2]:
                if atom != bond.atom1:
                    if atom > bond.atom1:
                        uniqueAngles.add((bond.atom1, bond.atom2, atom))
                    else:
                        uniqueAngles.add((atom, bond.atom2, bond.atom1))
        angles = sorted(list(uniqueAngles))
        theta0 = 109.5 * unit.degrees # TODO: Adapt based on number of bonds to each atom?
        sigma_theta = 10 * unit.degrees
        Ktheta = kT / sigma_theta**2
        angle_force = openmm.HarmonicAngleForce()
        for (atom1, atom2, atom3) in angles:
            angles.addAngle(atom1.index, atom2.index, atom3.index, theta0, Ktheta)
        system.addForce(angle_force)

        # Make a list of all unique proper torsions
        uniquePropers = set()
        for angle in angles:
            for atom in bondedToAtom[angle[0]]:
                if atom not in angle:
                    if atom < angle[2]:
                        uniquePropers.add((atom, angle[0], angle[1], angle[2]))
                    else:
                        uniquePropers.add((angle[2], angle[1], angle[0], atom))
            for atom in bondedToAtom[angle[2]]:
                if atom not in angle:
                    if atom > angle[0]:
                        uniquePropers.add((angle[0], angle[1], angle[2], atom))
                    else:
                        uniquePropers.add((atom, angle[2], angle[1], angle[0]))
        propers = sorted(list(uniquePropers))
        torsion_force = openmm.PeriodicTorsionForce()
        periodicity = 3
        phase = 0.0 * unit.degrees
        Kphi = 0.0 * kT
        for (atom1, atom2, atom3, atom4) in propers:
            torsion_force.add_torsion(atom1.index, atom2.index, atom3.index, atom4.index, periodicity, phase, Kphi)
        system.addForce(torsion_force)

        return system

class DummySystemGenerator(SystemGenerator):
    """
    Dummy SystemGenerator that employs a universal simple force field.

    """
    def __init__(self, **kwargs):
        """
        Create a DummySystemGenerator with universal simple force field.

        All parameters except 'barostat' are ignored.

        """
        # Call base class
        super(DummySystemGenerator).__init__(**kwargs)

        self._forcefield = DummyForceField()
        self._forcefield_xmls = list()
        self._forcefield_kwargs = dict()
        self._barostat = None
        if barostat is not None:
            pressure = barostat.getDefaultPressure()
            if hasattr(barostat, 'getDefaultTemperature'):
                temperature = barostat.getDefaultTemperature()
            else:
                temperature = barostat.getTemperature()
            frequency = barostat.getFrequency()
            self._barostat = (pressure, temperature, frequency)

    def create_system(self, topology):
        """
        Create a system from the specified topology.

        Parameters
        ----------
        topology : openmmtools.topology.Topology object
            The topology describing the system to be created

        Returns
        -------
        system : simtk.openmm.System
            A system object generated from the topology
        """
        try:
            system = self._forcefield.createSystem(new_topology, **self._forcefield_kwargs)
        except Exception as e:
            # TODO: Write topology that failed to parameterize to a PDB file
            from simtk import unit
            import numpy as np
            nparticles = sum([1 for atom in new_topology.atoms()])
            positions = unit.Quantity(np.zeros([nparticles,3], np.float32), unit.angstroms)
            # Write PDB file of failed topology
            from simtk.openmm.app import PDBFile
            outfile = open('BuildSystem-failure.pdb', 'w')
            pdbfile = PDBFile.writeFile(topology.to_openmm(), positions, outfile)
            outfile.close()
            msg = str(e)
            import traceback
            msg += traceback.format_exc()
            msg += "\n"
            msg += "PDB file written as 'BuildSystem-failure.pdb'"
            raise Exception(msg)

        # Postprocess system
        self.postprocess_system(system)

        return system
