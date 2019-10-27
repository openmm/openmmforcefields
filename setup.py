"""
openmmforcefields
Additional force fields for OpenMM.
"""
import sys, os
from setuptools import setup, find_packages
import versioneer

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])

try:
    # Symlink converted ffxml directories so we don't need to copy files
    os.makedirs('openmmforcefields/ffxml')
    os.symlink('../../amber/ffxml/', 'openmmforcefields/ffxml/amber')
    os.symlink('../../charmm/ffxml/', 'openmmforcefields/ffxml/charmm')

    setup(
        # Self-descriptive entries which should always be present
        name='openmmforcefields',
        author='Chodera lab @ MSKCC',
        author_email='john.chodera@choderalab.org',
        description=short_description[0],
        long_description=long_description,
        long_description_content_type="text/markdown",
        keywords = "molecular mechanics, forcefield, OpenMM, AMBER, CHARMM, GAFF",
        url = "http://github.com/choderalab/openmm-forcefields",
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        license='MIT',

        classifiers=[
            "Development Status :: 3 - Alpha",
            "Topic :: Utilities",
            "License :: OSI Approved :: MIT",
        ],

        # Don't install as an egg, since OpenMM can't find ffxml files if we do
        zip_safe=False,

        # Which Python importable modules should be included when your package is installed
        # Handled automatically by setuptools. Use 'exclude' to prevent some specific
        # subpackage(s) from being added, if needed
        packages=find_packages(),

        # Optional include package data to ship with your package
        # Customize MANIFEST.in if the general case does not suit your needs
        # Comment out this line to prevent the files from being packaged with your software
        include_package_data=True,

        # Allows `setup.py test` to work correctly with pytest
        setup_requires=[] + pytest_runner,

        entry_points={
            'console_scripts': [],
            'openmm.forcefielddir' : [
                'ffxml = openmmforcefields.utils:get_ffxml_path',
                ],
            },

        package_data={
            'openmmforcefields': ['ffxml/amber/*.xml', 'ffxml/charmm/*.xml']
            },
    )

finally:
    # Clean up temporary symlinks
    os.unlink('openmmforcefields/ffxml/amber')
    os.unlink('openmmforcefields/ffxml/charmm')
    os.removedirs('openmmforcefields/ffxml')
