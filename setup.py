import os
from os.path import relpath, join
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def find_package_data(data_root, package_root):
    files = []
    for root, dirnames, filenames in os.walk(data_root):
        for fn in filenames:
            files.append(relpath(join(root, fn), package_root))
    return files

setup(
    name = "openmm-forcefields",
    version = "0.0.0",
    author = "Rafal P. Wiewiora, Chaya Stern, Peter Eastman, and John D. Chodera",
    author_email = "john.chodera@choderalab.org",
    description = ("Biomolecular forcefields and small molecule support for OpenMM"),
    license = "MIT",
    keywords = "molecular mechanics, forcefield, OpenMM, AMBER, CHARMM, GAFF",
    url = "http://github.com/choderalab/openmm-forcefields",
    packages=[
        'openmmforcefields',
        ],
    long_description=read('README.md'),
    install_requires=[
        'openmm>=7.2',
    ]
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT",
    ],
    entry_points={
        'console_scripts': [],
        'openmm.forcefielddir' : [
            'amber = openmmforcefields.utils.get_amber_directory',
            'charmm = openmmforcefields.utils.get_charmm_directory'
            ],
        },
    package_data={'openmmforcefields': ['amber/ffxml/*.xml', 'charmm/ffxml/*.xml']},
)
