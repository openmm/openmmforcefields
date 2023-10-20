# Add imports here
# Handle versioneer
from ._version import get_versions
from .utils import get_ffxml_path

versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

from . import _version
__version__ = _version.get_versions()['version']
