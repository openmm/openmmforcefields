# Add imports here
# Handle versioneer
from openmmforcefields.utils import get_ffxml_path

from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
