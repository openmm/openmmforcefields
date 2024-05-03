# Add imports here
# Handle versioneer
from openmmforcefields.utils import get_ffxml_path

from . import _version

__version__ = _version.get_versions()["version"]
