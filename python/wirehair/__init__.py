"""Canonical Python import for the Wirehair FEC binding.

The historical :mod:`whirehair` import remains supported and exposes the same
objects.  This package intentionally contains no native shared library; see
the distribution README for the system-library discovery contract.
"""

from whirehair import *  # noqa: F401,F403
from whirehair import __all__ as __all__
