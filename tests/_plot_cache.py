import tempfile
from pathlib import Path


_CACHE_DIRECTORY = tempfile.TemporaryDirectory(prefix='gritic-plot-tests-')
PLOT_CACHE_ROOT = Path(_CACHE_DIRECTORY.name)
MPL_CONFIG_DIRECTORY = PLOT_CACHE_ROOT / 'matplotlib'
XDG_CACHE_DIRECTORY = PLOT_CACHE_ROOT / 'cache'

