# Python imports
import os
import subprocess
import shutil
import glob
import tempfile

# local imports
# need relative imports!
from .. import FLASH_SRC_DIR
from . import setup_parse
from . import setup_globals
from .setup_globals import gvars

TEMPDIR = tempfile.gettempdir()

def clean_main(lvl):
    """Cleans up the current run directory to a given level."""
    # Setup gvars
    gvars.init(FLASH_SRC_DIR)

    # Clean up runs & tempfiles
    if 1 <= lvl:
        run_dirs = glob.glob(gvars.run_dir_prefix + "*")
        vc_types = ['git', 'hg', 'svn', 'release']
        temp_dirs = [d for vc in vc_types for d in glob.glob(os.path.join(TEMPDIR, "*"+vc+"*"))]
        for d in run_dirs + temp_dirs:
            shutil.rmtree(d, ignore_errors=True)

    # Clean up build
    if 2 <= lvl:
        shutil.rmtree(gvars.project_build_dir, ignore_errors=True)

    # Clean up local src dir
    if 3 <= lvl:
        shutil.rmtree(gvars.project_setup_dir, ignore_errors=True)
        if os.path.exists(gvars.desc_filename):
            os.remove(gvars.desc_filename)


def main(ns, rc):
    """Removes previously created files."""
    clean_main(ns.level)


USAGE = """usage: flmake clean [<level>]

Cleans the current directory by removing the files
and directories generated by previous runs of flmake.
The <level> dictates the level at which to clean up to.
If no level is specified, a default level of 1 is taken.
Levels are defined as:

  1: run dirs
  2: run dirs, the object dirs & files
  3: run dirs, the object dirs & files, the src dir

Other levels may be defined in future version."""
