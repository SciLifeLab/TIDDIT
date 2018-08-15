from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "TIDDIT",
    ext_modules = cythonize(['TIDDIT_calling.py',"TIDDIT_coverage.py","TIDDIT_filtering.py","TIDDIT_signals.py","DBSCAN.py"]),
)
