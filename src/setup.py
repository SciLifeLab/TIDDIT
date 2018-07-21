from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "TIDDIT",
    ext_modules = cythonize(['TIDDIT_clustering.py',"DBSCAN.py"]),
)
