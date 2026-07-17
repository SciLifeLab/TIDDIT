from setuptools import setup
import numpy
import pyximport
import pysam
pyximport.install()


try:
    from Cython.Build import cythonize
    has_cython = True
except ImportError:
    has_cython = False

from setuptools import Extension

if has_cython:
    extensions = [
        Extension("tiddit.tiddit_signal", ["tiddit/tiddit_signal.pyx"]),
        Extension("tiddit.tiddit_coverage", ["tiddit/tiddit_coverage.pyx"]),
        Extension("tiddit.tiddit_cluster", ["tiddit/tiddit_cluster.pyx"]),
        Extension("tiddit.tiddit_coverage_analysis", ["tiddit/tiddit_coverage_analysis.pyx"]),
        Extension("tiddit.tiddit_gc", ["tiddit/tiddit_gc.pyx"]),
        Extension("tiddit.silverfish", ["tiddit/silverfish.pyx"]),
        Extension("tiddit.graphlib", ["tiddit/graphlib.pyx"]),
        Extension("tiddit.tiddit_variant", ["tiddit/tiddit_variant.pyx"]),
        Extension("tiddit.tiddit_contig_analysis", ["tiddit/tiddit_contig_analysis.pyx"]),
    ]
    ext_modules = cythonize(extensions)
else:
    ext_modules = []


setup(
    name = 'tiddit',
    version = '3.9.6',


    url = "https://github.com/SciLifeLab/TIDDIT",
    author = "Jesper Eisfeldt",
    author_email= "jesper.eisfeldt@scilifelab.se",
    ext_modules = ext_modules,
    extra_link_args=pysam.get_libraries(),
    define_macros=pysam.get_defines(),
    include_dirs=[numpy.get_include()]+pysam.get_include(),
    packages = ['tiddit'],
    install_requires = ['numpy','pysam','joblib'],
    entry_points = {'console_scripts': ['tiddit = tiddit.__main__:main']},

)
