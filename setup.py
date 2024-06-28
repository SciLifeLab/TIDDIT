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

if has_cython:
    ext_modules = cythonize([
        "tiddit/tiddit_signal.pyx",
        "tiddit/tiddit_coverage.pyx",
        "tiddit/tiddit_cluster.pyx",
        "tiddit/tiddit_coverage_analysis.pyx",
        "tiddit/tiddit_variant.pyx",
        "tiddit/tiddit_contig_analysis.pyx"])
else:
    ext_modules = []

setup(
    name = 'tiddit',
    version = '3.7.0',


    url = "https://github.com/SciLifeLab/TIDDIT",
    author = "Jesper Eisfeldt",
    author_email= "jesper.eisfeldt@scilifelab.se",
    ext_modules = ext_modules,
    extra_link_args=pysam.get_libraries(),
    define_macros=pysam.get_defines(),
    include_dirs=[numpy.get_include()]+pysam.get_include(),
    packages = ['tiddit'],
    install_requires = ['numpy','pysam'],
    entry_points = {'console_scripts': ['tiddit = tiddit.__main__:main']},

)
