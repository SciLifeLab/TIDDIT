from setuptools import setup
import numpy

try:
    from Cython.Build import cythonize
    has_cython = True
except ImportError:
    has_cython = False

if has_cython:
    ext_modules = cythonize([
        "tiddit_signal.pyx",
        "tiddit_coverage.pyx",
        "tiddit_coverage_analysis.pyx",
        "tiddit_variant.pyx",
        "tiddit_contig_analysis.pyx"])
else:
    ext_modules = []

setup(
    name = 'tiddit',
    version = '4.0.0',
    url = "https://github.com/J35P312/SVDB",
    author = "Jesper Eisfeldt",
    author_email= "jesper.eisfeldt@scilifelab.se",
    ext_modules = ext_modules,
    include_dirs=[numpy.get_include()],
    packages = ['tiddit'],
    install_requires = ['numpy'],
)
#    entry_points = {'console_scripts': ['svdb = svdb.__main__:main']},
