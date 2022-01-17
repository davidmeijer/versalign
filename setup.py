from setuptools import setup, find_packages

from monomer_aligner.version import __version__

setup(
    name='monomer_aligner',
    version=__version__,
    packages=find_packages(),
    author='David Meijer',
    author_email='david.meijer@wur.nl'
)
