from setuptools import setup, find_packages

setup(
    name = "pyGEMME",
    version = "0.1",
    install_requires = ['numpy', 'scipy', 'matplotlib', 'pandas',  'seaborn', 'biopython'],
    packages = find_packages('src'),
    package_dir = {'':'src'},
    url = "https://github.com/andrewcboardman/pyGEMME.git",
    license = "MIT",
    author = "Andrew Boardman",
    author_email = "acb95@cam.ac.uk",
    description = "A python reimplementation of the GEMME algorithm (Elodie Laine, 2019)"
)