from setuptools import setup, find_packages

setup(
    name = "pyGEMME",
    version = "0.1",
    install_requires = ['numpy', 'scipy', 'matplotlib', 'pandas',  'seaborn', 'biopython', 'requests','tqdm'],
    scripts= [
            'scripts/MMseqs_query.py',
            'scripts/run_pyGEMME.py'
            ],
    packages = ['pyGEMME'],
    url = "https://github.com/andrewcboardman/pyGEMME.git",
    license = "MIT",
    author = "Andrew Boardman",
    author_email = "acb95@cam.ac.uk",
    description = "A python reimplementation of the GEMME algorithm (copyright Elodie Laine, 2019)"
)