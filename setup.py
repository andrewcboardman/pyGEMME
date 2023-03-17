from setuptools import find_packages
from setuptools import setup
import os

def package_files(directory,ext):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            if filename.split('.')[-1] == ext:
                paths.append(os.path.join('/'.join(path.split('/')[1:]), filename))
    return paths

jet_class_files = package_files('pyGEMME/jet','class')
jet_jar_files = package_files('pyGEMME/jet','jar')
print(jet_class_files)
#jet_class_files = package_files('pyGEMME/jet','class')



setup(
    name='pyGEMME',
    version='0.1.0',
    install_requires=['certifi==2022.12.7', 'pandas', 'itertools-resources','scipy', 'biopython','tqdm','matplotlib'],
    packages=find_packages(),
    package_dir={},
    url='https://github.com/andrewcboardman/pyGEMME',
    license='MIT',
    author='Andrew Boardman',
    author_email='andrewcboardman6@gmail.com',
    description='Your main project',
    include_package_data=True,
    package_data={'': [
        'data/alphabets/*.txt',
        'data/jet_matrices/*',
        'data/matrix/blosum62p.txt',
        'data/*',
        'tools/muscle', 
        *jet_class_files,
        *jet_jar_files,
        'pred/*.R']}
)