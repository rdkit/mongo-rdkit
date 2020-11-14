import io
import os
import re

from setuptools import find_packages, setup

# Adapted from https://stackoverflow.com/a/39671214
this_directory = os.path.dirname(os.path.realpath(__file__))
version_matches = re.search(r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                            io.open(f'{this_directory}/mongordkit/__init__.py', encoding='utf_8_sig').read())
if version_matches is None:
    raise Exception('Could not determine mongordkit version from __init__.py')
__version__ = version_matches.group(1)

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name='mongordkit',
    version=__version__,
    author='chriswzou',
    author_email='chriswzou@example.com',
    description='Mongo-rdkit is an integration between MongoDB and RDKit',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/rdkit/mongo-rdkit',

    packages=find_packages(),

    install_requires=[
      'attrs==19.3.0',
      'importlib-metadata==1.6.1',
      'mongomock==3.19.0',
      'more-itertools==8.3.0',
      'packaging==20.4',
      'pluggy==0.13.1',
      'py==1.8.1',
      'pymongo==3.10.1',
      'pyparsing==2.4.7',
      'pytest==5.4.3',
      'sentinels==1.0.0',
      'wcwidth==0.2.4',
      'zipp==3.1.0',
    ],
    include_package_data=True,
    zip_safe=True,

    classifiers=(
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ),
)
