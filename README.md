# mongo-rdkit
[![Build Status](https://dev.azure.com/cwzou/mongo-rdkit/_apis/build/status/rdkit.mongo-rdkit?branchName=master)](https://dev.azure.com/cwzou/mongo-rdkit/_build/latest?definitionId=1&branchName=master)

<a href="https://github.com/rdkit/mongo-rdkit">Mongo-rdkit</a> is an integration between MongoDB, 
a NoSQL database platform, and RDKit, a collection of chemoinformatics and machine-learning software.
This package contains tools to create and manipulate a chemically-intelligent database, as well as
methods for high-performance searches on the database that leverage native MongoDB features.

Useful links: 
* [BSD License](https://github.com/rdkit/mongo-rdkit/blob/master/LICENSE) - a business friendly license for open-source.
* [Jupyter Notebooks](https://github.com/rdkit/mongo-rdkit/tree/master/docs) - resources for getting started.

## Documentation
Jupyter Notebooks and resources for getting started in the [docs](https://github.com/rdkit/mongo-rdkit/tree/master/docs) 
folder on GitHub

## Installation
As the package is not officially configured with a setup.py file or pushed onto PyPi, these are working install instructions.

### macOS and Linux:
Ensure that you have either Anaconda or Miniconda installed and that `conda` has been added to `PATH`.

Clone the repository into your desired directory.

Navigate so that your current working directory is `mongo-rdkit`.

Create a conda environment called mongo_rdkit that includes all dependencies needed for this package: 
```
conda env create --quiet --force --file env.yml
```
Activate said conda environment:
```
source activate mongo_rdkit
```
Add the package mongo-rdkit to `PYTHONPATH`:
```
export PYTHONPATH="$PWD"
```
Check that `your/path/here/mongo-rdkit` lies in `PYTHONPATH`:
```
echo $PYTHONPATH
```
You can now `import mongordkit` in your Python interpreter or run all tests using the `pytest` command.

### Windows:




## License
Code released under the BSD License.
