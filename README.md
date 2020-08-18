# mongo-rdkit
[![Build Status](https://dev.azure.com/cwzou/mongo-rdkit/_apis/build/status/rdkit.mongo-rdkit?branchName=master)](https://dev.azure.com/cwzou/mongo-rdkit/_build/latest?definitionId=1&branchName=master)

<a href="https://github.com/rdkit/mongo-rdkit">Mongo-rdkit</a> is an integration between MongoDB, 
a NoSQL database platform, and RDKit, a collection of cheminformatics and machine-learning software.
This package contains tools to create and manipulate a chemically-intelligent database, as well as
methods for high-performance searches on the database that leverage native MongoDB features.

Useful links: 
* [BSD License](https://github.com/rdkit/mongo-rdkit/blob/master/LICENSE) - a business friendly license for open-source.
* [Jupyter Notebooks](https://github.com/rdkit/mongo-rdkit/tree/master/docs) - walkthroughs for main functionality.
* [Testing Guide](https://github.com/rdkit/mongo-rdkit/blob/master/docs/testing.md) - walkthrough of running `mongordkit` tests.

## Documentation
Jupyter Notebooks and resources for getting started in the [docs](https://github.com/rdkit/mongo-rdkit/tree/master/docs) 
folder on GitHub.

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
Similarly, ensure that `conda` has been added to `PATH`.

Clone the repository into your desired directory and navigate into it.

Create a conda environment called mongo_rdkit that includes dependencies:
```
conda env create --quiet --force --file env.yml
```
Activate this conda environment: 
```
call activate mongo_rdkit
```
Check that you are able to import mongordkit:
```
python -c "import mongordkit"
```
If this fails, you may need to add the current directory manually to `PYTHONPATH`:
```
set PYTHONPATH=%PYTHONPATH%;C:.
```
You can now use `mongordkit` in your interpreter and run tests using `python -m pytest`.
## Package Contents
### Modules
`mongordkit` contains two main modules, each of which contains a variety of importable methods and classes.
`Database` contains functionality for writing and registering data. `Search` contains functionality for setting up and performing
substructure and similarity search. Detailed walkthroughs can be found in the notebooks, listed below.

### Notebooks
- **Creating and Writing to MongoDB**: documentation and demos for creating and modifying mongo-rdkit databases.
- **Similarity and Substructure Search**: documentation and demos for similarity and substructure search.
- **Similarity Benchmarking**: documentation for reproducing similarity benchmarking.
- **Substructure Benchmarking**: documentation for reproducing substructure benchmarking.

### Configuration
- **azure_pipelines.yml**: CI/CD pipeline configurations.
- **conftest.py**: `pytest` configurations.
- **env.yml**: required dependencies.

## License
Code released under the BSD License.
