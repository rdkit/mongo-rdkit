from bson import Binary
from mongomock import ObjectId
import sys
from pathlib import Path
import mongordkit
import pymongo
import rdkit
from rdkit import Chem
import mongomock
import random
from rdkit.Chem import AllChem
from mongordkit.Database import write
from mongordkit.Search import similarity
from mongordkit.Search.tests import utils

sys.path.append(Path('.').resolve().parent.parent.parent)


def test_pythonSetup():
    assert 200 == len(utils.setupPythonDB('data/test_data/first_200.props.sdf'))


def test_zeroThreshold():
    """
    If we search the database with a Tanimoto threshold of zero, we should get
    the entire contents of the database back.
    """
    db_python = utils.setupPythonDB('data/test_data/first_200.props.sdf')
    db_mongo = utils.setupMongoDB()
    write.writeFromSDF(db_mongo, 'data/test_data/first_200.props.sdf', 'test')
    similarity.addMorganFingerprints(db_mongo)

    mol = Chem.Mol(db_python[0]['rdmol'])
    assert 200 == len(utils.similaritySearchPython(mol, db_python, 0))
    assert 200 == len(similarity.similaritySearchAggregate(mol, db_mongo, 0))
    assert 200 == len(similarity.similaritySearch(mol, db_mongo, 0))


def test_similarityAccuracy():
    """
    Tests for basic accuracy against a brute-force constructed Python 'database'
    at thresholds 0.2, 0.4, 0.6, 0.8, and 1. This test is relatively long and
    will modify your local MongoDB instance.
    """
    db_python = utils.setupPythonDB('data/test_data/first_200.props.sdf')
    db_mongo = utils.setupMongoDB()
    write.writeFromSDF(db_mongo, 'data/test_data/first_200.props.sdf', 'test')
    similarity.addMorganFingerprints(db_mongo)
    thresholds = [0.2, 0.4, 0.6, 0.8, 1]
    for t in thresholds:
        for i in range(200):
            mol = Chem.Mol(db_python[i]['rdmol'])
            search_python = utils.similaritySearchPython(mol, db_python, t)
            search_mongo_aggregate = similarity.similaritySearchAggregate(mol, db_mongo, t)
            search_mongo = similarity.similaritySearch(mol, db_mongo, t)
            assert sorted(search_python) == sorted(search_mongo) == sorted(search_mongo_aggregate)


def test_similarityProgression():
    db_python = utils.setupPythonDB('data/test_data/first_200.props.sdf')
    db_mongo = utils.setupMongoDB()
    write.writeFromSDF(db_mongo, 'data/test_data/first_200.props.sdf', 'test')
    similarity.addMorganFingerprints(db_mongo)
    thresholds = [1, 0.8, 0.6, 0.4, 0.2]
    for i in range(200):
        mol = Chem.Mol(db_python[i]['rdmol'])
        last = []
        for t in thresholds:
            search_mongo = similarity.similaritySearch(mol, db_mongo, t)
            assert len(search_mongo) >= len(last)
            assert (all(l in search_mongo for l in last))
            last = search_mongo

