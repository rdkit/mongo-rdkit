import sys
import mongomock
import pytest
from pathlib import Path
from rdkit import Chem
from mongordkit.Database import write
from mongordkit.Search import similarity
from mongordkit.Search.tests import utils
import time

sys.path.append(Path('.').resolve().parent.parent.parent)


def test_pythonSetup():
    assert 200 == len(utils.setupPythonDB('data/test_data/first_200.props.sdf'))


def test_zeroThreshold():
    """
    If we search the database with a Tanimoto threshold of zero, we should get
    the entire contents of the database back.
    """
    db_python = utils.setupPythonDB('data/test_data/first_200.props.sdf')
    db_mongo = utils.setupMockDB()
    write.WriteFromSDF(db_mongo.molecules, 'data/test_data/first_200.props.sdf')
    similarity.AddMorganFingerprints(db_mongo.molecules, db_mongo.mfp_counts)

    mol = Chem.Mol(db_python[0]['rdmol'])
    assert 200 == len(utils.similaritySearchPython(mol, db_python, 0))
    assert 200 == len(similarity.SimSearchAggregate(mol, db_mongo.molecules, db_mongo.mfp_counts, 0))
    assert 200 == len(similarity.SimSearch(mol, db_mongo.molecules, db_mongo.mfp_counts, 0))


def test_similarityAccuracy():
    """
    Tests for basic accuracy against a brute-force constructed Python 'database'
    at thresholds 0.2, 0.4, 0.6, 0.8, and 1. This test is implemented using MongoMock.
    """
    db_python = utils.setupPythonDB('data/test_data/first_200.props.sdf')
    db_mongo = utils.setupMockDB()
    write.WriteFromSDF(db_mongo.molecules, 'data/test_data/first_200.props.sdf')
    similarity.AddMorganFingerprints(db_mongo.molecules, db_mongo.mfp_counts)
    thresholds = [0.2, 0.4, 0.6, 0.8, 1]
    for t in thresholds:
        for i in range(200):
            mol = Chem.Mol(db_python[i]['rdmol'])
            search_python = utils.similaritySearchPython(mol, db_python, t)
            search_mongo = similarity.SimSearch(mol, db_mongo.molecules, db_mongo.mfp_counts, t)
            assert sorted(search_python) == sorted(search_mongo)


@pytest.mark.mongo
def test_similarityAccuracyAggregate(mongoURI):
    """
    Tests for basic accuracy against a brute-force constructed Python 'database'
    at thresholds 0.2, 0.4, 0.6, 0.8, and 1. This test is relatively long and
    will modify your local MongoDB instance.
    """
    db_python = utils.setupPythonDB('data/test_data/first_200.props.sdf')
    if mongoURI == 'local':
        db_mongo = utils.setupMongoDB()
    else:
        db_mongo = utils.setupMongoDB(mongoURI)
    write.WriteFromSDF(db_mongo.molecules, 'data/test_data/first_200.props.sdf')
    similarity.AddMorganFingerprints(db_mongo.molecules, db_mongo.mfp_counts)
    thresholds = [0.2, 0.4, 0.6, 0.8, 1]
    counter = 0
    for t in thresholds:
        for i in range(200):
            mol = Chem.Mol(db_python[i]['rdmol'])
            search_python = utils.similaritySearchPython(mol, db_python, t)
            search_mongo_aggregate = similarity.SimSearchAggregate(mol, db_mongo.molecules, db_mongo.mfp_counts, t)
            assert sorted(search_python) == sorted(search_mongo_aggregate)
            print(counter)
            counter += 1


def test_similarityProgression():
    """
    Tests that decreasing similarity thresholds return increasing result lists.
    This test is implemented using MongoMock.
    """
    db_python = utils.setupPythonDB('data/test_data/first_200.props.sdf')
    db_mongo = utils.setupMockDB()
    write.WriteFromSDF(db_mongo.molecules, 'data/test_data/first_200.props.sdf')
    similarity.AddMorganFingerprints(db_mongo.molecules, db_mongo.mfp_counts)
    thresholds = [1, 0.8, 0.6, 0.4, 0.2]
    for i in range(200):
        mol = Chem.Mol(db_python[i]['rdmol'])
        last = []
        for t in thresholds:
            search_mongo = similarity.SimSearch(mol, db_mongo.molecules, db_mongo.mfp_counts, t)
            assert len(search_mongo) >= len(last)
            assert (all(l in search_mongo for l in last))
            last = search_mongo


@pytest.mark.mongo
def test_similarityAggregateProgression(mongoURI):
    """
    Tests that decreasing similarity thresholds return increasing result lists. This
    test will modify your local MongoDB instance.
    """
    db_python = utils.setupPythonDB('data/test_data/first_200.props.sdf')
    if mongoURI == 'local':
        db_mongo = utils.setupMongoDB()
    else:
        db_mongo = utils.setupMongoDB(mongoURI)
    write.WriteFromSDF(db_mongo.molecules, 'data/test_data/first_200.props.sdf')
    similarity.AddMorganFingerprints(db_mongo.molecules, db_mongo.mfp_counts)
    thresholds = [1, 0.8, 0.6, 0.4, 0.2]
    for i in range(200):
        mol = Chem.Mol(db_python[i]['rdmol'])
        last = []
        for t in thresholds:
            search_mongo = similarity.SimSearchAggregate(mol, db_mongo.molecules, db_mongo.mfp_counts, t)
            assert len(search_mongo) >= len(last)
            assert (all(l in search_mongo for l in last))
            last = search_mongo


@pytest.mark.mongo
def test_similarity_accuracy_LSH(mongoURI):
    db_python = utils.setupPythonDB('data/test_data/first_200.props.sdf')
    if mongoURI == 'local':
        db_mongo = utils.setupMongoDB()
    else:
        db_mongo = utils.setupMongoDB(mongoURI)
    write.WriteFromSDF(db_mongo.molecules, 'data/test_data/first_200.props.sdf')
    similarity.AddMorganFingerprints(db_mongo.molecules, db_mongo.mfp_counts)
    similarity.AddRandPermutations(db_mongo.permutations)
    similarity.AddLocalityHashes(db_mongo.molecules, db_mongo.permutations, 25)
    similarity.AddHashCollections(db_mongo, db_mongo.molecules)
    thresholds = [1, 0.8, 0.6, 0.4, 0.2]
    counter = 0
    for t in thresholds:
        for i in range(200):
            mol = Chem.Mol(db_python[i]['rdmol'])
            search_python = [result[1] for result in utils.similaritySearchPython(mol, db_python, t)]
            search_mongo_LSH = [result[1] for result in
                                similarity.SimSearchLSH(mol, db_mongo, db_mongo.molecules, db_mongo.permutations, t)]
            assert set(search_mongo_LSH).issubset(search_python)
            print(counter)
            counter += 1