from pathlib import Path
import sys
import os
import pytest
import mongomock

sys.path.append(Path('.').resolve().parent.parent)

from mongordkit.Database import write

def setupDB():
    client = mongomock.MongoClient()
    return client.db

def test_writeCount():
    db = setupDB()
    assert 200 == write.writeFromSDF(db, 'data/test_data/first_200.props.sdf', 'test')

def test_invalidIndex():
    with pytest.raises(ValueError):
        db = setupDB()
        write.writeFromSDF(db, 'data/test_data/first_200.props.sdf', 'test', 'standard_setting', 'canonica_smiles')

def test_hashes():
    db = setupDB()
    assert 200 == write.writeFromSDF(db, 'data/test_data/first_200.props.sdf', 'test', 'standard_setting', 'canonical_smiles')
    db = setupDB()
    assert 200 == write.writeFromSDF(db, 'data/test_data/first_200.props.sdf', 'test', 'standard_setting', 'het_atom_tautomer')

def test_uniqueInsertion():
    db = setupDB()
    write.writeFromSDF(db, 'data/test_data/first_200.props.sdf', 'test')
    assert 0 == write.writeFromSDF(db, 'data/test_data/first_200.props.sdf', 'test')
    assert 200 == write.writeFromSDF(db, 'data/test_data/first_200.props.sdf', 'test', 'standard_setting', 'canonical_smiles')