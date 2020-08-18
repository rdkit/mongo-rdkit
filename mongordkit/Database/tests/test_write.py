from pathlib import Path
import sys
import os
import pytest
import mongomock
from rdkit import Chem

sys.path.append(Path('.').resolve().parent.parent)

from mongordkit.Database import write
from mongordkit.Database import registration

def setupDB():
    client = mongomock.MongoClient()
    return client.db

class TestWrite:

    def test_writeCount(self):
        data_scheme = registration.MolDocScheme()
        db = setupDB()
        assert 200 == write.WriteFromSDF(db.molecules, 'data/test_data/first_200.props.sdf', data_scheme)

    def test_invalidIndex(self):
        db = setupDB()
        data_scheme = registration.MolDocScheme()
        with pytest.raises(Exception):
            data_scheme.set_index('moo')

    def test_hashes(self):
        db = setupDB()
        data_scheme = registration.MolDocScheme()
        data_scheme.set_index("CanonicalSmiles")
        assert 200 == write.WriteFromSDF(db.molecules, 'data/test_data/first_200.props.sdf', data_scheme)
        data_scheme.set_index("inchi_standard")
        assert 200 == write.WriteFromSDF(db.molecules, 'data/test_data/first_200.props.sdf', data_scheme)

    def test_uniqueInsertion(self):
        db = setupDB()
        data_scheme = registration.MolDocScheme()
        assert 200 == write.WriteFromSDF(db.molecules, 'data/test_data/first_200.props.sdf', data_scheme)
        assert 0 == write.WriteFromSDF(db.molecules, 'data/test_data/first_200.props.sdf', data_scheme)

    def test_writeLimit(self):
        db = setupDB()
        data_scheme = registration.MolDocScheme()
        assert 10 == write.WriteFromSDF(db.molecules, 'data/test_data/first_200.props.sdf', data_scheme, limit=10)

    def test_WriteMolListCount(self):
        db = setupDB()
        data_scheme = registration.MolDocScheme()
        f = open('data/zinc.frags.500.q.smi')
        frags = [Chem.MolFromSmiles(line.split()[0]) for line in f]
        f.close()
        write.WriteFromMolList(db.molecules, frags, data_scheme)
        assert 499 == db.molecules.count_documents({})