from mongordkit.Database import write
from mongordkit.Search import substructure, similarity
from mongordkit.Search.tests import utils
from rdkit import Chem


def test_addPatternFingerprints():
    db = utils.setupMockDB()
    write.WriteFromSDF(db.molecules, 'data/test_data/first_200.props.sdf')
    substructure.AddPatternFingerprints(db.molecules)
    counter = 0
    assert db.molecules.count_documents({"fingerprints.pattern_fp": {"$exists": True}}) == 200


def test_SubSearchAccuracy():
    db_mock = utils.setupMockDB()
    write.WriteFromSDF(db_mock.molecules, 'data/test_data/first_200.props.sdf')
    substructure.AddPatternFingerprints(db_mock.molecules)
    db_python = utils.setupPythonDB('data/test_data/first_200.props.sdf')
    for i in range(200):
        moldoc = db_python[i]
        pattern = Chem.Mol(moldoc['rdmol'])
        results_python = utils.SubSearchPython(pattern, db_python)
        results = substructure.SubSearch(pattern, db_mock.molecules)
        assert sorted(results_python) == sorted(results)

