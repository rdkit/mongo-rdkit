from bson import Binary
from mongomock import ObjectId
import sys
from pathlib import Path
import mongordkit
import pymongo
import pymongo.errors
import rdkit
from rdkit import Chem
import mongomock
from rdkit.Chem import AllChem
from mongordkit.Database import write, utils
from mongordkit.Search import similarity
from rdkit.Chem.rdmolops import PatternFingerprint


def similaritySearchPython(query_mol, molecules, threshold):
    """
    Takes in a list of molecules and performs a brute force
    similarity search on them with query_mol as query.
    """
    results = []
    qfp = list(AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048).GetOnBits())
    for mol in molecules:
        mfp = list(AllChem.GetMorganFingerprintAsBitVect(Chem.Mol(mol['rdmol']), 2, nBits=2048).GetOnBits())
        tanimoto = calc_tanimoto(qfp, mfp)
        if calc_tanimoto(qfp, mfp) >= threshold:
            results.append([tanimoto, str(mol['index'])])
    return results


def SubSearchPython(pattern, molecules):
    """
    Takes in a list of molecules MOLECULES and performs a brute force
    substructure search on them with pattern PATTERN as query.
    """
    results = []
    for moldoc in molecules:
        mol = Chem.Mol(moldoc['rdmol'])
        if mol.HasSubstructMatch(pattern):
            results.append(moldoc['index'])
    return results


def setupPythonDB(sdf):
    """
    Inputs documents from SDF file into a Python set.
    :param sdf: path to an SDF file.
    :return: The resulting set.
    """
    data = []
    number_added = 0
    for rdmol in Chem.ForwardSDMolSupplier(sdf):
        if rdmol is None:
            continue
        hash = utils.HASH_FUNCTIONS['inchikey']
        document = {
            'index': hash(rdmol),
            'smiles': Chem.MolToSmiles(rdmol),
            'rdmol': Binary(rdmol.ToBinary()),
            'registration_setting': 'standard_setting'
        }
        data.append(document)
        number_added += 1
    print("{} molecules successfully imported".format(number_added))
    return data


def setupMongoDB(mongoURI=None):
    """
    WARNING: THIS DIRECTLY MODIFIES YOUR LOCAL MONGODB INSTANCE.
    """
    client = pymongo.MongoClient(mongoURI)
    client.drop_database('pytest_db')
    db = client['pytest_db']
    return db


def checkMongoDB(mongoURI=None):
    try:
        pymongo.MongoClient(mongoURI, serverSelectionTimeoutMS=50).server_info()
        return True
    except pymongo.errors.ServerSelectionTimeoutError:
        return False


def setupMockDB():
    """
    Set up a mock database for testing using MongoMock.
    """
    client = mongomock.MongoClient()
    db = client.db
    db.molecules.drop()
    db.mfp_counts.drop()
    return client.db


def calc_tanimoto(Na, Nb):
    Nab = len(set(Na).intersection((set(Nb))))
    return float(Nab) / (len(Na) + len(Nb) - Nab)