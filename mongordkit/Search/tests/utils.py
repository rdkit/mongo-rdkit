from bson import Binary
from mongomock import ObjectId
import sys
from pathlib import Path
import mongordkit
import pymongo
import rdkit
from rdkit import Chem
import mongomock
from rdkit.Chem import AllChem
from mongordkit.Database import write, utils
from mongordkit.Search import similarity


def similaritySearchPython(query_mol, molecules, threshold):
    """
    Takes in a list of molecules and performs a brute force
    similarity search on them with query_mol as query.
    :param query_mol: The query_molecule, as an rdkit mol.
    :param molecules: A list of molecules, represented as dictionaries.
    :return: The resulting similar molecules.
    """
    results = []
    qfp = list(AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=1024).GetOnBits())
    for mol in molecules:
        mfp = list(AllChem.GetMorganFingerprintAsBitVect(Chem.Mol(mol['rdmol']), 2, nBits=1024).GetOnBits())
        tanimoto = calc_tanimoto(qfp, mfp)
        if calc_tanimoto(qfp, mfp) >= threshold:
            results.append([tanimoto, mol['smiles']])
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


def setupMongoDB():
    """
    WARNING: THIS TEST DIRECTLY MODIFIES YOUR LOCAL MONGODB INSTANCE.
    """
    client = pymongo.MongoClient()
    db = client.db
    db.molecules.drop()
    db.mfp_counts.drop()
    return client.db


def calc_tanimoto(Na, Nb):
    Nab = len(set(Na).intersection((set(Nb))))
    return float(Nab) / (len(Na) + len(Nb) - Nab)