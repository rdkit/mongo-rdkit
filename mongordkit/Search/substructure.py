import pymongo
import rdkit
import math
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolops import PatternFingerprint
from mongordkit import Database


def SubSearchNaive(pattern, db, chirality=False):
    """
    Super naive implementation of substructure search, without
    any fingerprint screening whatsoever.
    :param pattern: An rdmol object that represents a desired substructure pattern.
    :param db: The database over which to search. Must have a molecules collection.
    :param chilarity: Whether or not to include chilarity in search. Defaults to False.
    :return: A list of SMILES that have the desired substructure pattern.
    """
    results = []
    for molDoc in db.molecules.find():
        rdmol = Chem.Mol(molDoc['rdmol'])
        if rdmol.HasSubstructMatch(pattern, useChirality=chirality):
            results.append(molDoc['smiles'])
    return results


def AddPatternFingerprints(db, length=2048):
    """
    Adds pattern fingerprints to each molecule-representing document in DB.molecules.
    :param db: A MongoDB database instance that contains a molecules collection.
    :param length: The length of the uncompressed Pattern Fingerprint inserted.
    :return: None
    """
    for moldoc in db.molecules.find():
        mol = Chem.Mol(moldoc['rdmol'])
        bit_vector = list(PatternFingerprint(mol, length).GetOnBits())
        count = len(bit_vector)
        db.molecules.update_one({'_id': moldoc['_id']}, {'$set': {'pattern_fp': {'bits': bit_vector,
                                                                                 'count': count}}})
    return


def SubSearch(pattern, db, chirality=False):
    results = []
    query_fp = list(PatternFingerprint(pattern).GetOnBits())
    qfp_len = len(query_fp)
    for molDoc in db.molecules.find({'pattern_fp.count': {'$gte': qfp_len},
                                     'pattern_fp.bits': {'$all': query_fp}
                                     }):
        rdmol = Chem.Mol(molDoc['rdmol'])
        if rdmol.HasSubstructMatch(pattern, useChirality=chirality):
            results.append(molDoc['smiles'])
    return results


def SubSearchAggregate(pattern, db, chirality=False):
    return
    