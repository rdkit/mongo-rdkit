import pymongo
import rdkit
import math
from rdkit import Chem
from rdkit.Chem import AllChem
from mongordkit import Database


def substructureSearchNaive(pattern, db, chilarity=False):
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
        if rdmol.HasSubstructMatch(pattern, useChilarity=chilarity):
            results.append([molDoc['smiles']])
