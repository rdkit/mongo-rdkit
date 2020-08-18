import pymongo, pickle
from bson import Binary
from rdkit import Chem
import rdkit
from rdkit.Chem import rdMolHash
from rdkit.Chem import rdinchi
from .registration import MolDocScheme, HASH_FUNCTIONS


def WriteFromSDF(mol_collection, sdf, scheme=MolDocScheme(),
                 reg_collection=None, chunk_size=100, limit=None, warnings=False):
    """
    Writes the contents of SDF to MOL_COLLECTION and creates
    an index on the index specificed in SCHEME. Optional parameters:
    - customize document structure by specifying a SCHEME MolDocScheme object.
    - write the scheme object into a separate collection by specifying REG_COLLECTION.
    - customize how many molecules are inserted at a time by setting CHUNK_SIZE.
    - limit the number of molecules written in by setting LIMIT.
    :param mol_collection: A MongoDB collection.
    :param sdf: A Python File object of the desired SDF file.
    :param scheme: A registration.MolDocScheme() object.
    :param reg_collection: A MongoDB collection.
    :param chunk_size: Integer indicating how many molecules inserted at a time.
    :param limit: Integer indicating how many molecules to insert.
    """
    if not warnings:
        rdkit.RDLogger.DisableLog('rdApp.*')
    molecules = mol_collection
    print('populating mongodb collection with compounds from SDF...')
    chunk = []
    inserted = 0
    duplicates = 0
    if reg_collection:
        reg_collection.insert_one({scheme.scheme_name: pickle.dumps(scheme)})
    for rdmol in Chem.ForwardSDMolSupplier(sdf):
        if limit is not None and inserted >= limit:
            break
        if rdmol is None:
            continue
        index = scheme.get_index_value(rdmol)
        if mol_collection.count_documents({"index": index}) != 0:
            duplicates += 1
            continue
        document = scheme.generate_mol_doc(rdmol)
        chunk.append(document)
        inserted += 1
        if len(chunk) == chunk_size:
            molecules.insert_many(chunk)
            chunk = []
    if len(chunk) != 0:
        for i in chunk:
            if limit is not None and inserted >= limit:
                break
            else:
                molecules.insert_one(i)
                inserted += 1
    print("{} molecules successfully imported".format(inserted))
    print("{} duplicates skipped".format(duplicates))
    mol_collection.create_index('index')
    return inserted


def WriteFromMolList(mol_collection, list, scheme=MolDocScheme(),
                     reg_collection=None, chunk_size=100, limit=None):
    """
    Writes the contents of LIST to MOL_COLLECTION and creates
    an index on the index specificed in SCHEME. Optional parameters:
    - customize document structure by specifying a SCHEME MolDocScheme object.
    - write the scheme object into a separate collection by specifying REG_COLLECTION.
    - customize how many molecules are inserted at a time by setting CHUNK_SIZE.
    - limit the number of molecules written in by setting LIMIT.
    :param mol_collection: A MongoDB collection.
    :param list: A Python list of rdmol objects.
    :param scheme: A registration.MolDocScheme() object.
    :param reg_collection: A MongoDB collection.
    :param chunk_size: Integer indicating how many molecules inserted at a time.
    :param limit: Integer indicating how many molecules to insert.
    """
    molecules = mol_collection
    print('populating mongodb collection with compounds from list...')
    chunk = []
    inserted = 0
    duplicates = 0
    if reg_collection:
        reg_collection.insert_one({scheme.scheme_name: pickle.dumps(scheme)})
    for rdmol in list:
        if limit is not None and inserted >= limit:
            break
        if rdmol is None:
            continue
        index = scheme.get_index_value(rdmol)
        if mol_collection.count_documents({"index": index}) != 0:
            duplicates += 1
            continue
        document = scheme.generate_mol_doc(rdmol)
        chunk.append(document)
        inserted += 1
        if len(chunk) == chunk_size:
            molecules.insert_many(chunk)
            chunk = []
    if len(chunk) != 0:
        for i in chunk:
            if limit is not None and inserted >= limit:
                break
            else:
                molecules.insert_one(i)
                inserted += 1
    print("{} molecules successfully imported".format(inserted))
    print("{} duplicates skipped".format(duplicates))
    mol_collection.create_index('index')
    return inserted


def addSetting(db, setting):
    """
    Placeholder code for a function that writes a setting called
    SETTING to the collection 'registration' in the database DB.
    :param db: a MongoDB database instance.
    :param setting: a dictionary that with setting attributes.
    :return:
    """


def listAvailableSettings(db):
    """
    Placeholder code for a function that lists available
    registration settings in a database DB.
    :param db: a MongoDB database instance.
    :return: A list of available settings.
    """

### UTILITY FUNCTIONS
