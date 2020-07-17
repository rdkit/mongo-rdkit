import pymongo
from bson import Binary
from rdkit import Chem
from rdkit.Chem import rdMolHash
from rdkit.Chem import rdinchi
from .utils import *

def writeFromSDF(db, sdf, src, reg_option="standard_setting", index_option="inchikey", chunk_size=100, limit=None):
    """
    Writes to database instance DB the
    contents of an SDF file SDF from SRC in a collection called molecules.
    To limit the number of calls to the database,
    chunks the molecules into groups of CHUNK_SIZE,
    which defaults to 100. Returns the number of molecules
    inserted into the collection during its call.
    :param db: A mongoDB database instance.
    :param sdf: A Python File object of the desired SDF file.
    :param src: For the user's utility, a string that indicates where the file originates.
    :param reg_option: Allows control over the inserted document structure. Defaults to a standard
    setting that includes an index, SMILES, an rdkit Molecule, and a registration setting.
    :param index_option: Allows control over indexing settings. Defaults to generating an Inchikey for
    each molecule.
    :param chunk_size: How many documents are inserted into the database instance at a time.
    :return: The total number of molecules inserted into the collection.
    """
    molecules = db.molecules
    print('populating mongodb collection with compounds from chembl...')
    # This is placeholder code for when more registration options exist.
    if reg_option is not "standard_setting" and "registration" in db.list_collection_names():
        settings = db.registration.find_one({"$name": reg_option})

    if reg_option is not "standard_setting" and "registration" not in db.collection_names():
        print("The specified setting does not exist. Will only insert default molecules")

    if index_option not in VALID_IDS:
        options = ', '.join(VALID_IDS)
        raise ValueError("id_option must be one of {}".format(options))
    else:
        hash = HASH_FUNCTIONS[index_option]
    chunk = []
    inserted = 0
    for rdmol in Chem.ForwardSDMolSupplier(sdf):
        if rdmol is None:
            continue
        index = hash(rdmol)
        if db.molecules.count_documents({"index": index}) != 0:
            continue
        # Placeholder for where molecule standardization might take place.
        document = {
            'index': hash(rdmol),
            'smiles': Chem.MolToSmiles(rdmol),
            'rdmol': Binary(rdmol.ToBinary()),
            'registration_setting': reg_option
        }
        # Placeholder for adding setting specific fields to the document.
        chunk.append(document)
        if limit is not None and inserted > limit:
            break
        if len(chunk) == chunk_size:
            molecules.insert_many(chunk)
            inserted += chunk_size
            chunk = []
    if inserted + len(chunk) <= limit:
        molecules.insert_many(chunk)
    inserted += len(chunk)
    print("{} molecules successfully imported".format(inserted))
    return inserted

def WriteMolList(db, list, src, reg_option="standard_setting", index_option="inchikey", chunk_size=100, limit=None):
    """
    Writes to database instance DB a list of rdmols in the standard format into a molecules collection.
    :param db: A mongoDB database instance.
    :param sdf: A Python list of rdmol objects.
    :param src: For the user's utility, a string that indicates where the file originates.
    :param reg_option: Allows control over the inserted document structure. Defaults to a standard
    setting that includes an index, SMILES, an rdkit Molecule, and a registration setting.
    :param index_option: Allows control over indexing settings. Defaults to generating an Inchikey for
    each molecule.
    :param chunk_size: How many documents are inserted into the database instance at a time.
    :return: The total number of molecules inserted into the collection.
    """
    molecules = db.molecules
    print('populating mongodb collection with compounds from chembl...')
    # This is placeholder code for when more registration options exist.
    if reg_option is not "standard_setting" and "registration" in db.list_collection_names():
        settings = db.registration.find_one({"$name": reg_option})

    if reg_option is not "standard_setting" and "registration" not in db.collection_names():
        print("The specified setting does not exist. Will only insert default molecules")

    if index_option not in VALID_IDS:
        options = ', '.join(VALID_IDS)
        raise ValueError("id_option must be one of {}".format(options))
    else:
        hash = HASH_FUNCTIONS[index_option]
    chunk = []
    inserted = 0
    for rdmol in list:
        if rdmol is None:
            continue
        index = hash(rdmol)
        if db.molecules.count_documents({"index": index}) != 0:
            continue
        # Placeholder for where molecule standardization might take place.
        document = {
            'index': hash(rdmol),
            'smiles': Chem.MolToSmiles(rdmol),
            'rdmol': Binary(rdmol.ToBinary()),
            'registration_setting': reg_option
        }
        # Placeholder for adding setting specific fields to the document.
        chunk.append(document)
        if limit is not None and inserted > limit:
            break
        if len(chunk) == chunk_size:
            molecules.insert_many(chunk)
            inserted += chunk_size
            chunk = []
    if len(chunk) != 0:
        if limit is None or (inserted + len(chunk) <= limit):
            molecules.insert_many(chunk)
    inserted += len(chunk)
    print("{} molecules successfully imported".format(inserted))
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