import pymongo
from bson import Binary
from rdkit import Chem
from rdkit.Chem import rdMolHash
from .utils import RegistrationUtils

def addSetting(db, setting):
    """
    Writes a setting called SETTING to the collection 'registration' in
    the database DB.
    :param db: a MongoDB database instance.
    :param setting: a dictionary that with setting attributes.
    :return:
    """

def listAvailableSettings(db):
    """
    Lists available write settings in a database DB.
    :param db: a MongoDB database instance.
    :return: A list of available settings.
    """


def writeFromSDF(db, sdf, src, reg_option='standard_setting', chunk_size=100):
    """Writes to database instance DB the
    contents of an SDF file SDF from SRC in a collection called molecules.
    To limit the number of calls to the database,
    chunks the molecules into groups of CHUNK_SIZE,
    which defaults to 100. Returns the number of molecules
    now in the collection.

    TODO:
    Add standardization process:
    Strip salt option
    RDKit standardization of mol
    Add optional representations of molecules from list
    of other non-default registration options

    Add fields to the document:
    Internal index (use RDKit Mol Hash)
    Source information
    Generalize chembl_id as 'external index'
    Standardized rdmol
    """
    if 'molecules' not in db.collection_names():
        molecules = db.molecules
    if reg_option is not 'standard_setting' and 'registration' not in db.collection_names():
        print('The specified setting does not exist. Will only insert default molecules')
    chunk = []
    #we can hard code a bunch of premade settings into the RegistrationSettings (thinking of rename) class.



    for rdmol in Chem.ForwardSDMolSupplier(sdf):
        if rdmol is None:
            continue
        #here we will insert the standard molecule. What should our specifications be?
        document = {
            'internal_id': rdMolHash.MolHash(rdmol), #what is the default hash function here?
            'smiles': Chem.MolToSmiles(rdmol, isomericSmiles=True),
            'src_id': rdmol.GetProp('chembl_id'),
            'rdmol': Binary(rdmol.ToBinary()),
            'registration_setting': reg_option
        }
        settings = db.registration.findOne({"$name": reg_option})
        if reg_option is not 'standard_setting' and settings:
            #here we will execute statements corresponding to a particular setting.
            if settings['salt_stripped']:
                #execute salt stripping
                print('salt stripped')
            print('foo')
        chunk.append(document)
        if len(chunk) == chunk_size:
            molecules.insert_many(chunk)
            chunk = []
    molecules.insert_many(chunk)
    return db.molecules.count()

