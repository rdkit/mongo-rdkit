import pymongo
from bson import Binary
from rdkit import Chem


def populateFromSDF(db, sdf, chunk_size=100):
    """Populates a database instance DB with the
    contents of an SDF file SDF in a collection called molecules.
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
    global molecules
    if 'molecules' not in db.collection_names():
        molecules = db.molecules
    chunk = []
    for rdmol in Chem.ForwardSDMolSupplier(sdf):
        if rdmol is None:
            continue
        document = {
            'smiles': Chem.MolToSmiles(rdmol, isomericSmiles=True),
            'chembl_id': rdmol.GetProp('chembl_id'),
            'rdmol': Binary(rdmol.ToBinary()),
        }
        chunk.append(document)
        if len(chunk) == chunk_size:
            molecules.insert_many(chunk)
            chunk = []
    print()
    return db.molecules.count()

