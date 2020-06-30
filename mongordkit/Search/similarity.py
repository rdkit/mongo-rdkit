import pymongo
import rdkit
import math
from rdkit import Chem
from rdkit.Chem import AllChem
from mongordkit import Database


def addMorganFingerprints(db, radius=2, length=1024):
    """Add Morgan fingerprints with radius RADIUS and
    LENGTH bits to all molecules in a database DB. Creates
    indices on bits and counts in order to optimize query speed.
    Also creates a count collection that stores how common
    each Morgan fingerprint is. This aids in similarity search.
    """
    for m in db.molecules.find():
        bit_vector = list((AllChem.GetMorganFingerprintAsBitVect(Chem.Mol(m['rdmol']), 2, nBits=1024).GetOnBits()))
        count = len(bit_vector)
        db.molecules.update_one({'_id': m['_id']}, {'$set': {'morgan_fp': {'bits': bit_vector,
                                                                           'count': count}}})
    counts = {}
    chunk_size = 100
    for m in db.molecules.find():
        for bit in m['morgan_fp']['bits']:
            counts[bit] = counts.get(bit, 0) + 1
    for k, v in counts.items():
        db.mfp_counts.insert_one({'_id': k, 'count': v})
    db.molecules.create_index('morgan_fp.bits')
    db.molecules.create_index('morgan_fp.count')
    return None


def similaritySearch(mol, db, threshold=0.8):
    """Perform a similarity search for an rdmol MOL
    against a MongoDB database DB with threshold THRESHOLD.
    THRESHOLD defaults to 0.8. Searching with a threshold of
    0 returns the entire database.
    """
    results = []
    qfp = list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024).GetOnBits())
    if threshold == 0:
        for mol in db.molecules.find():
            results.append([calc_tanimoto(qfp, mol['morgan_fp']['bits']), mol['smiles']])
        return results
    qfp_count = len(qfp)
    fp_min = int(math.ceil((threshold * qfp_count)))
    try:
        fp_max = int(qfp_count / threshold)
    except ZeroDivisionError:
        fp_max = float('inf')
    req_common_count = qfp_count - fp_min + 1
    if 'mfp' in db.list_collection_names():
        req_common_bits = [count['_id'] for count in db.mfp_counts.find(
            {'_id': {'$in': qfp}}).sort('count', 1).limit(req_common_count)]
    else:
        req_common_bits = qfp[:req_common_count]
    for mol in db.molecules.find({'morgan_fp.count': {'$gte': fp_min, '$lte': fp_max},
                                 'morgan_fp.bits': {'$in': req_common_bits}}):


        tanimoto = calc_tanimoto(qfp, mol['morgan_fp']['bits'])
        if tanimoto >= threshold:
            results.append([tanimoto, mol['smiles']])
    return results


def similaritySearchAggregate(mol, db, threshold=0.8):
    """Perform a similarity search for an rdmol MOL
    against a MongoDB database DB with threshold THRESHOLD.
    THRESHOLD defaults to 0.8. This method uses a MongoDB aggregation pipeline.
    Searching with a threshold of 0 returns the entire database.
    """
    results = []
    qfp = list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024).GetOnBits())
    if threshold == 0:
        for mol in db.molecules.find():
            results.append([calc_tanimoto(qfp, mol['morgan_fp']['bits']), mol['smiles']])
        return results
    qfp_count = len(qfp)
    fp_min = int(math.ceil((threshold * qfp_count)))
    try:
        fp_max = int(qfp_count / threshold)
    except ZeroDivisionError:
        fp_max = float('inf')
    req_common_count = qfp_count - fp_min + 1
    if 'mfp' in db.list_collection_names():
        req_common_bits = [count['_id'] for count in db.mfp_counts.find(
            {'_id': {'$in': qfp}}).sort('count', 1).limit(req_common_count)]
    else:
        req_common_bits = qfp[:req_common_count]
    aggregate = [
        {'$match': {'morgan_fp.count': {'$gte': fp_min, '$lte': fp_max},
                    'morgan_fp.bits': {'$in': req_common_bits}}},
        {'$project':
             {'tanimoto': {'$let': {'vars': {'common': {'$size': {'$setIntersection': ['$morgan_fp.bits', qfp]}}},
                                    'in': {'$divide': ['$$common', {'$add': [qfp_count, {'$subtract': ['$morgan_fp.count', '$$common']}]}]}
                                    }
                           },
              'smiles': 1
            }
         },
        {'$match': {'tanimoto': {'$gte': threshold}}}
    ]
    response = db.molecules.aggregate(aggregate)
    return [[r['tanimoto'], r['smiles']] for r in response]


def similaritySearchLSH(mol, db, threshold=0.8):
    return None


def calc_tanimoto(Na, Nb):
    Nab = len(set(Na).intersection((set(Nb))))
    return float(Nab) / (len(Na) + len(Nb) - Nab)
