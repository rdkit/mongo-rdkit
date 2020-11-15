import pymongo, rdkit, math
from bson import ObjectId
from rdkit import Chem
from rdkit.Chem import AllChem
from mongordkit import Database
import numpy as np
import functools

# Default configurations for a variety of constants.
DEFAULT_THRESHOLD = 0.8

# Morgan fingerprints
DEFAULT_MORGAN_RADIUS = 2
DEFAULT_MORGAN_LEN = 2048

# LSH constants
DEFAULT_BIT_N = 2048
DEFAULT_BUCKET_N = 25
DEFAULT_PERM_LEN = 2048
DEFAULT_PERM_N = 100

# PyMongo configurations
DEFAULT_BATCH_SIZE = 100

def SimSearchNaive(mol, mol_collection, threshold=DEFAULT_THRESHOLD):
    """
    Searches MOL_COLLECTION for molecules with Tanimoto
    similarity to MOL greater than or equal to THRESHOLD.
    :param mol: An rdmol object.
    :param mol_collection: A MongoDB collection that meets requirements.
    :param threshold: Tanimoto threshold for similarity. Defaults to 0.8
    :return: A list of SMILES that fulfill threshold, along with their tanimoto scores.
    """
    results = []
    qfp = list(AllChem.GetMorganFingerprintAsBitVect(mol, DEFAULT_MORGAN_RADIUS,
                                                     nBits=DEFAULT_MORGAN_LEN).GetOnBits())
    for molDoc in mol_collection.find():
        mfp = list((AllChem.GetMorganFingerprintAsBitVect(Chem.Mol(molDoc['rdmol']), DEFAULT_MORGAN_RADIUS,
                                                          nBits=DEFAULT_MORGAN_LEN).GetOnBits()))
        tanimoto = calc_tanimoto(qfp, mfp)
        if tanimoto >= threshold:
            results.append(["tanimoto: " + str(tanimoto)] +
                           ["index" + ": " + str(mol["index"])])
    return results


def AddMorganFingerprints(mol_collection, count_collection):
    """
    Adds Morgan fingerprint bit vectors and counts with RADIUS and
    LENGTH bits to all documents in MOL_COLLECTION. Inserts bit frequency information
    into COUNT_COLLECTION.
    :param mol_collection: A MongoDB collection storing molecules.
    :param count_collection: A MongoDB collection storing frequency of bits.
    :param radius: Radius of the desired Morgan fingerprints. Defaults to 2.
    :param length: NBits in desired Morgan fingerprints. Defaults to 2048.
    """
    for m in mol_collection.find(batch_size=DEFAULT_BATCH_SIZE):
        bit_vector = list((AllChem.GetMorganFingerprintAsBitVect(Chem.Mol(m['rdmol']),
                                                                 DEFAULT_MORGAN_RADIUS,
                                                                 nBits=DEFAULT_MORGAN_LEN).GetOnBits()))
        count = len(bit_vector)
        mol_collection.update_one({'_id': m['_id']}, {'$set': {'fingerprints.morgan_fp': {'bits': bit_vector,
                                                                           'count': count}}})
    counts = {}
    for m in mol_collection.find(batch_size=DEFAULT_BATCH_SIZE):
        for bit in m['fingerprints']['morgan_fp']['bits']:
            counts[bit] = counts.get(bit, 0) + 1
    for k, v in counts.items():
        count_collection.insert_one({'_id': k, 'count': v})
    a = count_collection.find_one()
    mol_collection.create_index('fingerprints.morgan_fp.bits')
    mol_collection.create_index('fingerprints.morgan_fp.count')
    return None


def SimSearch(mol, mol_collection, count_collection=None, threshold=DEFAULT_THRESHOLD):
    """
    Searches MOL_COLLECTION for molecules with Tanimoto
    similarity to MOL greater than or equal to THRESHOLD.
    :param mol: An rdmol object.
    :param mol_collection: A MongoDB collection that meets requirements.
    :param threshold: Tanimoto threshold for similarity. Defaults to 0.8
    :param return_fields: A list of strings that indicate fields to return in addition
    to the Tanimoto threshold. Defaults to only canonical smiles.
    :return: A list of SMILES that fulfill threshold, along with their tanimoto scores.
    """
    results = []
    a = count_collection.find_one()
    qfp = list(AllChem.GetMorganFingerprintAsBitVect(mol, DEFAULT_MORGAN_RADIUS,
                                                     nBits=DEFAULT_MORGAN_LEN).GetOnBits())
    if threshold == 0:
        for mol in mol_collection.find():
            tanimoto = [str(calc_tanimoto(qfp, mol['fingerprints']['morgan_fp']['bits']))]
            results.append(tanimoto + [str(mol['index'])])
        return results
    qfp_count = len(qfp)
    fp_min = int(math.ceil((threshold * qfp_count)))
    try:
        fp_max = int(qfp_count / threshold)
    except ZeroDivisionError:
        fp_max = float('inf')
    req_common_count = qfp_count - fp_min + 1
    if count_collection:
        req_common_bits = [count['_id'] for count in count_collection.find(
            {'_id': {'$in': qfp}}).sort('count', 1).limit(req_common_count)]
    else:
        req_common_bits = qfp[:req_common_count]
    for mol in mol_collection.find({'fingerprints.morgan_fp.count': {'$gte': fp_min, '$lte': fp_max},
                                 'fingerprints.morgan_fp.bits': {'$in': req_common_bits}}):
        tanimoto = calc_tanimoto(qfp, mol['fingerprints']['morgan_fp']['bits'])
        if tanimoto >= threshold:
            results.append([tanimoto] +
                       [str(mol["index"])])
    return results


def SimSearchAggregate(mol, mol_collection, count_collection=None, threshold=DEFAULT_THRESHOLD):
    """
    Searches MOL_COLLECTION for molecules with Tanimoto
    similarity to MOL greater than or equal to THRESHOLD.
    This method uses a MongoDB aggregation pipeline.
    :param mol: An rdmol object.
    :param mol_collection: A MongoDB collection that meets requirements.
    :param threshold: Tanimoto threshold for similarity. Defaults to 0.8
    :return: A list of SMILES that fulfill threshold, along with their tanimoto scores.
    """
    results = []
    qfp = list(AllChem.GetMorganFingerprintAsBitVect(mol, DEFAULT_MORGAN_RADIUS,
                                                     nBits=DEFAULT_MORGAN_LEN).GetOnBits())
    if threshold == 0:
        for mol in mol_collection.find():
            tanimoto = [str(calc_tanimoto(qfp, mol['fingerprints']['morgan_fp']['bits']))]
            results.append(tanimoto + [str(mol['index'])])
        return results
    qfp_count = len(qfp)
    fp_min = int(math.ceil((threshold * qfp_count)))
    try:
        fp_max = int(qfp_count / threshold)
    except ZeroDivisionError:
        fp_max = float('inf')
    req_common_count = qfp_count - fp_min + 1
    if count_collection:
        req_common_bits = [count['_id'] for count in count_collection.find(
            {'_id': {'$in': qfp}}).sort('count', 1).limit(req_common_count)]
    else:
        req_common_bits = qfp[:req_common_count]
    aggregate = [
        {'$match': {'fingerprints.morgan_fp.count': {'$gte': fp_min, '$lte': fp_max},
                    'fingerprints.morgan_fp.bits': {'$in': req_common_bits}}},
        {'$project':
             {'tanimoto': {'$let': {'vars': {'common': {'$size': {'$setIntersection': ['$fingerprints.morgan_fp.bits', qfp]}}},
                                    'in': {'$divide': ['$$common', {'$add': [qfp_count, {'$subtract': ['$fingerprints.morgan_fp.count', '$$common']}]}]}
                                    }
                           },
              'index': 1
            }
         },
        {'$match': {'tanimoto': {'$gte': threshold}}}
    ]
    response = mol_collection.aggregate(aggregate)
    return [[r['tanimoto'], r['index']] for r in response]


def AddRandPermutations(perm_collection, len=DEFAULT_PERM_LEN, num=DEFAULT_PERM_N):
    """
    Uses the function get_permutations to generate NUM random permutations
    of bits of length LEN and saves each in COLLECTION as a separate document.
    :param collection: A MongoDB collection.
    :param len: Length of fingerprints to generate permutations for.
    :param num: Number of random permutations to save.
    :return: None
    """
    perm_collection.insert_many([{'_id':i, 'permutation': perm.tolist()}
                                     for i, perm in enumerate(get_permutations(len, num))])


def AddLocalityHashes(mol_collection, perm_collection, nBuckets=DEFAULT_BUCKET_N):
    """
    Adds locality-sensitive hash values to each document in MOL_COLLECTION
    based on permutations in PERM_COLLECTION. This method requires documents
    in PERM_COLLECTION to have a 'permutation' field and documents in
    MOL_COLLECTION to have a 'rdmol' field.
    :param mol_collection: A MongoDB collection.
    :param perm_collection: A MongoDB collection.
    :return: None
    """
    length = len(perm_collection.find_one()['permutation'])
    permutations = []
    for doc in perm_collection.find(batch_size=DEFAULT_BATCH_SIZE):
        if len(doc['permutation']) != length:
            Exception('Permutations must split evenly into buckets')
            return
        permutations.append(doc['permutation'])
    for moldoc in mol_collection.find(batch_size=DEFAULT_BATCH_SIZE):
        mol = Chem.Mol(moldoc['rdmol'])
        min_hash = get_min_hash(mol, permutations)
        hash_groups = hash_to_buckets(min_hash, nBuckets, length)
        mol_collection.update_one({'_id': moldoc['_id']}, {'$set': {'LSH': hash_groups}})


def AddHashCollections(db, mol_collection):
    """
    Creates different collections in DB that each store a subset of molecules
    in MOL_COLLECTION hashed through LSH.
    :param mol_collection: a MongoDB collection.
    :return: None
    """
    for moldoc in mol_collection.find(batch_size=DEFAULT_BATCH_SIZE):
        hash_groups = moldoc['LSH']
        for number, hash in enumerate(hash_groups):
            db['LSHash_' + str(number)].update_one({'_id': hash},
                                                   {'$push': {'molecules': moldoc['_id']}}, True)


def SimSearchLSH(mol, db, mol_collection, perm_collection, count_collection, threshold=DEFAULT_THRESHOLD):
    """
    Conducts a similarity search for query molecule MOL in MOL_COLLECTION
    with Tanimoto threshold THRESHOLD.
    TODO: this method should check whether or not the collection and database have been prepared properly.
    :param mol:
    :param db:
    :param threshold:
    :return:
    """
    results = []
    qfp = list(AllChem.GetMorganFingerprintAsBitVect(mol, DEFAULT_MORGAN_RADIUS,
                                                     nBits=DEFAULT_MORGAN_LEN).GetOnBits())
    qfp_bits = [int(n) for n in qfp]
    if threshold == 0:
        for mol in db.molecules.find():
            tanimoto = [str(calc_tanimoto(qfp, mol['fingerprints']['morgan_fp']['bits']))]
            results.append(tanimoto + [str(mol['index'])])
        return results
    qfp_count = len(qfp)
    fp_min = int(math.ceil((threshold * qfp_count)))
    try:
        fp_max = int(qfp_count / threshold)
    except ZeroDivisionError:
        fp_max = float('inf')
    req_common_count = qfp_count - fp_min + 1
    if count_collection:
        req_common_bits = [count['_id'] for count in count_collection.find(
            {'_id': {'$in': qfp}}).sort('count', 1).limit(req_common_count)]
    else:
        req_common_bits = qfp[:req_common_count]
    permutations = [p['permutation'] for p in perm_collection.find()]
    min_hash = get_min_hash(mol, permutations)
    hash_groups = hash_to_buckets(min_hash)
    nested_res = []
    cursors = [db['LSHash_' + str(i)].find({'_id': h}, {'molecules': 1}) for i, h in enumerate(hash_groups)]
    for c in cursors:
        cursor = list(c)
        if len(cursor) == 0:
            continue
        else:
            nested_res.append(cursor[0]['molecules'])
    hashed_ids = [ObjectId(x) for x in (set([str(item) for sublist in nested_res for item in sublist]))]
    aggregate = [
        {'$match': {'_id':{'$in': hashed_ids}, 'fingerprints.morgan_fp.count': {'$gte': fp_min, '$lte': fp_max},
                    'fingerprints.morgan_fp.bits': {'$in': req_common_bits}}},
        {'$project':
             {'tanimoto': {'$let': {'vars': {'common': {'$size': {'$setIntersection': ['$fingerprints.morgan_fp.bits', qfp]}}},
                                    'in': {'$divide': ['$$common', {
                                        '$add': [qfp_count, {'$subtract': ['$fingerprints.morgan_fp.count', '$$common']}]}]}
                                    }
                           },
              'index': 1
              }
         },
        {'$match': {'tanimoto': {'$gte': threshold}}}
    ]
    response = mol_collection.aggregate(aggregate)
    return [[r['tanimoto'], r['index']] for r in response]


def get_min_hash(mol, permutations):
    qfp = [int(n) for n in list(AllChem.GetMorganFingerprintAsBitVect(mol, DEFAULT_MORGAN_RADIUS,
                                                                      nBits=DEFAULT_MORGAN_LEN))]
    min_hash = []
    for perm in permutations:
        for idx, i in enumerate(perm):
            if qfp[i]:
                min_hash.append(idx)
                break
    return min_hash


def hash_to_buckets(min_hash, num_buckets=DEFAULT_BUCKET_N, nBits=DEFAULT_BIT_N):
    if len(min_hash) % num_buckets:
        print('Length of min_hash must be divisible by number of buckets.')
        return
    buckets = []
    hash_per_bucket = int(len(min_hash) / num_buckets)
    num_bits = (nBits-1).bit_length()
    for b in range(num_buckets):
        buckets.append(functools.reduce(lambda x,y: (x << num_bits) + y, min_hash[b:(b + hash_per_bucket)]))
    return buckets


def get_permutations(len_permutations=DEFAULT_PERM_LEN, num_permutations=DEFAULT_PERM_N):
    """Gets NUM_PERMUTATIONS random permutations of numbers of length LEN_PERMUTATIONS each."""
    return map(lambda _: np.random.permutation(len_permutations), range(num_permutations))


def calc_tanimoto(Na, Nb):
    """Calculates the Tanimoto similarity coefficient between two sets NA and NB."""
    Nab = len(set(Na).intersection((set(Nb))))
    return float(Nab) / (len(Na) + len(Nb) - Nab)
