from .substructure import *
from .similarity import *


def PrepareForSearch(db, mol_collection, count_collection, perm_collection):
    print("Preparing database and collections for search...")
    AddPatternFingerprints(mol_collection)
    AddMorganFingerprints(mol_collection, count_collection)
    similarity.AddRandPermutations(perm_collection)
    similarity.AddLocalityHashes(mol_collection, perm_collection, DEFAULT_BUCKET_N)
    similarity.AddHashCollections(db, mol_collection)
    print("Added pattern fps, morgan fps, and support for LSH.")