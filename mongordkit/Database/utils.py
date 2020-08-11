import rdkit
import rdkit.Chem
from rdkit.Chem import rdinchi
from rdkit.Chem import rdMolHash


def inchiKeyHash(m):
    return rdinchi.MolToInchiKey(m)

def tautomerHash(m):
    return rdMolHash.MolHash(m, rdMolHash.HashFunction.HetAtomTautomer)

def canonicalHash(m):
    return rdMolHash.MolHash(m, rdMolHash.HashFunction.CanonicalSmiles)

STANDARD_SETTING = {
        'name': 'STANDARD',
        'salt_stripped': True,
        'charge_status': 'neutral',
        'fragments': 'removed',
    }

VALID_HASHES = {'inchikey', 'canonical_smiles', 'het_atom_tautomer'}

HASH_FUNCTIONS = {'inchikey': rdinchi.MolToInchiKey,
                      'canonical_smiles': canonicalHash,
                      'het_atom_tautomer': tautomerHash}