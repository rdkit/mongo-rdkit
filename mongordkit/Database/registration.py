"""
Constructs document representations of molecules and
handles data registration settings.
"""
import rdkit
from bson import Binary
from rdkit import Chem
from rdkit.Chem import rdMolHash
import pickle


DEFAULT_SCHEME_NAME = 'default'
DEFAULT_AUTHOR = 'package-native'
DEFAULT_PREPROCESS = False
DEFAULT_INDEX = 'inchikey_standard'


RDKIT_HASH_FUNCTIONS = rdkit.Chem.rdMolHash.HashFunction.names
HASH_FUNCTIONS = {}
HASH_FUNCTIONS['MoleculeHashString'] = rdMolHash.GenerateMoleculeHashString
HASH_FUNCTIONS['inchi_standard'] = Chem.MolToInchi
HASH_FUNCTIONS['inchikey_standard'] = Chem.MolToInchiKey
HASH_FUNCTIONS['inchi_KET_15T'] = lambda rdmol: Chem.MolToInchi(rdmol, options='-KET -15T')
HASH_FUNCTIONS['inchikey_KET_15T'] = lambda rdmol: Chem.MolToInchiKey(rdmol, options='-KET -15T')
HASH_FUNCTIONS['noiso_smiles'] = lambda rdmol: Chem.MolToSmiles(rdmol, isomericSmiles=False)
HASH_FUNCTIONS['cx_smiles'] = Chem.MolToCXSmiles


class MolDocScheme():

    def __init__(self):
        self.scheme_name = DEFAULT_SCHEME_NAME
        self.author = DEFAULT_AUTHOR
        self.pre_processed = DEFAULT_PREPROCESS
        self.index_option = DEFAULT_INDEX
        self.rdkit_hashes = set(RDKIT_HASH_FUNCTIONS.keys())
        self.hashes = set(HASH_FUNCTIONS.keys())
        self.fingerprints = {}
        self.value_fields = {}

    def __repr__(self):
        return 'Molecule document representation schema. ' \
               '{Name: ' + str(self.scheme_name) + \
               ', author: ' + str(self.author) + \
               ', index' + str(self.index_option) + \
               '}'

    def set_index(self, new_index):
        if new_index not in HASH_FUNCTIONS.keys() and new_index not in RDKIT_HASH_FUNCTIONS.keys():
            raise Exception("Please add this hash first.")
        else:
            self.index_option = new_index

    def get_index_value(self, rdmol):
        if self.index_option in HASH_FUNCTIONS.keys():
            return HASH_FUNCTIONS[self.index_option](rdmol)
        elif self.index_option in RDKIT_HASH_FUNCTIONS.keys():
            return rdMolHash.MolHash(rdmol, RDKIT_HASH_FUNCTIONS[self.index_option])
        else:
            raise Exception("Specified index option does not exist.")

    def add_hash_field(self, field_name, field_method):
        self.hashes.add(field_name)
        HASH_FUNCTIONS[field_name] = field_method

    def add_value_field(self, field_name, field_value):
        self.value_fields[field_name] = field_value

    def remove_field(self, field_name):
        if field_name in self.hashes:
            self.hashes.remove(field_name)
            print(f'removed {field_name} from scheme')
        if field_name in self.value_fields.keys():
            self.value_fields.pop(field_name)
            print(f'removed {field_name} from scheme')
        if field_name in self.rdkit_hashes:
            self.rdkit_hashes.remove(field_name)
            print(f'removed {field_name} from scheme')

    def generate_mol_doc(self, rdmol):
        molDoc = {
            'rdmol': Binary(rdmol.ToBinary()),
            'index': self.get_index_value(rdmol),
            'smiles': Chem.MolToSmiles(rdmol),
            'scheme': self.scheme_name,
            'hashes': {hash_name: HASH_FUNCTIONS[hash_name](rdmol) for hash_name in self.hashes},
            'rdkit_hashes': {hash_name: rdMolHash.MolHash(rdmol, RDKIT_HASH_FUNCTIONS[hash_name])
                             for hash_name in self.rdkit_hashes},
            'fingerprints': {fp: fp_method(rdmol) for fp, fp_method in self.fingerprints.items()},
            'value_data': {field_name: value for field_name, value in self.value_fields.items()}
        }
        return molDoc

