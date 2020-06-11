from pathlib import Path
import sys
import os
import mongomock

sys.path.append(Path('.').resolve().parent.parent)

import Database.utils

def test_registrationPresent():
    client = mongomock.MongoClient()
    db = client.db
    collection = db['registration']
    collection.insert_one(Database.utils.RegistrationUtils.standard_setting)
    assert 'registration' in db.list_collection_names()



