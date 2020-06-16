from pathlib import Path
import sys
import os
import mongomock

sys.path.append(Path('.').resolve().parent.parent)

import mongordkit.Database

def test_registrationPresent():
    client = mongomock.MongoClient()
    db = client.db
    collection = db['registration']
    collection.insert_one(mongordkit.Database.utils.STANDARD_SETTING)
    assert 'registration' in db.list_collection_names()



