from pathlib import Path
import sys
import os
import mongomock
from mongordkit.Database import create

sys.path.append(Path('.').resolve().parent.parent)

import mongordkit.Database

def test_registrationPresent():
    client = mongomock.MongoClient()
    db = client.db
    collection = db['registration']
    collection.insert_one(mongordkit.Database.utils.STANDARD_SETTING)
    assert 'registration' in db.list_collection_names()

def test_create():
    hdb = create.createFromHostPort('MyDatabase', host=None, port=None)
    adb = create.createFromURL('MyDatabase', url=None)



